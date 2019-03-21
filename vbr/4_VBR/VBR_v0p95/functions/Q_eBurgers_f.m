function[VBR]=Q_eBurgers_f(VBR)
%% ========================================================
%% Initialization
%% ========================================================

% Read in State Variables
   f_vec = VBR.in.SV.f ;
   phi =  VBR.in.SV.phi ;
   %Mu = VBR.out.elastic.poro_Takei.Gu ; % Pa unrelaxed shear modulus
   if isfield(VBR.in.elastic,'poro_Takei')
     Mu = VBR.out.elastic.poro_Takei.Gu ;
   elseif isfield(VBR.in.elastic,'anharmonic')
     Mu = VBR.out.elastic.anharmonic.Gu ;
   end

   T_K_mat = VBR.in.SV.T_K ; % temperature [K]
   P_Pa_mat = VBR.in.SV.P_GPa.*1e9 ; % convert pressure GPa to Pa = GPa*1e9
   rho_mat = VBR.in.SV.rho ; % density
   d_mat = VBR.in.SV.dg_um ; % microns grain size

   w_vec = 2*pi.*f_vec ;
   Ju_mat = 1./Mu ;

% Settings
   method=0 ; % integration method (trapezoid, 0; quadrature, 1)

% Allocation
%  frequency is added as a new dimension at end of array.
%  e.g., if size(T_K_mat)=n1,n2 then the frequency dependent variables will
%        have size(Q)=n1,n2,nfreq.
   nfreq = numel(f_vec);
   Jz=zeros(size(Ju_mat));
   J1 = proc_add_freq_indeces(Jz,nfreq);
   J2 = J1; Q = J1; Qinv = J1; M = J1; V = J1;
   Vave = Jz;

% Read in reference values (Temp, Pressure and Grain Size) and parameters
  Burger_params=VBR.in.anelastic.eBurgers;
  bType=Burger_params.eBurgerMethod;
  TR = Burger_params.(bType).TR ;% Kelvins
  PR = Burger_params.(bType).PR *1e9; % convert pressure GPa to Pa = GPa*1e9
  dR = Burger_params.(bType).dR ; % microns grain size
  E = Burger_params.(bType).E ; % activation energy J/mol
  R = Burger_params.R ; % gas constant
  Vstar = Burger_params.(bType).Vstar ; % m^3/mol Activation Volume
  m_a = Burger_params.(bType).m_a ; % grain size exponent (anelastic)
  m_v = Burger_params.(bType).m_v ; % grain size exponent (viscous)
  alf = Burger_params.(bType).alf ;
  Delta = Burger_params.(bType).DeltaB ; % relaxation strength.
  Tau_LR = Burger_params.(bType).Tau_LR ;
  Tau_HR = Burger_params.(bType).Tau_HR ;
  Tau_MR = Burger_params.(bType).Tau_MR ;
  JF10G_UR = Burger_params.(bType).G_UR ;

  % peak parameters (will be 0 if eBurgType=='bg_only')
  DeltaP=Burger_params.(bType).DeltaP;
  sigma=Burger_params.(bType).sigma;
  Tau_PR=Burger_params.(bType).Tau_PR;

% Melt Effects
% sharper response than creep at phi_c, but less sensitive at higher phi (for HTB only)
  alpha = Burger_params.melt_alpha ; % post-critical melt fraction dependence
  phi_c = Burger_params.phi_c ; % critical melt fraction
  x_phi_c = Burger_params.x_phi_c ;% melt enhancement factor

% Calculate Scaling Matrix:
scale.visc=((d_mat./dR).^m_v).*exp((E/R).*(1./T_K_mat-1/TR)).*exp((Vstar/R).*(P_Pa_mat./T_K_mat-PR/TR));
scale.LHP=((d_mat./dR).^m_a).*exp((E/R).*(1./T_K_mat-1/TR)).*exp((Vstar/R).*(P_Pa_mat./T_K_mat-PR/TR));

% account for lack of truly melt free samples and the drop at the onset of melting.
if VBR.in.GlobalSettings.melt_enhacement==0
  x_phi_c=1;
else
  scale.visc = scale.visc .* x_phi_c ;
  scale.LHP = scale.LHP .* x_phi_c ;
end

% add melt effects
[scale_mat_prime] = sr_melt_enhancement(phi,alpha,x_phi_c,phi_c) ;
scale.visc = scale.visc ./ scale_mat_prime ;
scale.LHP = scale.LHP ./ scale_mat_prime ;

% ====================================================
% LOOP over the spatial DOMAIN of the state variables
% ====================================================

% use linear indexing!! will loop over n-dimensions of Ju_mat.
n_th = numel(Ju_mat); % number of thermodynamic states
for x1 = 1:n_th; % loop using linear index!

%   pull out variables at current index
    scale_visc = scale.visc(x1);
    scale_LHP = scale.LHP(x1);
    Ju = Ju_mat(x1) ;
    rho = rho_mat(x1) ;

    % tau limits & integration vector
    Tau_L = Tau_LR.*scale_LHP ;
    Tau_H = Tau_HR.*scale_LHP ;
    Tau_P = Tau_PR.*scale_LHP ;

    ntau = 500 ;
    Tau_X_vec = logspace(log10(Tau_L),log10(Tau_H),ntau) ;

    % maxwell relaxation time (period)
    % ??? adjust the constant to go from JF10 Gu reference to the current Gu.
    % ??? Tau_MR_adjusted = Tau_MR * JF10G_UR / Mu(x1);
    Tau_M = Tau_MR.*scale_visc ;

%   loop over frequency
    for i=1:nfreq
        i_glob = x1 + (i - 1) * n_th; % the linear index of the arrays with
                                         % a frequency index
        w = w_vec(i) ;

        if method==0
        %% MY METHOD--
            D_vec = (alf.*Tau_X_vec.^(alf-1))./(Tau_H^alf - Tau_L^alf) ;

            int_J1 = trapz(Tau_X_vec,(D_vec./(1+w^2.*Tau_X_vec.^2))) ;
            J1(i_glob) = (1+Delta.*int_J1) ;

            int_J2 = trapz(Tau_X_vec,((Tau_X_vec.*D_vec)./(1+w^2.*Tau_X_vec.^2))) ;
            J2(i_glob) = (w*Delta*int_J2 + 1/(w*Tau_M)) ;

        elseif method==1
        %% GEOFF's METHOD-- works !
        % GA: try iterative quadrature
        % (still seems to be a problem with high freq, w*tauH > 1 say)
            Tau_fac = alf.*Delta./(Tau_H.^alf - Tau_L.^alf);

            FINT1 = @(x) (x.^(alf-1))./(1+(w.*x).^2);
            int1 = Tau_fac.*quadl(FINT1, Tau_L, Tau_H);

            FINT2 = @(x) (x.^alf)./(1+(w.*x).^2);
            int2 = w.*Tau_fac.*quadl(FINT2, Tau_L, Tau_H);

            J1(i_glob) = (1 + int1);
            J2(i_glob) = (int2 + 1./(w.*Tau_M));
        end

        % add on peak if it's being used. May trigger warning, this integral
        % is not easy. 
        if DeltaP>0
          FINT2 = @(x) (exp(-(log(x./Tau_P)/sigma).^2/2)./(1+(w.*x).^2));
          int2a = quadgk(FINT2, 0, inf);
          J2(i_glob)=J2(i_glob)+DeltaP*w*(int2a)/(sigma*sqrt(2*pi));

          FINT1 = @(x) ( 1./x .* exp(-(log(x./Tau_P)/sigma).^2/2)./(1+(w.*x).^2));
          int1 = quadgk(FINT1, 0, inf);
          J1(i_glob)=J1(i_glob)+DeltaP*int1 / (sigma*sqrt(2*pi)) ;
        end

        % multiply on the unrelaxed compliance
        J1(i_glob)=Ju.*J1(i_glob);
        J2(i_glob)=Ju.*J2(i_glob);

        % See McCarthy et al, 2011, Appendix B, Eqns B6 !
        J2_J1_frac=(1+sqrt(1+(J2(i_glob)./J1(i_glob)).^2))/2;
        Qinv(i_glob) = J2(i_glob)./J1(i_glob).*(J2_J1_frac.^-1);
        Q(i_glob) = 1./Qinv(i_glob);

        %Q(i_glob) = J1(i_glob)./J2(i_glob) ;
        %Qinv(i_glob) = 1./Q(i_glob) ; % J2 / J1

        M(i_glob) = (J1(i_glob).^2 + J2(i_glob).^2).^(-0.5) ;
        V(i_glob) = sqrt(M(i_glob)./rho) ;

        Vave(x1) = Vave(x1) + V(i_glob); % add them all, divide by nfreq later
    % end loop over frequency
    end
% end the loop(s) over spatial dimension(s)
end

%% WRITE VBR
   VBR.out.anelastic.eBurgers.J1 = J1;
   VBR.out.anelastic.eBurgers.J2 = J2;
   VBR.out.anelastic.eBurgers.Q = Q;
   VBR.out.anelastic.eBurgers.Qinv = Qinv;
   VBR.out.anelastic.eBurgers.M=M;
   VBR.out.anelastic.eBurgers.V=V;
   VBR.out.anelastic.eBurgers.Vave = Vave./nfreq;

end
