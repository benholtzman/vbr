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
  TR = Burger_params.TR ;% Kelvins
  PR = Burger_params.PR *1e9; % convert pressure GPa to Pa = GPa*1e9
  dR = Burger_params.dR ; % microns grain size

  E = Burger_params.E ; % J/mol
  R = Burger_params.R ;
  Vstar = Burger_params.Vstar ; % m^3/mol (Activation Volume? or molar volume?)
  m = Burger_params.m ;

% Jackson n Faul 2010, table 1 :
  alf = Burger_params.alf ; % is this the same as n in Andrade ?
  Delta = Burger_params.Delta ;%1.4 ;% ; % relaxation strength..
  Tau_LR = Burger_params.Tau_LR ;
  Tau_HR = Burger_params.Tau_HR ;
  Tau_MR = Burger_params.Tau_MR ;

% EFFECT OF MELT, hypothetical:
% sharper response than creep at phi_c, but less sensitive at higher phi (for HTB only)
  alpha = Burger_params.melt_alpha ;
  phi_c = Burger_params.phi_c ;
  x_phi_c = Burger_params.x_phi_c ;


% ====================================================
% LOOP over the spatial DOMAIN of the state variables
% ====================================================
scale_mat = ((d_mat./dR).^m).*exp((E/R).*(1./T_K_mat-1/TR)).*exp((Vstar/R).*(P_Pa_mat./T_K_mat-PR/TR)) ; % more like viscosity, so divide by melt factor

%(Xtilde is like strain rate!, where you multiply by rate factor)
scale_mat = scale_mat.*x_phi_c ; % to account for lack of truly melt free samples and the drop at the onset of melting.

% melt enhancement
[scale_mat_prime] = sr_melt_enhancement(phi,alpha,x_phi_c,phi_c) ;

scale_mat = scale_mat./scale_mat_prime ;

% use linear indexing!! will loop over n-dimensions of Ju_mat.
n_th = numel(Ju_mat); % number of thermodynamic states
for x1 = 1:n_th; % loop using linear index!

%   pull out variables at current index
    scale = scale_mat(x1);
    Ju = Ju_mat(x1) ;
    rho = rho_mat(x1) ;

%   loop over frequency
    for i=1:nfreq
        i_glob = x1 + (i - 1) * n_th; % the linear index of the arrays with
                                         % a frequency index
        w = w_vec(i) ;
        Tau_L = Tau_LR.*scale ;
        Tau_H = Tau_HR.*scale ;

        ntau = 200 ;
        Tau_X_vec = logspace(log10(Tau_L),log10(Tau_H),ntau) ;

        % maxwell relaxation time (period)
        Tau_M = Tau_MR.*scale ;

        if method==0
        %% MY METHOD--
            D_vec = (alf.*Tau_X_vec.^(alf-1))./(Tau_H^alf - Tau_L^alf) ;

            int_J1 = trapz(Tau_X_vec,(D_vec./(1+w^2.*Tau_X_vec.^2))) ;
            J1(i_glob) = Ju.*(1+Delta.*int_J1) ;

            int_J2 = trapz(Tau_X_vec,((Tau_X_vec.*D_vec)./(1+w^2.*Tau_X_vec.^2))) ;
            J2(i_glob) = Ju.*(w*Delta*int_J2 + 1/(w*Tau_M)) ;

        elseif method==1
        %% GEOFF's METHOD-- works !
        % GA: try iterative quadrature
        % (still seems to be a problem with high freq, w*tauH > 1 say)
            Tau_fac = alf.*Delta./(Tau_H.^alf - Tau_L.^alf);

            FINT1 = @(x) (x.^(alf-1))./(1+(w.*x).^2);
            int1 = Tau_fac.*quad(FINT1, Tau_L, Tau_H);

            FINT2 = @(x) (x.^alf)./(1+(w.*x).^2);
            int2 = w.*Tau_fac.*quad(FINT2, Tau_L, Tau_H);

            J1(i_glob) = Ju.*(1 + int1);
            J2(i_glob) = Ju.*(int2 + 1./(w.*Tau_M));
        end

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
