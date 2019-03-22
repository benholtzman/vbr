function[VBR]=Q_eBurgers_f(VBR)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %
% extended burgers model after JF10
%
% JF10:
% Jackson & Faul, "Grainsize-sensitive viscoelastic relaxation in
% olivine: Towards a robust laboratory-based model for seismological
% application," Physics of the Earth and Planetary Interiors 183 (2010) 151–163
%
% includes high temperature background, optional dissipation peak.
% see Projects/vbr_core_examples/CB_test_tmep.m for using with/without peak
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %

  % angular frequency
  f_vec = VBR.in.SV.f ;
  w_vec = 2*pi.*f_vec ;

  % unrelaxed compliance and density
  if isfield(VBR.in.elastic,'poro_Takei')
   Mu = VBR.out.elastic.poro_Takei.Gu ;
  elseif isfield(VBR.in.elastic,'anharmonic')
   Mu = VBR.out.elastic.anharmonic.Gu ;
  end
  Ju_mat = 1./Mu ;
  rho_mat = VBR.in.SV.rho ; % density

  % allocation (frequency is added as a new dimension at end of array.)
  nfreq = numel(f_vec);
  J1 = proc_add_freq_indeces(zeros(size(Ju_mat)),nfreq);
  J2 = J1; Q = J1; Qinv = J1; M = J1; V = J1;
  Vave = zeros(size(Ju_mat));

  % Calculate maxwell time, integration limits and location of peak:
  tau=MaxwellTimes(VBR,Mu);

  % Read in parameters needed for integration
  Burger_params=VBR.in.anelastic.eBurgers;
  bType=Burger_params.eBurgerMethod;
  alf = Burger_params.(bType).alf ;
  DeltaB = Burger_params.(bType).DeltaB ; % relaxation strength of background
  DeltaP=Burger_params.(bType).DeltaP; % relaxation strength of peak
  sig=Burger_params.(bType).sig;
  HTB_int_meth=Burger_params.integration_method ; % (trapezoidal, 0; quadrature, 1)
  ntau = Burger_params.tau_integration_points ;

  % ============================================================================
  % loop over state varialbes, calculate anelatic effects at every frquency
  % ============================================================================

  n_th = numel(Ju_mat); % number of thermodynamic states
  for x1 = 1:n_th;

    Ju = Ju_mat(x1) ; % unrelaxed compliance
    rho = rho_mat(x1) ; % density

    % maxwell times
    Tau_M = tau.maxwell(x1);
    Tau_L = tau.L(x1);
    Tau_H = tau.H(x1);
    Tau_P = tau.P(x1);
    if HTB_int_meth == 0
      Tau_X_vec = logspace(log10(Tau_L),log10(Tau_H),ntau) ;
    end

    % loop over frequency
    for i=1:nfreq
      i_glob = x1 + (i - 1) * n_th; % the linear index of the arrays with a frequency index
      w = w_vec(i);

      if HTB_int_meth==0 %% trapezoidal integration --
          D_vec = (alf.*Tau_X_vec.^(alf-1))./(Tau_H^alf - Tau_L^alf);

          int_J1 = trapz(Tau_X_vec,(D_vec./(1+w^2.*Tau_X_vec.^2)));
          J1(i_glob) = (1+DeltaB.*int_J1);

          int_J2 = trapz(Tau_X_vec,((Tau_X_vec.*D_vec)./(1+w^2.*Tau_X_vec.^2)));
          J2(i_glob) = (w*DeltaB*int_J2 + 1/(w*Tau_M));

      elseif HTB_int_meth==1 % iterative quadrature
          Tau_fac = alf.*DeltaB./(Tau_H.^alf - Tau_L.^alf);

          FINT1 = @(x) (x.^(alf-1))./(1+(w.*x).^2);
          int1 = Tau_fac.*quadl(FINT1, Tau_L, Tau_H);

          FINT2 = @(x) (x.^alf)./(1+(w.*x).^2);
          int2 = w.*Tau_fac.*quadl(FINT2, Tau_L, Tau_H);

          J1(i_glob) = (1 + int1);
          J2(i_glob) = (int2 + 1./(w.*Tau_M));
      end

      % add on peak if it's being used.
      % May trigger warning, this integral is not easy.
      if DeltaP>0
        FINT2 = @(x) (exp(-(log(x./Tau_P)/sig).^2/2)./(1+(w.*x).^2));
        int2a = quadgk(FINT2, 0, inf);
        J2(i_glob)=J2(i_glob)+DeltaP*w*(int2a)/(sig*sqrt(2*pi));

        FINT1 = @(x) ( 1./x .* exp(-(log(x./Tau_P)/sig).^2/2)./(1+(w.*x).^2));
        int1 = quadgk(FINT1, 0, inf);
        J1(i_glob)=J1(i_glob)+DeltaP*int1 / (sig*sqrt(2*pi)) ;
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
    end % end loop over frequency
  end % end the loop(s) over spatial dimension(s)
  % ============================================================================

  % Store relevant values
  VBR.out.anelastic.eBurgers.J1 = J1;
  VBR.out.anelastic.eBurgers.J2 = J2;
  VBR.out.anelastic.eBurgers.Q = Q;
  VBR.out.anelastic.eBurgers.Qinv = Qinv;
  VBR.out.anelastic.eBurgers.M=M;
  VBR.out.anelastic.eBurgers.V=V;
  VBR.out.anelastic.eBurgers.Vave = Vave./nfreq;

end

function tau=MaxwellTimes(VBR,Mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculatues the maxwell time & limits for extended burgers model:
% tau.maxwell = steady state viscous maxwell time (i.e., eta / Gu)
% tau.L = lower limit of integration for high temp background
% tau.H = upper limit of integration for high temp background
% tau.P = center period of dissipation peak (if being used)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Burger_params=VBR.in.anelastic.eBurgers;
  bType=Burger_params.eBurgerMethod;

  % state variables for either maxwell time or integration limits, peak loc:
  phi =  VBR.in.SV.phi ;
  T_K_mat = VBR.in.SV.T_K ; % temperature [K]
  P_Pa_mat = VBR.in.SV.P_GPa.*1e9 ; % convert pressure GPa to Pa = GPa*1e9
  d_mat = VBR.in.SV.dg_um ; % microns grain size

  % Scaling values from JF10
  TR = Burger_params.(bType).TR ;% Kelvins
  PR = Burger_params.(bType).PR *1e9; % convert pressure GPa to Pa = GPa*1e9
  dR = Burger_params.(bType).dR ; % microns grain size
  E = Burger_params.(bType).E ; % activation energy J/mol
  R = Burger_params.R ; % gas constant
  Vstar = Burger_params.(bType).Vstar ; % m^3/mol Activation Volume
  m_a = Burger_params.(bType).m_a ; % grain size exponent (anelastic)
  m_v = Burger_params.(bType).m_v ; % grain size exponent (viscous)

  % maxwell time calculation
  [visc_exists,missing]=checkStructForField(VBR,{'in','viscous','methods_list'},0);
  if Burger_params.useJF10visc || visc_exists==0
    % use JF10's exact relationship
    scale=((d_mat./dR).^m_v).*exp((E/R).*(1./T_K_mat-1/TR)).*exp((Vstar/R).*(P_Pa_mat./T_K_mat-PR/TR));
    scale=addMeltEffects(phi,scale,VBR.in.GlobalSettings,Burger_params);
    Tau_MR = Burger_params.(bType).Tau_MR ;
    tau.maxwell=Tau_MR.*scale ; % steady state viscous maxwell time
  else
    % use diffusion viscosity from VBR to get maxwell time
    visc_method=VBR.in.viscous.methods_list{1};
    eta_diff = VBR.out.viscous.(visc_method).diff.eta ; % viscosity for maxwell relaxation time
    tau.maxwell = eta_diff./ Mu ; % maxwell relaxtion time
  end

  % integration limits and peak location
  LHP=((d_mat./dR).^m_a).*exp((E/R).*(1./T_K_mat-1/TR)).*exp((Vstar/R).*(P_Pa_mat./T_K_mat-PR/TR));
  LHP=addMeltEffects(phi,LHP,VBR.in.GlobalSettings,Burger_params);
  tau.L = Burger_params.(bType).Tau_LR * LHP;
  tau.H = Burger_params.(bType).Tau_HR * LHP;
  tau.P = Burger_params.(bType).Tau_PR * LHP;

end

function scaleMat=addMeltEffects(phi,scaleMat,GlobalSettings,Burger_params)
  % adds on Melt Effects.

  % sharper response than creep at phi_c, but less sensitive at higher phi (for HTB only)
  alpha = Burger_params.melt_alpha ; % post-critical melt fraction dependence
  phi_c = Burger_params.phi_c ; % critical melt fraction
  x_phi_c = Burger_params.x_phi_c ;% melt enhancement factor

  % x_phi_c adjustment ("nominally melt free" to truly melt free)
  if GlobalSettings.melt_enhacement==0
    x_phi_c=1;
  else
    scaleMat = scaleMat.* x_phi_c ;
  end

  % add melt effects
  [scale_mat_prime] = sr_melt_enhancement(phi,alpha,x_phi_c,phi_c) ;
  scaleMat = scaleMat ./ scale_mat_prime;
end
