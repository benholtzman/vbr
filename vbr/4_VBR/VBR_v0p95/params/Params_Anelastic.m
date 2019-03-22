function params = Params_Anelastic(method)

%% ========================================================================
%% Anelastic Properties ===================================================
%% ========================================================================


%===========  extended BURGERS parameters =================================

  if strcmp(method,'eBurgers')

      params.method='PointWise'; % keep at 'PointWise' until 'FastBurger' fixed
      params.nTauGlob=3000; % points for global Tau discretization ('FastBurger' ONLY)
      params.R = 8.314 ; % gas constant
      params.eBurgerMethod='bg_only'; % 'bg_only' or 'bg_peak'
      params.useJF10visc=1; % if 1, will use the scaling from JF10 for maxwell time.
      params.integration_method=0; % 0 for trapezoidal, 1 for quadrature.
      params.tau_integration_points = 500 ; % number of points for integration of high-T background if trapezoidal

      % values from Table 2 of JF10 all melt-free samples
      % JF10:
      % Jackson & Faul, "Grainsize-sensitive viscoelastic relaxation in
      % olivine: Towards a robust laboratory-based model for seismological
      % application," Physics of the Earth and Planetary Interiors 183 (2010) 151–163

      % best high-temp background only fit:
      params.bg_only.TR=1173; % ref temp [K]
      params.bg_only.PR = 0.2; % ref confining pressure of experiments, GPa
      params.bg_only.dR = 13.4; % ref grain size in microns
      params.bg_only.G_UR = 62.5 ; % GPa, unrel. G, reference val.
      params.bg_only.E = 303000 ; % J/mol
      params.bg_only.Vstar = 10e-6 ; % m^3/mol (Activation Volume? or molar volume?)
      params.bg_only.m_a = 1.19 ; % grain size exponent for tau_i, i in (L,H,P)
      params.bg_only.m_v = 3 ; % viscous grain size exponent for maxwell time
      params.bg_only.alf = 0.257 ; % is this the same as n in Andrade ?
      params.bg_only.DeltaB = 1.13 ;% relaxation strength..
      params.bg_only.Tau_LR = 1e-3 ; % Relaxation time lower limit reference
      params.bg_only.Tau_HR = 1e7 ; % Relaxation time higher limit reference
      params.bg_only.Tau_MR = 10^6.95 ; % Reference Maxwell relaxation time
      params.bg_only.DeltaP=0; % no peak, set to 0
      params.bg_only.sig=0;% no peak, set to 0
      params.bg_only.Tau_PR=0;% no peak, set to 0

      % best high-temp background + peak fit:
      params.bg_peak.DeltaP=0.057; % relaxation strength of peak
      params.bg_peak.sig=4; % sigma, peak breadth
      params.bg_peak.Tau_PR=10^-3.4;
      params.bg_peak.TR=1173; % ref temp [K]
      params.bg_peak.PR = 0.2; % ref confining pressure of experiments, GPa
      params.bg_peak.dR = 13.4; % ref grain size in microns
      params.bg_peak.G_UR = 66.5 ; % GPa, unrel. G, reference val.
      params.bg_peak.E = 360000 ; % J/mol
      params.bg_peak.Vstar = 10e-6 ; % m^3/mol (Activation Volume? or molar volume?)
      params.bg_peak.m_a = 1.31 ; % grain size exponent for tau_i, i in (L,H,P)
      params.bg_peak.m_v = 3 ; % viscous grain size exponent for maxwell time
      params.bg_peak.alf = 0.274 ; % is this the same as n in Andrade ?
      params.bg_peak.DeltaB = 1.13 ;% relaxation strength of background.
      params.bg_peak.Tau_LR = 1e-3 ; % Relaxation time lower limit reference
      params.bg_peak.Tau_HR = 1e7 ; % Relaxation time higher limit reference
      params.bg_peak.Tau_MR = 10^7.48 ; % Reference Maxwell relaxation time

      %% melt effects (don't adjust, use VBR.in.GlobalSettings.melt_enhacement=0 to turn off)
      % use diffusion creep values
      HK2003 = Params_Viscous('HK2003'); % viscous parameters
      params.melt_alpha = HK2003.diff.alf ;
      params.phi_c = HK2003.diff.phi_c ;
      params.x_phi_c = HK2003.diff.x_phi_c ;

  end

%========= ANDRADE (pseudoperiod scaling) parameters (JF10) ===============
% citation:
     if strcmp(method,'AndradePsP')
        params.n = 0.33 ; % 1/3 ;
        params.Beta = 0.020;
        params.Tau_MR = 10^5.3 ;
        params.E = 303e3 ; %J/mol

%       reference values (uses eBurgers values from above)
        params.TR = 1173;% Kelvins
        params.PR = 0.2; % confining pressure of experiments, GPa
        params.dR = 5 ; % 3.1 microns grain size
        params.E = 303e3 ; % J/mol
        params.R = 8.314 ; % gas constant
        params.Vstar = 10e-6 ; % m^3/mol (Activation Volume? or molar volume?)
        params.m = 1 ;

%%      melt effects
%       use diffusion creep values
        HK2003 = Params_Viscous('HK2003'); % viscous parameters
        params.melt_alpha = HK2003.diff.alf ;
        params.phi_c = HK2003.diff.phi_c ;
        params.x_phi_c = HK2003.diff.x_phi_c ;

%       the BUMP
        params.Te = 0.1 ;
        params.Tgbs = 0.0833 ;% sec
        params.Delta = 0.3 ; % 0.43 ; % Relaxation strength
        params.FUDGE = 1 ; %0.7 ;
     end

%========= ANDRADE parameters (from Sundberg+Cooper)=======================
% citation:

     if strcmp(method,'Andrade')
        params.n = 1/3 ; % 1/3 ;

%    scaling option:
%       1= experimental params,
%       2= Marshall's email
%       3= pseudoperiod
        params.scaling_opt=1 ;

        % for 1200, 1250, 1300 C
        %params.A_vec =  [ 5.1e-12 9.83e-12 2.7e-11 ] ;
        %params.eta_ss_vec = [ 3.42e12 8.14e11 1.87e11 ];
        %params.A_vec =  [ 6.37e-12 1.34e-11 2.8e-11 ] ;
        params.A =  2.8e-11 ;
        params.eta_ss =  1.86e11 ;

   %    Bump on (1) or off (0) --
        params.bump = 1 ;
        params.Te = 0.17 ; % not sure what this is
        params.Tgbs = 0.05 ;% sec
        params.Delta = 0.5 ; % 0.43 ; % Relaxation strength
        params.FUDGE = 0.7 ;

     end

%========= YT_maxwell parameters =======================
% citation:
  if strcmp(method,'YT_maxwell')
    params.beta1 = 0.32 ;
    params.beta2 = 1853.0 ;
    params.alpha2 = 0.5 ;
    params.Alpha_a=0.39 ;
    params.Alpha_b=0.28;
    params.Alpha_c=2.6;
    params.Alpha_taun=0.1;
    params.reference='citation';
  end

end
%% =================== END OF Params_Anelastic.m ========================
