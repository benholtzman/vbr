function params = Params_Anelastic(method)

%% ========================================================================
%% Anelastic Properties ===================================================
%% ========================================================================


%===========  extended BURGERS parameters =================================
% citation:
  if strcmp(method,'eBurgers')
      % REFERENCE Temp, Pressure and d (grain size)
      params.TR = 1173;% Kelvins
      params.PR = 0.2; % confining pressure of experiments, GPa
      params.dR = 3.1 ; % 3.1 microns grain size

      % REFERENCE Modulus (not the same as the reference Gu_0 for unrelaxed
      % calculations, it is the reference modulus for scaling experiments)
      params.G_UR = 62 ; % GPa, unrel. G, reference val. # NOT USED RIGHT NOW

      params.E = 303000 ; % J/mol
      params.R = 8.314 ; % gas constant
      params.Vstar = 10e-6 ; % m^3/mol (Activation Volume? or molar volume?)
      params.m = 1 ;
      params.method='PointWise'; % 'FastBurger' or 'PointWise' or 'None'
      params.nTauGlob=3000; % points for global Tau discretization ('FastBurger' ONLY)
      
      %% THESE ARE THE THINGS FROM JF10
      %   %gsr = 1.34E-5; % reference grain size in m
      %   %deltaB = 1.04; % background relaxation strength,
      %   %alpha = 0.274; % background frequency exponent
      %   %%reference values for relaxation times
      %   %tauLo = 1E-3; tauHo = 1E7; tauMo = 3.02E7;
      %   %ma = 1.31;  % anelastic grain size exponent
      %   %mv = 3;     % viscous grain size exponent
      %   %EB = 3.6E5; % activation energy for background and peak
      %   %AV = 1E-5;  % activation volume [m3/mol]
      %   %% peak parameters:
      %   %tauPo = 3.98E-4; % reference peak relaxation time,
      %   %deltaP = 0.057;  % peak relaxation strength pref
      %   %sig = 4;         % peak width
      %   %cp = deltaP * (2*pi)^(-0.5)/sig; %peak integration const.

      % Jackson n Faul 2010, table 1 :
      params.alf = 0.33 ; % is this the same as n in Andrade ?
      params.Delta = 1.4 ;% ; % relaxation strength..

      params.Tau_LR = 1e-2 ; % Relaxation time lower limit reference
      params.Tau_HR = 1e6 ; % Relaxation time higher limit reference
      params.Tau_MR = 10^5.2 ; % Reference Maxwell relaxation time

      %% melt effects
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
