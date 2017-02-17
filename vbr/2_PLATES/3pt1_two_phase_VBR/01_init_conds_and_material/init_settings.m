function [settings]=init_settings
    
% melt migration physics
      settings.Vbg = 0; % [cm/yr]  
%     exponents
      settings.nxi = -1; % compaction viscosity porosity exponent, xi ~ phi^nxi
      settings.np = 3; % permeability porosity exponent, k ~ phi^np
      settings.n1 = 1; % flag for 1-phi approximation. set to 0 for 1-phi = 1.
%     physical quantities     
      settings.mufo = 1; % fluid shear viscosity [Pa s]
      settings.muso = 1e19; % Solid shear viscosity at phi = 0 [Pa s]
      settings.grain0 = 0.01; % grain size [m]
      settings.C = 300; % tortuosity
      settings.alpha = 40; % shear viscosity melt decay constant, exp(-alpha phi).
      settings.y_sl = 1*5; % solid-liquid surface energy [J/m2] (value is a guess)

% crustal thickness
  settings.Z_moho_km = 20; % Moho depth [km]
  settings.Moho_thickness_km = 6; % lengthscale for the gradient from crust to mantle [km]
  settings.Z_LAB_max = 100000; % will quit when zLAB exceeds this value [km]
  
% density settings (THERMAL PROBLEM USES VARIABLE DENSITY, MELT MIGRATION 
% USES UNIFORM SOLID DENSITY RIGHT NOW). 
  settings.rhos = 3300; % solid density [kg/m3]
  settings.rhos_crust = 2800; % crustal density [km/m3]
  settings.rhof = 2800; % fluid density [kg/m3]

% chemical transport 
  settings.kd_H2O = 1e-2; % partition coefficent for H2O, kd = Cs / Cf
  settings.Cs0 = 300 * (1e-6 * 1e2); % initial water concentration [wt %]
  settings.Cs_Df = 1e-8; % dispersion coefficient in fluid [m2/s]
  settings.Cs_Ds = 1e-12; % diffusivity in solid [m2/s]
  
  settings.kd_CO2 = 1e-4; % partition coefficent for CO2, kd = Cs / Cf
  settings.Cs0_CO2 = 0 * (1e-6 * 1e2); % initial CO2 concentration [wt %]    
  
  settings.kd_parg = 10; % effective partition coefficient for H2O in 
                          % pargasite (made up so that water prefers to 
                          % go into the solid when amphibole is around). 
  settings.parg_file='pargasite_grabit_fit.mat';
% thermal properties at reference state
  settings.Kc_olv = 4.17; % thermal conductivity [W/m/K]
  settings.Kc_crust = 4.17; % thermal conductivity [W/m/K]
  settings.Cp_olv = 1100; % heat capacity [J/kg/K]
  settings.Cp_crust = 1100; % heat capacity [J/kg/K]
  settings.L = 500 * 1e3; % latent heat of crystallization [J/kg]

% elasticity (reference unrelaxed values, all other parameters
% related to elasticity are set in the VBR param files)
  settings.Gu_crust = 44.0 ; %GPa
  settings.Gu_olv = 66.5 ; %GPa
     
% other physics
  settings.g = 9.8; % gravtational acceleration [m/s2]
  settings.gz = 1; % gravity flag          
  settings.P0 = 1e5; % pressure at surface [Pa]

% adiabat (could be calculated self consistently...) 
  settings.dTdz_ad = 0.5*1e-3; % adiabatic gradient [K/m]
  %settings.dTdP_ad = settings.dTdz_ad/settings.rhos/settings.g; % adiabtaic gradient [K/Pa]
  
% used for diking flux
  settings.Sd_coefficient = 1 * (settings.rhos - settings.rhof) * settings.g ...
                            * settings.grain0^2 / settings.mufo / settings.C;

% used for other things....
  settings.DBL_m = 5e3; % [m] mechancial boundary layer for upwelling V
  settings.phi_init = 0.005;
  settings.Cs0 = 0 * 1e-4; % H2O concentration [wt % -- needs to be wt % for solidus]
  settings.Cs0_CO2 = 0 *1e-4; % CO2 concentration [wt %]
  settings.y_sl = 1;  
  settings.Sd_coefficient = 1e-5;    
  settings.DikingBL_km = 10; % [km] diking process zone boundary layer     
  settings.XtalFactor = 1; % fraction of infiltrated melt that freezes
  settings.Tpot=1325; % [C]
  settings.Tpot_excess=1325; % [C]


%%% disequilbrium melting (experimental... need to uncomment section in TwoPhase.m to use this) 
%     settings.Dam0 = 10; % true Damkohler number
%     to = 3600 * 24 * 365 * 1e6; % time scale = 1 Myr? 
%     settings.Dam = settings.Dam0 / to; % dimensional Damkohler number

%  Mesh settings
   settings.Zinfo.dz = .5; % node spacing in [km]
   settings.Zinfo.zmin =0; % min depth for model domain [km]
   settings.Zinfo.zmax = 100; % max z depth [km]   

%  Computational settings 
   settings.nt =20; % max time steps (switch to a max time)
   settings.outk =5; % output frequency
   settings.phimin = 1e-10; % minimum phi
   settings.phimax = 0.1; % maximum phi
   settings.sstol = 1e-20; % steady state target residual  
   
%  all the flags with default settings 
%  Flags
%    .problem   'One_Phase_T' = single phase energy conservation only
%               'Two_Phase_Sd' = dike-flux melt infiltration
%    .PropType  specifies dependencies of Kc, rho and Cp.           
%               'con'     Constant rho, Cp and Kc
%               'P_dep'   pressure dependent rho, constant Cp and Kc
%               'T_dep'   temperature dependent rho, Cp and Kc
%               'PT_dep'  temperature and pressure dependent rho, Cp and Kc
%    .LABdef    sets the method for calculating the "LAB." The melt 
%               fraction evolution always uses the phi-LAB, this flags is 
%               for the upwelling velocity.              
%               'visc'    eta/eta_astheno > 10
%               'phi'     first point above solidus
%    .VbzFlag   sets the method for caculating or setting Vbg
%               'constant'   constant value (i.e., V +S = Vo)
%               'variable'   tapers to 0 at LAB, recalculates as LAB moves
%               'var_z_con_t' tapers to 0 at LAB, doesn't recalculate
%    .H2O       'constant' or 'variable' 
%    .XtalFactor 'constant' or 'variable'
%    .problemkill  'stop_on_no_melt' or 'keepgoing_on_no_melt'
%    .phikill    'stop_on_max_phi' or 'keepgoing'

   settings.Flags.VbzFlag = 'constant'; % 'constant' or 'variable'  or 'var_z_con_t'
   settings.Flags.PropType='PT_dep';  % for density, thermal exp., conductivity
   settings.Flags.H2O='constant';
   settings.Flags.progress_plot = 'no';
   settings.Flags.parg = 'no';
   settings.Flags.problem='Two_Phase_Sd';
   settings.Flags.LABdef='phi';
   settings.Flags.XtalFactor = 'variable';
   settings.Flags.problemkill='stop_on_no_melt';
   settings.Flags.phikill='stop_on_max_phi';


  
end
