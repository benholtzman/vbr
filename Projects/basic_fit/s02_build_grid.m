%% ===================================================================== %%
%%                     fit_scalar.m
%% ===================================================================== %%
%  simple grid search for scalar velo
%% ===================================================================== %%

  clear

% put VBR in the path
  path_to_top_level_vbr='../../';
  addpath(path_to_top_level_vbr)
  vbr_init

  VBR_file='./VBR_output.mat';

%% set grid search
   T_K=1300:25:1600 + 273; %
   phi=logspace(-8,-1.3,20); % up to 0.05
   d_um=logspace(-4,-1.5,15)*1e6; % up to 3 cm
   [VBR.in.SV.T_K,VBR.in.SV.phi,VBR.in.SV.dg_um]=meshgrid(T_K,phi,d_um);
   sz=size(VBR.in.SV.T_K); % temperature [K]
   rho=3300;
   P_GPa=80*1e3*9.8*rho / 1e9;

 %  frequencies to calculate at
     periods = [50, 100, 150];
    VBR.in.SV.f = 1./periods;

%% ====================================================
%% Load and set VBR parameters ========================
%% ====================================================


%  write method list (these are the things to calculate)
%  all methods will end up as output like:
%      VBR.out.elastic.anharmonic, VBR.out.anelastic.eBurgers, etc.
   VBR.in.elastic.methods_list={'anharmonic'};
   VBR.in.anelastic.methods_list={'eBurgers';'AndradePsP';'YT_maxwell'};

   % uncomment to use VBR viscosity for maxwell time calculation:
   % VBR.in.anelastic.eBurgers=Params_Anelastic('eBurgers');
   % VBR.in.viscous.methods_list={'LH2012'};
   % VBR.in.anelastic.eBurgers.useJF10visc=0;

%  load anharmonic parameters, adjust Gu_0_ol
%  all paramss in ../4_VBR/VBR_version/params/ will be loaded in call to VBR spine,
%  but you can load them here and adjust any one of them (rather than changing those
%  parameter files).
   VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity

   % JF10 have Gu_0=62.5 GPa, but that's at 900 Kelvin and 0.2 GPa,
   % so set Gu_0_ol s.t. it ends up at 62.5 at those conditions
   dGdT=VBR.in.elastic.anharmonic.dG_dT;
   dGdP=VBR.in.elastic.anharmonic.dG_dP;
   Tref=VBR.in.elastic.anharmonic.T_K_ref;
   Pref=VBR.in.elastic.anharmonic.P_Pa_ref/1e9;
   VBR.in.elastic.anharmonic.Gu_0_ol = 62.5 - (900+273-Tref) * dGdT/1e9 - (0.2-Pref)*dGdP; % olivine reference shear modulus [GPa]

   % VBR.in.anelastic.eBurgers=Params_Anelastic('eBurgers');
   % VBR.in.anelastic.eBurgers.eBurgerMethod='bg_peak';
   % VBR.in.elastic.anharmonic.Gu_0_ol = 66.5 - (900+273-Tref) * dGdT/1e9 - (0.2-Pref)*dGdP;



%% ====================================================
%% Define the Remainig SVs
%% ====================================================


   VBR.in.SV.rho = rho * ones(sz); % density [kg m^-3]
   VBR.in.SV.P_GPa = P_GPa * ones(sz); % pressure [GPa]
   VBR.in.SV.sig_MPa = 0.1 * ones(sz); % differential stress [MPa]
   VBR.in.SV.chi=ones(sz); % composition fraction  1 = olivine, 0 = crust
   VBR.in.SV.Ch2o = 0 * ones(sz) ; % water concentration

%% ====================================================
%% CALL THE VBR CALCULATOR ============================
%% ====================================================

   [VBR] = VBR_spine(VBR) ;

  save(VBR_file,'VBR')
