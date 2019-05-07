%% ===================================================================== %%
%%                     CB_007_elastic.m
%% ===================================================================== %%
%  Calculates anharmonic velocities, SLB parametrization 
%% ===================================================================== %%
   clear

%% ====================================================
%% Load and set VBR parameters ========================
%% ====================================================

% put VBR in the path
  path_to_top_level_vbr='../../';
  addpath(path_to_top_level_vbr)
  vbr_init

%  write method list (these are the things to calculate)
%  all methods will end up as output like:
%      VBR.out.elastic.anharmonic, VBR.out.anelastic.eBurgers, etc.
   VBR.in.elastic.methods_list={'anharmonic','poro_Takei','SLB2005'};

%% ====================================================
%% Define the Thermodynamic State =====================
%% ====================================================

%  size of the state variable arrays. arrays can be any shape
%  but all arays must be the same shape.
   rho=3300;
   z=linspace(0,300,100)*1e3;
   VBR.in.SV.P_GPa = rho*9.8*z / 1e9;

   zPlate=100*1e3;
   dTdz=0.6 / 1000 ; % deg/m
   VBR.in.SV.T_K = z/zPlate * 1300;
   VBR.in.SV.T_K(z>zPlate)=1300+(z(z>zPlate)-zPlate)*dTdz;
   VBR.in.SV.T_K=VBR.in.SV.T_K+273;

   sz=size(z);
   VBR.in.SV.phi=zeros(sz);
   VBR.in.SV.phi(z>zPlate)=0.01;

%  intensive state variables (ISV)
   VBR.in.SV.rho = rho * ones(sz); % density [kg m^-3]
   VBR.in.SV.sig_MPa = 10 * ones(sz); % differential stress [MPa]

%  compositional state variables (CSV)
   VBR.in.SV.Ch2o = 0 * ones(sz) ; % water concentration [wt%]
   VBR.in.SV.chi = 1 * ones(sz) ; %

%% ====================================================
%% CALL THE VBR CALCULATOR ============================
%% ====================================================

   [VBR] = VBR_spine(VBR) ;
