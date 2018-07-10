%% ===================================================================== %%
%%                     CB_001_0D_scalar.m
%% ===================================================================== %%
%  Calls VBR using a single thermodynamic state
%% ===================================================================== %%
   clear

%% ====================================================
%% Load and set VBR parameters ========================
%% ====================================================

%  put VBR in the path
   VBR_version = 'VBR_v0p95';
   addpath(genpath(['../4_VBR/',VBR_version ])); % recursive add path

%  write method list (these are the things to calculate)
%  all methods will end up as output like:
%      VBR.out.elastic.anharmonic, VBR.out.anelastic.eBurgers, etc.
  %VBR.in.elastic.methods_list={'anharmonic';'poro_Takei';'SLB2005'};
   VBR.in.elastic.methods_list={'anharmonic'};
   VBR.in.viscous.methods_list={'HK2003'; 'LH2012'}; 
   VBR.in.anelastic.methods_list={'eBurgers';'AndradePsP'};

%  load anharmonic parameters, adjust Gu_0_ol
%  all paramss in ../4_VBR/VBR_version/params/ will be loaded in call to VBR spine,
%  but you can load them here and adjust any one of them (rather than changing those
%  parameter files).
   VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity



   % JF10 have Gu_0=62 GPa, but that's at 900 Kelvin and 0.2 GPa,
   % so set Gu_0_ol s.t. it ends up at 62 at those conditions
   dGdT=VBR.in.elastic.anharmonic.dG_dT;
   dGdP=VBR.in.elastic.anharmonic.dG_dP;
   Tref=VBR.in.elastic.anharmonic.T_K_ref;
   Pref=VBR.in.elastic.anharmonic.P_Pa_ref/1e9;
   VBR.in.elastic.anharmonic.Gu_0_ol = 66.5 - (900+273-Tref) * dGdT/1e9 - (0.2-Pref)*dGdP; % olivine reference shear modulus [GPa]


%  frequencies to calculate at
   VBR.in.SV.f = 1./logspace(-2,4,100);

%% ====================================================
%% Define the Thermodynamic State =====================
%% ====================================================

%  size of the state variable arrays. arrays can be any shape
%  but all arays must be the same shape.
   VBR.in.SV.T_K=700:50:1200;
   VBR.in.SV.T_K=VBR.in.SV.T_K+273;
   sz=size(VBR.in.SV.T_K); % temperature [K]

%  intensive state variables (ISV)
   VBR.in.SV.dg_um=3.1*ones(sz);
   VBR.in.SV.P_GPa = 0.2 * ones(sz); % pressure [GPa]
   VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
   VBR.in.SV.sig_MPa = 10 * ones(sz); % differential stress [MPa]
   VBR.in.SV.chi=ones(sz); % composition fraction  1 = olivine, 0 = crust

%  structural state variables (SSV)
   VBR.in.SV.phi = 0.0 * ones(sz); % melt fraction

%  compositional state variables (CSV)
   VBR.in.SV.Ch2o = 0 * ones(sz) ; % water concentration

%% ====================================================
%% CALL THE VBR CALCULATOR ============================
%% ====================================================

   [VBR] = VBR_spine(VBR) ;

%% ====================================================
%% Display some things ================================
%% ====================================================

close all;
figure;
subplot(1,2,1)
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.eBurgers.M(1,:,:)/1e9));
ylabel('M [GPa]'); xlabel('period [s]')
ylim([0,80])

subplot(1,2,2)
loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.eBurgers.Qinv(1,:,:)));
ylabel('Q^-1'); xlabel('period [s]')
ylim([3e-3,3])
