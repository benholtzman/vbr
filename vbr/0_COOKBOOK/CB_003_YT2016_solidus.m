%% ===================================================================== %%
%%                     CB_003_YT2016_solidus.m
%% ===================================================================== %%
%  Calls VBR using YT2016_solidus method from:
%  Hatsuki Yamauchi and Yasuko Takei, JGR 2016, "Polycrystal anelasticity at
%  near-solidus temperatures,"
%
%  sets elastic parameters to match their results
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
   VBR.in.elastic.methods_list={'anharmonic'};
   VBR.in.anelastic.methods_list={'YT2016_solidus'};

%  load anharmonic parameters, adjust Gu_0_ol and derivatives to match YT2016
   VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
   VBR.in.elastic.anharmonic.Gu_0_ol=72.45; %[GPa]
   VBR.in.elastic.anharmonic.dG_dT = -10.94*1e6; % Pa/C    (equivalent ot Pa/K)
   VBR.in.elastic.anharmonic.dG_dP = 1.987; % GPa / GPa

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

%  this method requires the solidus
%  you should write your own function for the solidus that takes all the other
%  state variables as input. This is just for illustration
   dTdz=0.5 ; % solidus slope [C/km]
   dTdP=dTdz / 3300 / 9.8 / 1000 * 1e9; % [C/GPa ]
   VBR.in.SV.Tsolidus_K=1000+dTdP*VBR.in.SV.P_GPa;

%% ====================================================
%% CALL THE VBR CALCULATOR ============================
%% ====================================================

   [VBR] = VBR_spine(VBR) ;

%% ====================================================
%% Display some things ================================
%% ====================================================

close all;
figure;
subplot(1,3,1)
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT2016_solidus.M(1,:,:)/1e9));
ylabel('M [GPa]'); xlabel('period [s]')
ylim([0,80])

subplot(1,3,2)
loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT2016_solidus.Qinv(1,:,:)));
ylabel('Q^-1'); xlabel('period [s]')
ylim([1e-3,.1])

subplot(1,3,3)
semilogx(1./VBR.in.SV.f,1e-3*squeeze(VBR.out.anelastic.YT2016_solidus.V(1,:,:)));
ylabel('V_s [km/s]'); xlabel('period [s]')
