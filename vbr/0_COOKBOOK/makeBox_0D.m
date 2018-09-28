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

filename = 'VBR_GIA_LUT_lab.mat'
%filename = 'VBR_T_gs_melt_LUT.mat'
%% ====================================================
%% Load and set VBR parameters ========================
%% ====================================================

%  put VBR in the path
   VBR_version = 'VBR_v0p95';
   addpath(genpath(['../4_VBR/',VBR_version ])); % recursive add path

%  write method list (these are the things to calculate)
%  all methods will end up as output like:
%      VBR.out.elastic.anharmonic, VBR.out.anelastic.eBurgers, etc.
   VBR.in.elastic.methods_list={'anharmonic';'poro_Takei'}; %;'SLB2005'};
   VBR.in.viscous.methods_list={'HK2003'; 'LH2012'};
   VBR.in.anelastic.methods_list={'eBurgers';'AndradePsP';'YT2016_solidus';'YT_maxwell'};

%  load anharmonic parameters, adjust Gu_0_ol and derivatives to match YT2016
   VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
   VBR.in.elastic.anharmonic.Gu_0_ol=72.45; %[GPa]
   VBR.in.elastic.anharmonic.dG_dT = -10.94*1e6; % Pa/C    (equivalent ot Pa/K)
   VBR.in.elastic.anharmonic.dG_dP = 1.987; % GPa / GPa

%  frequencies to calculate at
fmin = 1/(1e4*pi*1e7) % 10,000 years period
fmax = 1/1 % 1 Hz !
f_vec = logspace(log10(fmin),log10(fmax),100)
VBR.in.SV.f = f_vec ;

%% ====================================================
%% define variable vectors =====================
%% ====================================================
T_C_vec = 800:50:1500 ;
gs_um_vec = logspace(log10(1),log10(1e4),20);
%phi_vec = linspace(0,0.05,25) ;

VBR.in.SV_vectors.T_K_vec_dim1 = T_C_vec+273 ;
VBR.in.SV_vectors.gs_um_vec_dim2 = gs_um_vec ;
%VBR.in.SV_vectors.phi_vec_dim3 = phi_vec ;

%[T_K_ra,gs_um_ra] = ndgrid(T_C_vec+273, gs_um_vec) ;
[T_K_ra,gs_um_ra] = ndgrid(T_C_vec+273,gs_um_vec) ;
%T_K_ra = T_K_ra'
oneses = ones(size(T_K_ra)) ; %,len(gs_um_vec),len(phi_vec));
sz=size(oneses)  %

%% ====================================================
%% Define the Thermodynamic State ARRAYS=====================
%% ====================================================
% indexes:
% 1 = temperature
% 2 = grain size
% 3 = melt fraction
% add stress?
% 4 = frequency

% for grid search, add melt fraction, grain size

%  size of the state variable arrays. arrays can be any shape
%  but all arays must be the same shape.
   %VBR.in.SV.T_C = T_K_ra;
   VBR.in.SV.T_K= T_K_ra ; % VBR.in.SV.T_C+273; % temperature [K]
   VBR.in.SV.dg_um = gs_um_ra;

%  intensive state variables (ISV)
% laboratory conditions
   VBR.in.SV.P_GPa = 0.3 * oneses; % pressure [GPa]
   VBR.in.SV.rho = 3300 * oneses; % density [kg m^-3]
   VBR.in.SV.sig_MPa = 10 * oneses; % differential stress [MPa]

% asthenosphere conditions
   %VBR.in.SV.dg_um = 5e3*oneses;
   % VBR.in.SV.P_GPa = 0.3 * oneses; % pressure [GPa]
   % VBR.in.SV.rho = 3300 * oneses; % density [kg m^-3]
   % VBR.in.SV.sig_MPa = 0.5 * oneses; % differential stress [MPa]

   VBR.in.SV.chi = oneses; % composition fraction  1 = olivine, 0 = crust

%  structural state variables (SSV)
   VBR.in.SV.phi = 0.0 * oneses; % melt fraction

%  compositional state variables (CSV)
   VBR.in.SV.Ch2o = 0 * oneses ; % water concentration

%  this method requires the solidus
%  you should write your own function for the solidus that takes all the other
%  state variables as input. This is just for illustration
   dTdz=0.5 ; % solidus slope [C/km]
   dTdP=dTdz / 3300 / 9.8 / 1000 * 1e9; % [C/GPa ]
   VBR.in.SV.Tsolidus_K=1000+dTdP*VBR.in.SV.P_GPa;

%% ====================================================
%% CALL THE VBR CALCULATOR ============================
%% ====================================================

   [VBR] = VBR_spine(VBR)

   % for the python version ! (because it will not read the reserved word 'in')
      VBR.input = VBR.in

save(filename,'VBR')
%% ====================================================
%% Display some things ================================
%% ====================================================
i1_T = 5 ;
i2_gs = 3 ;

close all;
figure;
subplot(1,3,1)
%semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT2016_solidus.M(1,:,:)/1e9), 'k'); hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT_maxwell.M(i1_T,:,:)/1e9), 'r');hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.AndradePsP.Ma(i1_T,:,:)/1e9), 'b' );hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.eBurgers.M(i1_T,:,:)/1e9), 'g' );hold on;
ylabel('M [GPa]'); xlabel('period [s]')
ylim([0,80])
axis('tight')

subplot(1,3,2)
%loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT2016_solidus.Qinv(1,:,:)), 'k');hold on;
loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT_maxwell.Qinv(i1_T,:,:)), 'r');hold on;
loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.AndradePsP.Qinv(i1_T,:,:)), 'b');hold on;
loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.eBurgers.Qinv(i1_T,:,:)), 'g');hold on;
ylabel('Q^-1'); xlabel('period [s]')
ylim([1e-3,.1])
axis('tight')

subplot(1,3,3)
%semilogx(1./VBR.in.SV.f,1e-3*squeeze(VBR.out.anelastic.YT2016_solidus.V(1,:,:)), 'k');hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT_maxwell.V(i1_T,:,:))./1e3, 'r');hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.AndradePsP.Va(i1_T,:,:)./1e3), 'b');hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.eBurgers.V(i1_T,:,:)./1e3), 'g');hold on;
ylabel('V_s [km/s]'); xlabel('period [s]')
axis('tight')
