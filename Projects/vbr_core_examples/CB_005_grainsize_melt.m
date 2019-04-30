%% ===================================================================== %%
%%                     CB_005_grainsize_melt.m
%% ===================================================================== %%
%  Compares YT2016_solidus and AndradePsp at a range of grain size and
%  melt fractions for a single period (100s) at astheno conditions.
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
   VBR.in.elastic.methods_list={'anharmonic'};
   VBR.in.anelastic.methods_list={'AndradePsP','YT2016_solidus'};

%  frequencies to calculate at
   VBR.in.SV.f = 0.01;

%% ====================================================
%% Define the Thermodynamic State =====================
%% ====================================================

%  grid of dg and phi
   dg_m=logspace(-3,-1.3,50); % grain size [m] (1 mm to 5 cm)
   phi = logspace(-8,-1,20);
   [VBR.in.SV.dg_um,VBR.in.SV.phi]=meshgrid(dg_m*1e6,phi);

%  size of the state variable arrays. arrays can be any shape
%  but all arays must be the same shape.
   sz=size(VBR.in.SV.dg_um);

%  intensive state variables (ISV)
   VBR.in.SV.T_K=(1350+273) * ones(sz); % temperature [K]
   VBR.in.SV.P_GPa = 3.2 * ones(sz); % pressure [GPa]
   VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
   VBR.in.SV.sig_MPa = 0.1 * ones(sz); % differential stress [MPa]
   VBR.in.SV.chi=ones(sz); % composition fraction  1 = olivine, 0 = crust

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
contourf(log10(dg_m),log10(phi),VBR.out.anelastic.YT2016_solidus.V/1e3,100,'LineStyle','none');
colorbar
title('Vs [km/s], YT2016 solidus')
xlabel('d [um]'); ylabel('log10(phi)')

subplot(1,3,2)
contourf(log10(dg_m),log10(phi),VBR.out.anelastic.AndradePsP.V/1e3,100,'LineStyle','none');
colorbar
title('Vs [km/s], AndradePsP')
xlabel('log10(d [um])');

subplot(1,3,3)
dV=(VBR.out.anelastic.YT2016_solidus.V-VBR.out.anelastic.AndradePsP.V) ./ ...
     VBR.out.anelastic.AndradePsP.V;
contourf(log10(dg_m),log10(phi),dV,100,'LineStyle','none');
colorbar
title('( YT2016 - AndradePSP ) / AndradePSP')
xlabel('log10(d [um])');
