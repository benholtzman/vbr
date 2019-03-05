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

%filename = 'VBR_GIA_LUT_asth.mat'
filename = 'VBR_LUT_labdata_y190304.mat'
%filename = 'VBR_T_gs_melt_LUT.mat'
%% ====================================================
%% Load and set VBR parameters ========================
%% ====================================================

%  put VBR in the path
   VBR_version = 'VBR_v0p95';
   addpath(genpath(['../../../vbr/4_VBR/',VBR_version ])); % recursive add path

%  write method list (these are the things to calculate)
%  all methods will end up as output like:
%      VBR.out.elastic.anharmonic, VBR.out.anelastic.eBurgers, etc.
 VBR.in.elastic.methods_list={'anharmonic';'poro_Takei'}; %;'SLB2005'};
 VBR.in.viscous.methods_list={'HK2003'; 'LH2012'};
 VBR.in.anelastic.methods_list={'eBurgers';'AndradePsP';'YT2016_solidus';'YT_maxwell'};
 VBR.in.GlobalSettings.melt_enhacement=0; % turn off critcal melt fraction effect 

%  load anharmonic parameters, adjust Gu_0_ol and derivatives to match YT2016
% FLOWCHART !! WHERE DO THESE PROPAGATE, vs SCALING IN FJ ????
 VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
 %VBR.in.elastic.anharmonic.Gu_0_ol=72.45; %[GPa]

% are these different that our standard? if so, why?
 %VBR.in.elastic.anharmonic.dG_dT = -10.94*1e6; % Pa/C    (equivalent ot Pa/K)
 %VBR.in.elastic.anharmonic.dG_dP = 1.987; % GPa / GPa

 % JF10 have Gu_0=66.5 GPa, but that's at 900 C and 0.2 GPa,
 % so set Gu_0_ol s.t. it ends up at 66.5 at those conditions
   % calculate dGdT from their figure.
   G_900=66.5; % from plot at log10(period)=-2
   G_1200=58; % from plot at log10(period)=-2
   dGdT=(G_1200 - G_900)/(1200-900)
   VBR.in.elastic.anharmonic.dG_dT=dGdT * 1e9; % Pa / C
   dGdT=VBR.in.elastic.anharmonic.dG_dT;

   % set JF10 ref modulus and ref T/P
   Gu0_x=G_900;
   T_ref_JF10=900+273;
   P_ref_JF10=0.2;

   % back out ref Modlus at STP.
   dGdT=VBR.in.elastic.anharmonic.dG_dT;
   dGdP=VBR.in.elastic.anharmonic.dG_dP;
   Tref=VBR.in.elastic.anharmonic.T_K_ref;
   Pref=VBR.in.elastic.anharmonic.P_Pa_ref/1e9;
   VBR.in.elastic.anharmonic.Gu_0_ol =  Gu0_x - (T_ref_JF10-Tref) * dGdT/1e9 ...
                 - (P_ref_JF10-Pref)*dGdP; % olivine reference shear modulus [GPa]


%  frequencies to calculate at
fmin = 0.001 ; %
fmax = 2.0 ; %
f_vec = logspace(log10(fmin),log10(fmax),20) ;
VBR.in.SV.f = f_vec ;
% FJax: f = [ 0.0056000 0.0100000 0.017800 0.031600 0.056200 0.10000 0.17780 0.31600 0.56230 1.0000 ]
%% ====================================================
%% define variable vectors ============================
%% ====================================================
% for LABORATORY CONDITIONS !

%T_C_vec = 800:50:1500 ;
T_C_vec = 900:100:1300 ;
gs_um_vec = linspace(1,20,20) ;
P_GPa_vec = [1e-4,0.100,0.200,0.300,0.400] ; %linspace(1e-4,400,40) ;
% room pressure (STP) is 101 kPa = 0.1 MPa = 1e-4 GPa
%phi_vec = linspace(0,0.05,25) ;

VBR.in.SV_vectors.T_K_vec_dim1 = T_C_vec+273 ;
VBR.in.SV_vectors.gs_um_vec_dim2 = gs_um_vec ;
VBR.in.SV_vectors.P_GPa_vec_dim3 = P_GPa_vec ;
%VBR.in.SV_vectors.phi_vec_dim3 = phi_vec ;

%[T_K_ra,gs_um_ra] = ndgrid(T_C_vec+273, gs_um_vec) ;
[T_K_ra,gs_um_ra,P_GPa_ra] = ndgrid(T_C_vec+273,gs_um_vec,P_GPa_vec) ;
%T_K_ra = T_K_ra'
oneses = ones(size(T_K_ra)) ; %,len(gs_um_vec),len(phi_vec));
sz=size(oneses) ; %

%% ====================================================
%% Define the Thermodynamic State ARRAYS===============
%% ====================================================

% indexes:
% 1 = temperature
% 2 = grain size
% 3 = pressure
% 4 = frequency
% (melt fraction, grain size, etc)

%  size of the state variable arrays. arrays can be any shape
%  but all arays must be the same shape.
VBR.in.SV.T_K= T_K_ra ; % VBR.in.SV.T_C+273; % temperature [K]
VBR.in.SV.dg_um = gs_um_ra;

%  intensive state variables (ISV)
VBR.in.SV.P_GPa = P_GPa_ra; % pressure [GPa]
VBR.in.SV.rho = 3300 * oneses; % density [kg m^-3]
VBR.in.SV.sig_MPa = 0.1 * oneses; % differential stress [MPa]

% asthenosphere conditions
   %VBR.in.SV.dg_um = 5e3*oneses;
   % VBR.in.SV.P_GPa = 0.3 * oneses; % pressure [GPa]
   % VBR.in.SV.rho = 3300 * oneses; % density [kg m^-3]
  %VBR.in.SV.sig_MPa = 0.5 * oneses; % differential stress [MPa]

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

[VBR] = VBR_spine(VBR)  ;

% for the python version ! (because it will not read the reserved word 'in')
VBR.input = VBR.in ;

save(filename,'VBR')


return
%% ====================================================
%% Display some things ================================
%% ====================================================
i1_T = [8:13] ;
i2_gs = 4 ;
disp('sigma (not varied) [MPa] = ');
disp(string(VBR.in.SV.sig_MPa(1,1,1)));
% disp('Temp [C] = ');
% disp(string(VBR.in.SV_vectors.T_K_vec_dim1(i1_T)-273));
% disp(string(VBR.in.SV.T_K(i1_T,1,1)-273));
disp('grain size (um) = ');
disp(string(VBR.in.SV_vectors.gs_um_vec_dim2(i2_gs )));
disp(string(VBR.in.SV.dg_um(1,i2_gs,1)));

% ==========================
close all;
figure;
subplot(1,3,1)
%semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT2016_solidus.M(1,:,:)/1e9), 'k'); hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT_maxwell.M(i1_T,i2_gs,:)/1e9), 'r');hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.AndradePsP.Ma(i1_T,i2_gs,:)/1e9), 'b' );hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.eBurgers.M(i1_T,i2_gs,:)/1e9), 'g' );hold on;
ylabel('M [GPa]'); xlabel('period [s]'); hold on;
title(['gs [um] = ', string(VBR.in.SV.dg_um(1,i2_gs,1)), 'stress [MPa]=', string(VBR.in.SV.sig_MPa(1,1,1)) ]); hold on;
text(1e1,10,'red=YTmxw')
text(1e1,6,'blue=AndradePsP')
text(1e1,2,'grn=eBurgers(PsP)')
ylim([0,80])
axis('tight')

subplot(1,3,2)
%loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT2016_solidus.Qinv(1,:,:)), 'k');hold on;
loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT_maxwell.Qinv(i1_T,i2_gs,:)), 'r');hold on;
loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.AndradePsP.Qinv(i1_T,i2_gs,:)), 'b');hold on;
loglog(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.eBurgers.Qinv(i1_T,i2_gs,:)), 'g');hold on;
ylabel('Q^-1'); xlabel('period [s]')
ylim([1e-3,.1])
axis('tight')

subplot(1,3,3)
%semilogx(1./VBR.in.SV.f,1e-3*squeeze(VBR.out.anelastic.YT2016_solidus.V(1,:,:)), 'k');hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.YT_maxwell.V(i1_T,i2_gs,:))./1e3, 'r');hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.AndradePsP.Va(i1_T,i2_gs,:)./1e3), 'b');hold on;
semilogx(1./VBR.in.SV.f,squeeze(VBR.out.anelastic.eBurgers.V(i1_T,i2_gs,:)./1e3), 'g');hold on;
ylabel('V_s [km/s]'); xlabel('period [s]')
title(['T [C] = ',VBR.in.SV_vectors.T_K_vec_dim1(min(i1_T))-273, ':', string(VBR.in.SV_vectors.T_K_vec_dim1(max(i1_T))-273) ])
axis('tight')
