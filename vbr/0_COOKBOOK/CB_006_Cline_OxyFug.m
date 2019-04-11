%% ===================================================================== %%
%%                     CB_006_Cline_OxyFug.m
%% ===================================================================== %%
%  WILL NOT WORK: was for an experimental version that used Cline's
%  maxwell-time oxygen fugacity dependence. See
%       vbr/4_VBR/VBR_v0p95/functions/addOxyFugacityEffects.m
%% ===================================================================== %%
   clear

%% ====================================================
%% Load and set VBR parameters ========================
%% ====================================================

% put VBR in the path
  path_to_top_level_vbr='../../';
  addpath(path_to_top_level_vbr)
  vbr_init

% write method list (these are the things to calculate)
  VBR.in.elastic.methods_list={'anharmonic'};
  VBR.in.viscous.methods_list={'LH2012'};
  VBR.in.anelastic.methods_list={'eBurgers','YT_maxwell'};

% frequencies to calculate at
  VBR.in.SV.f = 0.01;

%% ====================================================
%% Define the Thermodynamic State =====================
%% ====================================================

% oxygen fugacity
  VBR.in.SV.fO2_bar=logspace(-3,-0.4,30);
  VBR.in.anelastic.eBurgers=Params_Anelastic('eBurgers');
  VBR.in.anelastic.YT_maxwell=Params_Anelastic('YT_maxwell');
  VBR.in.elastic.anharmonic=Params_Elastic('anharmonic');

  VBR.in.anelastic.eBurgers.useJF10visc=0;
  VBR.in.GlobalSettings.melt_enhacement=0; % turn off critcal melt fraction effect

  dGdT=VBR.in.elastic.anharmonic.dG_dT;
  dGdP=VBR.in.elastic.anharmonic.dG_dP;
  Tref=VBR.in.elastic.anharmonic.T_K_ref;
  Pref=VBR.in.elastic.anharmonic.P_Pa_ref/1e9;
  VBR.in.elastic.anharmonic.Gu_0_ol = 66.5 - (900+273-Tref) * dGdT/1e9 - (0.2-Pref)*dGdP; % olivine reference shear modulus [GPa]

% size of the state variable arrays. arrays can be any shape
% but all arays must be the same shape.
  sz=size(VBR.in.SV.fO2_bar);

% intensive state variables (ISV)
  VBR.in.SV.dg_um=25 * ones(sz); % grain size [um]
  VBR.in.SV.phi=zeros(sz); % melt fraction
  VBR.in.SV.T_K=(1200+273) * ones(sz); % temperature [K]
  VBR.in.SV.P_GPa = 0.2 * ones(sz); % pressure [GPa]
  VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
  VBR.in.SV.sig_MPa = 10 * ones(sz); % differential stress [MPa]
  VBR.in.SV.chi=ones(sz); % composition fraction  1 = olivine, 0 = crust

%% ====================================================
%% CALL THE VBR CALCULATOR ============================
%% ====================================================

  [VBR] = VBR_spine(VBR) ;

  VBR_bai=VBR;
  VBR_bai.in.anelastic.eBurgers.m_fO2=-0.4; % Bai et al
  VBR_bai.in.anelastic.YT_maxwell.m_fO2=-0.4; % Bai et al
  [VBR_bai] = VBR_spine(VBR_bai) ;

%% ====================================================
%% Display some things ================================
%% ====================================================

close all;
figure;
plot(log10(VBR.in.SV.fO2_bar),log10(VBR.out.anelastic.eBurgers.tau_M),'k','DisplayName','m=-1.2','LineWidth',2);
xlabel('log f_O_2 [bar]','FontSize', 16)
ylabel('log tau_M [s]','FontSize', 16)

if exist('../../Data/ClineEtAl2018/cline_fig4.csv')
  Cline=csvread('../../Data/ClineEtAl2018/cline_fig4.csv');
  hold on
  plot(Cline(:,1),Cline(:,2),'.k','MarkerSize',12)
end
xlim([-4,0])
ylim([4,9])

hold on
plot(log10(VBR_bai.in.SV.fO2_bar),log10(VBR_bai.out.anelastic.eBurgers.tau_M),'--k','DisplayName','m=-0.4','LineWidth',2);
xlabel('log fO_2 [bar]','FontSize', 16)
ylabel('log tau_M [s]','FontSize', 16)
title('1200^oC, 25 {\mu}m','FontSize', 16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
set(gca,'XTickLabelMode','auto')
lgd=legend('Location','NorthEast');

figure;
subplot(1,3,1)
plot(log10(VBR.in.SV.fO2_bar),log10(VBR.out.anelastic.eBurgers.Qinv),'k','DisplayName','m=-1.2, eBurgers','LineWidth',2);
hold on
plot(log10(VBR_bai.in.SV.fO2_bar),log10(VBR_bai.out.anelastic.eBurgers.Qinv),'--k','DisplayName','m=-0.4, eBurgers','LineWidth',2);
plot(log10(VBR.in.SV.fO2_bar),log10(VBR.out.anelastic.YT_maxwell.Qinv),'b','DisplayName','m=-1.2, maxwell','LineWidth',2);
plot(log10(VBR_bai.in.SV.fO2_bar),log10(VBR_bai.out.anelastic.YT_maxwell.Qinv),'--b','DisplayName','m=-0.4, maxwell','LineWidth',2);
xlabel('log fO_2 [bar]','FontSize', 16)
ylabel('log Q^{-1}','FontSize', 16)
xlim([-4,0])
title('1200^oC, 100s, 25 {\mu}m','FontSize', 16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
set(gca,'XTickLabelMode','auto')

subplot(1,3,2)
set(gca,'FontSize',16)
plot(log10(VBR.in.SV.fO2_bar),log10(VBR.out.anelastic.eBurgers.Vave),'k','DisplayName','m=-1.2, eBurgers','LineWidth',2);
hold on
plot(log10(VBR_bai.in.SV.fO2_bar),log10(VBR_bai.out.anelastic.eBurgers.Vave),'--k','DisplayName','m=-0.4, eBurgers','LineWidth',2);
plot(log10(VBR.in.SV.fO2_bar),log10(VBR.out.anelastic.YT_maxwell.Vave),'b','DisplayName','m=-1.2, maxwell','LineWidth',2);
plot(log10(VBR_bai.in.SV.fO2_bar),log10(VBR_bai.out.anelastic.YT_maxwell.Vave),'--b','DisplayName','m=-0.4, maxwell','LineWidth',2);
xlabel('log fO_2 [bar]','FontSize', 16)
ylabel('V^s [km/s]','FontSize', 16)
xlim([-4,0])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
set(gca,'XTickLabelMode','auto')

hSub = subplot(1,3,3);
plot([1 2],[nan nan],'k','DisplayName','m=-1.2, eBurgers','LineWidth',2);
hold on
plot([1 2],[nan nan],'--k','DisplayName','m=-0.4, eBurgers','LineWidth',2);
plot([1 2],[nan nan],'b','DisplayName','m=-1.2, maxwell','LineWidth',2);
plot([1 2],[nan nan],'--b','DisplayName','m=-0.4, maxwell','LineWidth',2);
set(hSub, 'Visible', 'off');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
set(gca,'XTickLabelMode','auto')
legend(hSub, 'Location', 'north');
