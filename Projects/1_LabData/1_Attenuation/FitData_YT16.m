% PLOT Experimental data from
% Faul and Jackson, 2015 (Ann. Rev.), compilation of other data datasets
% and test fit to eburgers_psp (and andrade_psp)
% clear all ; clf;
clear
% put VBR in the path
  path_to_top_level_vbr='../../../';
  addpath(path_to_top_level_vbr)
  vbr_init
  addpath('./functions')
% LOAD THE DATA
if ~exist('ExptData.mat')
  Make_DATA ;
end
load('ExptData.mat');
data = Data ;

% ===========================================
% LOAD THE VBR or Run it here... not as LUT.

% runVBRwhere='here' % 'LUT'
%
% if strcmp(runVBRwhere,'LUT')==1
%   load('VBR_LUT_labdata');
% end

% if strcmp(runVBRwhere,'here')==1
for i=1:length(data.YT16)
  T_C_vec(i) = data.YT16(i).exptCond.T_C  ;
end
T_C_vec



VBR.in.elastic.methods_list={'anharmonic';'anh_poro'};
VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity


[E,dEdT,dEdT_ave]= YT16_E(T_C_vec);
[E_o,~]=YT16_E(0);
poifac=(2 *( 1 + VBR.in.elastic.anharmonic.nu));
Gu_o=E_o / poifac;
dGdT_ave=dEdT_ave /poifac;

VBR.in.elastic.anharmonic.Gu_0_ol = Gu_o;%2.5 ; %2.6 estimated from Fig 6 (YT16)
VBR.in.elastic.anharmonic.dG_dT = dGdT_ave;%-5.0e4 ; % -13600000 % Below about -5e4, this makes no difference..
% so it must be anelastic.
%VBR.in.elastic.anharmonic.dG_dP = 1.800 ;
%VBR.in.elastic.anharmonic.T_K_ref = 25 ;

VBR.in.viscous.methods_list={'xfit_premelt'}; %VBR.in.viscous.methods_list={'xfit_premelt'};
VBR.in.viscous.xfit_premelt=setBorneolParams();
% VBR.in.viscous.xfit_premelt.H=
VBR.in.anelastic.methods_list={'xfit_premelt'};
% VBR.in.anelastic.eburgers_psp=Params_Anelastic('eburgers_psp');
% fit_type='bg_peak';
% VBR.in.anelastic.eburgers_psp.eBurgerFit=fit_type; % 'bg_only' or 'bg_peak'



% ==================================================
%  frequencies to calculate at
f_data = data.YT16(1).exptCond.f ;

VBR.in.SV.f = logspace(log10(min(f_data)),log10(max(f_data)),30);

%  size of the state variable arrays. arrays can be any shape
%  but all arays must be the same shape.
VBR.in.SV.T_K = T_C_vec+273 ;
VBR.in.SV_vectors.T_K_vec_dim1 = VBR.in.SV.T_K ;
sz=size(VBR.in.SV.T_K) ; % temperature [K]

%  remaining state variables (ISV)
VBR.in.SV.dg_um= data.YT16(1).exptCond.dg_0 .* ones(sz);
VBR.in.SV.P_GPa = data.YT16(1).exptCond.P_GPa .* ones(sz); % pressure [GPa]
VBR.in.SV.rho = data.YT16(1).exptCond.rho .* ones(sz); % density [kg m^-3]
VBR.in.SV.sig_MPa = (data.YT16(1).exptCond.sig_0 .* ones(sz))./1e6; % differential stress [MPa]
VBR.in.SV.phi = data.YT16(1).exptCond.phi_0 .* ones(sz); % melt fraction
VBR.in.SV.Tsolidus_K = 43.0 + 273 ;

% run VBR
[VBR] = VBR_spine(VBR) ;
% end

% ===================================
% PLOTTING
% ===================================

% left bottom width height
figure()
W = 0.33 ;
H = 0.36 ;
plot_row1_A = [0.1 0.60 W H ] ;
plot_row1_B = [0.55 0.60 W H ] ;
plot_row2_C = [0.1 0.12 W H ] ;
plot_row2_D = [0.55 0.12 W H ] ;


LBLFNT = 14 ;
LineW = 2 ;
%LineW_vec = linspace(1,3,nlines);
dotsize = 12;
dotsize_D = 14 ;


plot_vs_freq=1;
if plot_vs_freq
  xlabel_text='log_{10} frequency';
else
  xlabel_text='log_{10} period';
end

%% PLOT =======================================================
nlines = length(data.YT16) ; %length(VBR_sols(1).T_params) ;
%cool to warm:
colorscale(:,1) = linspace(0.5,1,nlines) ;
colorscale(:,2) = linspace(0,0,nlines) ;
colorscale(:,3) = linspace(1,0,nlines) ;

%%  ==================================================
%%  finding and PLOTTING !
%%  ==================================================
f_vec = VBR.in.SV.f ;

%%  Q vs FREQUENCY (BURGERS) ==================================================
axes('Position', plot_row1_A);

for iT = 1:nlines
    %LineW = LineW_vec(j);
    clr = colorscale(iT,:) ;

    % Plot model
    Q = squeeze(VBR.out.anelastic.xfit_premelt.Q(1,iT,:)) ;
    if plot_vs_freq
      loglog((f_vec),(1./Q),'k-','LineWidth', LineW, 'Color', clr); hold on;
    else
      plot(log10(1./f_vec),log10(1./Q),'k-','LineWidth', LineW, 'Color', clr); hold on;
    end

  % PLOT DATA
  data_log10_Qinv = data.YT16(iT).Results.log10_Qinv ;
  if plot_vs_freq
    data_freq = data.YT16(iT).exptCond.f ;
    loglog((data_freq),10.^data_log10_Qinv,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
  else
    data_logPer = data.YT16(iT).exptCond.logPer
    plot(data_logPer,data_log10_Qinv,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
  end

end

axis tight
% xlim([1.8e1 3e2])
%ylim([-2.1,.5])
title(['Yamauchi+Takei 2016 data (sample 40), using xfit\_premelt '],'fontname','Times New Roman','fontsize',LBLFNT);
xlabel(xlabel_text, 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('log_{10} Q^{-1}, attenuation', 'fontname','Times New Roman','fontsize', LBLFNT)
%ylabel('log_{10} Q^{-1}, (J_1/J_2)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);


%%  G vs FREQUENCY (BURGERS) ==================================================
axes('Position', plot_row1_B);

for iT = 1:nlines

  %LineW = LineW_vec(j);
  clr = colorscale(iT,:) ;

  % Plot model
  M = poifac * squeeze(VBR.out.anelastic.xfit_premelt.M(1,iT,:)./1e9) ;
  if plot_vs_freq
    semilogx((f_vec),M,'k-','LineWidth', LineW, 'Color', clr); hold on;
  else
    plot(log10(1./f_vec),M,'k-','LineWidth', LineW, 'Color', clr); hold on;
  end

  % PLOT DATA
  data_M = data.YT16(iT).Results.G ;
  if plot_vs_freq
    data_freq = data.YT16(iT).exptCond.f ;
    semilogx((data_freq),data_M,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
  else
    data_logPer = data.YT16(iT).exptCond.logPer
    plot(data_logPer,data_M,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
  end


end

axis tight
%ylim([0,80])
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])

xlabel(xlabel_text, 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('Modulus, M (GPa)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);
%set(gca,'XTickLabel', [])


% % ===================
% return
% % ===================
