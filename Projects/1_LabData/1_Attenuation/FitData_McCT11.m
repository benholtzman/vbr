% PLOT Experimental data from
% Faul and Jackson, 2015 (Ann. Rev.), compilation of other data datasets
% and test fit to eburgers_psp (and andrade_psp)
clear all ; clf;

% put VBR in the path
  path_to_top_level_vbr='../../../';
  addpath(path_to_top_level_vbr)
  vbr_init

% LOAD THE DATA
if ~exist('ExptData.mat')
  Make_DATA ;
end
load('ExptData.mat');
data = Data ;

% ===========================================
% Run VBR

% if strcmp(runVBRwhere,'here')==1
for i=1:length(data.McCT11)
  d_vec(i) = data.McCT11(i).exptCond.dg  ;
end
d_vec


VBR.in.elastic.methods_list={'anharmonic'};
VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
VBR.in.elastic.anharmonic.Gu_0_ol = 2.5 ; %2.6 estimated from Fig 6 (YT16)
VBR.in.elastic.anharmonic.dG_dT = -5.0e4 ; % -13600000 % Below about -5e4, this makes no difference..
% so it must be anelastic.
%VBR.in.elastic.anharmonic.dG_dP = 1.800 ;
%VBR.in.elastic.anharmonic.T_K_ref = 25 ;

VBR.in.viscous.methods_list={'xfit_premelt'}; %VBR.in.viscous.methods_list={'xfit_premelt'};
VBR.in.viscous.xfit_premelt=SetBorneolParams();
VBR.in.viscous.xfit_premelt.eta_r = 150e12 ;  % 565e12 from figure 20 of reference paper Table 2 @ 8 degrees
% note that the Q_xfit_mxw takes the first viscosity method in the list, so use only that which you want it to use!

VBR.in.anelastic.methods_list={'xfit_mxw'};
% VBR.in.anelastic.eburgers_psp=Params_Anelastic('eburgers_psp');
% fit_type='bg_peak';
% VBR.in.anelastic.eburgers_psp.eBurgerFit=fit_type; % 'bg_only' or 'bg_peak'

VBR.in.GlobalSettings.melt_enhacement = 0 ;

% VBR.in.anelastic.eburgers_psp=Params_Anelastic('eburgers_psp');
% fit_type='bg_peak';
% VBR.in.anelastic.eburgers_psp.eBurgerFit=fit_type; % 'bg_only' or 'bg_peak'

% ===================================================
% rescale the reference modulus =====================

% % pull out anharmonic scaling
% dGdT=VBR.in.elastic.anharmonic.dG_dT;
% dGdP=VBR.in.elastic.anharmonic.dG_dP;
% Tref=VBR.in.elastic.anharmonic.T_K_ref;
% Pref=VBR.in.elastic.anharmonic.P_Pa_ref/1e9;

% JF10 ref modulus is for their T/P (900C,.2GPa):
% Gu0_x=VBR.in.anelastic.eburgers_psp.(fit_type).G_UR;
% T_ref_JF10=VBR.in.anelastic.eburgers_psp.(fit_type).TR;
% P_ref_JF10=VBR.in.anelastic.eburgers_psp.(fit_type).PR;

% % back out ref Modlus at STP.
% Gu_0_ol =  Gu0_x - (T_ref_JF10-Tref) * dGdT/1e9 - (P_ref_JF10-Pref)*dGdP
% VBR.in.elastic.anharmonic.Gu_0_ol = Gu_0_ol ;% olivine reference shear modulus [GPa]

% ==================================================
%  frequencies to calculate at
f_data = data.McCT11(1).exptCond.f ;
log10f = log10(f_data) ;
VBR.in.SV.f = logspace(min(log10f),max(log10f),30);

%  size of the state variable arrays. arrays can be any shape
%  but all arays must be the same shape.

VBR.in.SV.dg_um= d_vec ; %data.McCT11(1).exptCond.dg_0 .* ones(sz) ;
VBR.in.SV_vectors.d_vec_dim1 = VBR.in.SV.dg_um ;
sz=size(VBR.in.SV.dg_um) ; % grain size..

%  remaining state variables (ISV)
VBR.in.SV.T_K = data.McCT11(1).exptCond.T_C +273 .* ones(sz); % pressure [GPa]
VBR.in.SV.P_GPa = data.McCT11(1).exptCond.P_GPa .* ones(sz); % pressure [GPa]
VBR.in.SV.rho = data.McCT11(1).exptCond.rho .* ones(sz); % density [kg m^-3]
VBR.in.SV.sig_MPa = (data.McCT11(1).exptCond.sig_0 .* ones(sz))./1e6; % differential stress [MPa]
VBR.in.SV.phi = data.McCT11(1).exptCond.phi_0 .* ones(sz); % melt fraction

VBR.in.SV.Tsolidus_K = 43.0 + 273 ;

% run VBR
[VBR] = VBR_spine(VBR) ;
% end

% % adjust VBR input and get out eburgers_psp with background + peak
% VBR.in.anelastic.eburgers_psp=Params_Anelastic('eburgers_psp');
% VBR.in.anelastic.eburgers_psp.eBurgerFit='bg_peak';
% % Gu_0_ol = 62.5 - (900+273-Tref) * dGdT/1e9 - (0.2-Pref)*dGdP ;
% VBR.in.elastic.anharmonic.Gu_0_ol = 66.5 - (900+273-Tref) * dGdT/1e9 - (0.2-Pref)*dGdP ;
%
% [VBR_with_peak] = VBR_spine(VBR) ;

% ===================================
% PLOTTING
% ===================================

% left bottom width height
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
dotsize_D = 20 ;


plot_vs_freq=1;
if plot_vs_freq
  xlabel_text='log_{10} frequency';
else
  xlabel_text='log_{10} period';
end

%% PLOT =======================================================
nlines = length(data.McCT11) ; %length(VBR_sols(1).T_params) ;
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

  Q = squeeze(VBR.out.anelastic.xfit_mxw.Q(1,iT,:)) ;
  if plot_vs_freq
    plot(log10(f_vec),log10(1./Q),'k-','LineWidth', LineW, 'Color', clr); hold on;
  else
    plot(log10(1./f_vec),log10(1./Q),'k-','LineWidth', LineW, 'Color', clr); hold on;
  end

  % PLOT DATA
  data_log10_Qinv = data.McCT11(iT).Results.log10_Qinv ;
  if plot_vs_freq
    data_freq = data.McCT11(iT).exptCond.f ;
    plot(log10(data_freq),data_log10_Qinv,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
  else
    data_logPer = data.McCT11(iT).exptCond.logPer
    plot(data_logPer,data_log10_Qinv,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
  end

end

axis tight
% xlim([1.8e1 3e2])
ylim([-2.1,.5])
title(['McCarthy et al., 2011 data; xfit-maxwell fit'],'fontname','Times New Roman','fontsize',LBLFNT);
xlabel(xlabel_text, 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('log_{10} Q^{-1}, attenuation', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);



% ===================
return
% ===================

%%  G vs FREQUENCY (BURGERS) ==================================================
axes('Position', plot_row1_B);


for iT = 1:nlines
        %LineW = LineW_vec(j);
        % clr = colorscale(iT,:) ;
        % % %f_M = (VBR_sols(j).VBR.Gu_vec(1))/(VBR_sols(j).VBR.eta_total(1)) ;
        % % %I_fM = find(f>=f_M,1);
        %
        % state = data.FaulJax15(iT).exptCond ;
        % [i_T_d1, i_g_d2, i_P_d3] = find_index_f(VBR,state) ;
        % G = VBR.out.anelastic.eburgers_psp.M(i_T_d1, i_g_d2, i_P_d3,:) ;
        % G = squeeze(G) ;
        %
        % if plot_vs_freq
        %   plot(log10(f_vec),G./1e9,'k-','LineWidth', LineW, 'Color', clr); hold on;
        % else
        %   plot(log10(1./f_vec), G./1e9, 'k-','LineWidth', LineW, 'Color', clr); hold on;
        % end
        % % plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;
        %
        % % PLOT DATA
        %
        % data_G = data.FaulJax15(iT).Results.G ;
        % if plot_vs_freq
        %   data_freq = data.FaulJax15(iT).exptCond.f ;
        %   plot(log10(data_freq),data_G,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
        % else
        %   logPer = data.FaulJax15(iT).exptCond.logPer ;
        %   plot(logPer,data_G,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
        % end
        % %plot(data.TanJax.exptCond.logf,1./(data.TanJax.Results.Qinv),'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;

    %LineW = LineW_vec(j);
    clr = colorscale(iT,:) ;

  if strcmp(runVBRwhere,'LUT')==1
    state = data.FaulJax15(iT).exptCond ;
    [i_T_d1, i_g_d2, i_P_d3] = find_index_f(VBR,state) ;
    Ms = VBR.out.anelastic.eburgers_psp.M(i_T_d1, i_g_d2, i_P_d3,:)./1e9 ;
    M = squeeze(Ms) ;
  elseif strcmp(runVBRwhere,'here')==1
    M = squeeze(VBR.out.anelastic.eburgers_psp.M(1,iT,:)./1e9) ;
  end

  if plot_vs_freq
    plot(log10(f_vec),M,'k-','LineWidth', LineW, 'Color', clr); hold on;
  else
    plot(log10(1./f_vec),M,'k-','LineWidth', LineW, 'Color', clr); hold on;
  end

  % PLOT DATA
  data_M = data.FaulJax15(iT).Results.G ;
  if plot_vs_freq
    data_freq = data.FaulJax15(iT).exptCond.f ;
    plot(log10(data_freq),data_M,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
  else
    data_logPer = data.FaulJax15(iT).exptCond.logPer
    plot(data_logPer,data_M,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
  end


end

axis tight
ylim([0,80])
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])

xlabel(xlabel_text, 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('Modulus, M (GPa)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);
%set(gca,'XTickLabel', [])


function params=SetBorneolParams()
  % set the viscous parameters for borneol
  % near-solidus and melt effects
  params.alpha=0;
  params.T_eta=0.94; % eqn 17,18- T at which homologous T for premelting.
  params.gamma=5;
  % flow law constants for YT2016
  params.Tr_K=8+273; % p7817 of YT2016, second paragraph
  params.Pr_Pa=1000; % p7817 of YT2016, second paragraph
  params.eta_r=565e12; % figure 20 of reference paper Table 2 @ 8 degrees
  params.H=141*1e3; % activation energy [J/mol], figure 20 of YT2016
  params.V=0; % activation vol [m3/mol], figure 20 of YT2016
  params.R=8.314; % gas constant [J/mol/K]
  params.m=2.56; % grain size exponent
  params.dg_um_r=4.3 ; % caption of Fig 9. % 24.4; % reference grain size [um]
end
