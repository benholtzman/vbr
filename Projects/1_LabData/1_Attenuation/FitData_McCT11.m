% PLOT Experimental data from McCarthy and Takei, 2011
% and fit with xfit_mxw (empirical fit to data, scaled by the maxwell time for diffusion creep)

  clear all ; close all;

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
  % Run VBR
  for i=1:length(data.McCT11)
    d_vec(i) = data.McCT11(i).exptCond.dg  ;
  end
  d_vec


  VBR.in.elastic.methods_list={'anharmonic'};
  VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
  VBR.in.elastic.anharmonic.Gu_0_ol = 2.8 ; %2.6 estimated from Fig 6 (YT16)
  VBR.in.elastic.anharmonic.dG_dT = 0 ; % -13600000 % Below about -5e4, this makes no difference..
  % so it must be anelastic.
  VBR.in.elastic.anharmonic.dG_dP = 0;%1.800 ;
  %VBR.in.elastic.anharmonic.T_K_ref = 25 ;

  VBR.in.viscous.methods_list={'xfit_premelt'}; %VBR.in.viscous.methods_list={'xfit_premelt'};
  % small grain size viscosity parameters:
  VBR.in.viscous.xfit_premelt.Tr_K=23.6+273; %
  VBR.in.viscous.xfit_premelt.Pr_Pa=1000; % p7817 of YT2016, second paragraph
  VBR.in.viscous.xfit_premelt.eta_r=2.04*1e12; %
  VBR.in.viscous.xfit_premelt.H=85.4*1e3; % activation energy [J/mol], figure 20 of YT2016
  VBR.in.viscous.xfit_premelt.V=0; % activation vol [m3/mol], figure 20 of YT2016
  VBR.in.viscous.xfit_premelt.R=8.314; % gas constant [J/mol/K]
  VBR.in.viscous.xfit_premelt.m=3; % grain size exponent
  VBR.in.viscous.xfit_premelt.dg_um_r=3.35 ; % caption of Fig 9. % 24.4; % reference grain size [um]

  VBR.in.anelastic.methods_list={'xfit_mxw'};

  VBR.in.GlobalSettings.melt_enhacement = 0 ;

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



  %%  G vs FREQUENCY (BURGERS) ==================================================
  axes('Position', plot_row1_B);


  for iT = 1:nlines
      %LineW = LineW_vec(j);
    clr = colorscale(iT,:) ;

    E = squeeze(VBR.out.anelastic.xfit_mxw.M(1,iT,:)) ;
    if plot_vs_freq
      plot(log10(f_vec),E./1e9,'k-','LineWidth', LineW, 'Color', clr); hold on;
    else
      plot(log10(1./f_vec),E./1e9,'k-','LineWidth', LineW, 'Color', clr); hold on;
    end

    % PLOT DATA
    data_E = data.McCT11(iT).Results.E ;
    if plot_vs_freq
      data_freq = data.McCT11(iT).exptCond.f ;
      plot(log10(data_freq),data_E,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
    else
      data_logPer = data.McCT11(iT).exptCond.logPer
      plot(data_logPer,data_E,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
    end


  end

  axis tight
  %ylim([0,80])
  %xlim([1.8e1 3e2])
  %ylim([1e-6 5e-4])

  xlabel(xlabel_text, 'fontname','Times New Roman','fontsize', LBLFNT)
  ylabel('Youngs Modulus,E (GPa)', 'fontname','Times New Roman','fontsize', LBLFNT)
  set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
  set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);
  %set(gca,'XTickLabel', [])


% function params = SetBorneolParams()
%   % set the viscous parameters for borneol
%   % near-solidus and melt effects
%   params.alpha=0;
%   params.T_eta=0.94; % eqn 17,18- T at which homologous T for premelting.
%   params.gamma=5;
%   % flow law constants for YT2016
%   params.Tr_K=8+273; % p7817 of YT2016, second paragraph
%   params.Pr_Pa=1000; % p7817 of YT2016, second paragraph
%   params.eta_r=565e12; % figure 20 of reference paper Table 2 @ 8 degrees
%   params.H=141*1e3; % activation energy [J/mol], figure 20 of YT2016
%   params.V=0; % activation vol [m3/mol], figure 20 of YT2016
%   params.R=8.314; % gas constant [J/mol/K]
%   params.m=2.56; % grain size exponent
%   params.dg_um_r=4.3 ; % caption of Fig 9. % 24.4; % reference grain size [um]
% end
