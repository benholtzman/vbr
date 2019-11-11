% PLOT Experimental data from
% Faul and Jackson, 2015 (Ann. Rev.), compilation of other data datasets
% and test fit to eBurgersPsP (and AndradePsP)
clear all ; clf;

% LOAD THE DATA
if ~exist('ExptData.mat')
  Make_DATA ;
end
load('ExptData.mat');
data = Data ;

% ===========================================
% LOAD THE VBR or Run it here... not as LUT.

runVBRwhere='here' % 'LUT'

if strcmp(runVBRwhere,'LUT')==1
  load('VBR_LUT_labdata');
end

if strcmp(runVBRwhere,'here')==1
for i=1:length(data.FaulJax15)
  T_C_vec(i) = data.FaulJax15(i).exptCond.T_C  ;
end
T_C_vec ;

VBR.in.elastic.methods_list={'anharmonic'};
VBR.in.viscous.methods_list={'HK2003';'LH2011'};
VBR.in.anelastic.methods_list={'eBurgers';'AndradePsP'};
VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
VBR.in.anelastic.eBurgers.eBurgerMethod='bg_peak'; % 'bg_only' or 'bg_peak'
VBR.in.GlobalSettings.melt_enhacement = 0 ;


% ===================================================
% rescale the reference modulus =====================

dGdT=VBR.in.elastic.anharmonic.dG_dT;
dGdP=VBR.in.elastic.anharmonic.dG_dP;
Tref=VBR.in.elastic.anharmonic.T_K_ref;
Pref=VBR.in.elastic.anharmonic.P_Pa_ref/1e9;


% % method from CB_003_JF10:
% % JF10 have Gu_0=62.5 GPa, but that's at 900 Kelvin and 0.2 GPa,
% % so set Gu_0_ol s.t. it ends up at 62.5 at those conditions
% Gu_ref = data.FaulJax15(1).exptCond.Gu
% Gu_0_ol = Gu_ref - (900+273-Tref) * dGdT/1e9 - (0.2-Pref)*dGdP
% VBR.in.elastic.anharmonic.Gu_0_ol = Gu_0_ol ;   % olivine reference shear modulus [GPa]

% method from Make_LUT:
% JF10 have Gu_0=66.5 GPa, but that's at 900 C and 0.2 GPa,
% so set Gu_0_ol s.t. it ends up at 66.5 at those conditions
% calculate dGdT from their figure.
G_900=66.5; % from plot at log10(period)=-2
G_1200=58; % from plot at log10(period)=-2
dGdT=(G_1200 - G_900)/(1200-900);
VBR.in.elastic.anharmonic.dG_dT=dGdT * 1e9; % Pa / C
dGdT=VBR.in.elastic.anharmonic.dG_dT;

% set JF10 ref modulus and ref T/P
Gu0_x=G_900;
T_ref_JF10=900+273;
P_ref_JF10=0.2;
% back out ref Modlus at STP.
Gu_0_ol =  Gu0_x - (T_ref_JF10-Tref) * dGdT/1e9 - (P_ref_JF10-Pref)*dGdP
% olivine reference shear modulus [GPa]
VBR.in.elastic.anharmonic.Gu_0_ol = Gu_0_ol ;

% ==================================================
%  frequencies to calculate at
VBR.in.SV.f = 1./logspace(0,3,30);
% [0.7804 0.2616 0.1561 0.0863 0.0456 0.0214 0.0100 0.0048 0.0021 0.0010]

%  size of the state variable arrays. arrays can be any shape
%  but all arays must be the same shape.
VBR.in.SV.T_K = T_C_vec+273 ;
VBR.in.SV_vectors.T_K_vec_dim1 = VBR.in.SV.T_K ;
sz=size(VBR.in.SV.T_K) ; % temperature [K]

%  remaining state variables (ISV)
VBR.in.SV.dg_um= data.FaulJax15(1).exptCond.dg_0 .* ones(sz);
VBR.in.SV.P_GPa = data.FaulJax15(1).exptCond.P_GPa .* ones(sz); % pressure [GPa]
VBR.in.SV.rho = data.FaulJax15(1).exptCond.rho .* ones(sz); % density [kg m^-3]
VBR.in.SV.sig_MPa = (data.FaulJax15(1).exptCond.sig_0 .* ones(sz))./1e6; % differential stress [MPa]
VBR.in.SV.phi = data.FaulJax15(1).exptCond.phi_0 .* ones(sz); % melt fraction

% run VBR
[VBR] = VBR_spine(VBR) ;
end

% % adjust VBR input and get out eBurgers with background + peak
% VBR.in.anelastic.eBurgers=Params_Anelastic('eBurgers');
% VBR.in.anelastic.eBurgers.eBurgerMethod='bg_peak';
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
nlines = length(data.FaulJax15) ; %length(VBR_sols(1).T_params) ;
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

  if strcmp(runVBRwhere,'LUT')==1
    state = data.FaulJax15(iT).exptCond ;
    [i_T_d1, i_g_d2, i_P_d3] = find_index_f(VBR,state) ;
    Qs = VBR.out.anelastic.eBurgers.Q(i_T_d1, i_g_d2, i_P_d3,:) ;
    Q = squeeze(Qs) ;
  elseif strcmp(runVBRwhere,'here')==1
    Q(:) = squeeze(VBR.out.anelastic.eBurgers.Q(1,iT,:)) ;
  end

  if plot_vs_freq
    plot(log10(f_vec),log10(1./Q),'k-','LineWidth', LineW, 'Color', clr); hold on;
  else
    plot(log10(1./f_vec),log10(1./Q),'k-','LineWidth', LineW, 'Color', clr); hold on;
  end

  % PLOT DATA
  data_log10_Qinv = data.FaulJax15(iT).Results.log10_Qinv ;
  if plot_vs_freq
    data_freq = data.FaulJax15(iT).exptCond.f ;
    plot(log10(data_freq),data_log10_Qinv,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
  else
    data_logPer = data.FaulJax15(iT).exptCond.logPer
    plot(data_logPer,data_log10_Qinv,'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
  end

end

axis tight
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])
title(['Extended Burgers'],'fontname','Times New Roman','fontsize',LBLFNT);
xlabel(xlabel_text, 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('log_{10} Q^{-1}, attenuation', 'fontname','Times New Roman','fontsize', LBLFNT)
%ylabel('log_{10} Q^{-1}, (J_1/J_2)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);



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
        % G = VBR.out.anelastic.eBurgers.M(i_T_d1, i_g_d2, i_P_d3,:) ;
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
    Ms = VBR.out.anelastic.eBurgers.M(i_T_d1, i_g_d2, i_P_d3,:)./1e9 ;
    M = squeeze(Ms) ;
  elseif strcmp(runVBRwhere,'here')==1
    M(:) = squeeze(VBR.out.anelastic.eBurgers.M(1,iT,:)./1e9) ;
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
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])

xlabel(xlabel_text, 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('Modulus, M (GPa)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);
%set(gca,'XTickLabel', [])



return
%
% %%  Q vs FREQUENCY (AndradePsP) ==================================================
% axes('Position', plot_row1_B);
%
% f = VBR_sols(1).VBR.ISV.f ;
% % AndradePsP, two frequencies
% for j = 1:nlines
%         %LineW = LineW_vec(j);
%         clr = colorscale(j,:) ;
%        % f_M = (VBR_sols(j).VBR.Gu_vec(1))/(VBR_sols(j).VBR.eta_total(1)) ;
%        % I_fM = find(f>=f_M,1);
%         Q = VBR_sols(j).VBR.AndradePsP.Qa ;
%         Qcomp = VBR_sols(j).VBR.AndradePsP.Q_comp ;
%        % MAT(:,j) = Qs ;
%
%         plot(log10(f),log10(Q),'k--','LineWidth', LineW-1, 'Color', clr); hold on;
%         plot(log10(f),log10(Qcomp),'k-','LineWidth', LineW, 'Color', clr); hold on;
%         %plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;
%
%         % PLOT DATA:
%         %plot(data(j).logf,data(j).logQ,'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
%
% end
%
% axis tight
% %xlim([1.8e1 3e2])
% %ylim([1e-6 5e-4])
% title(['AndradePsP, P=' num2str(P) ' GPa, d=' num2str(d) ' \mum'],'fontname','Times New Roman','fontsize',LBLFNT);
% xlabel('log frequency', 'fontname','Times New Roman','fontsize', LBLFNT)
% ylabel('log Q, (J_1/J_2)', 'fontname','Times New Roman','fontsize', LBLFNT)
% set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
% set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);
%
%
%
%
% %%  G vs FREQUENCY (AndradePsP) ==================================================
% axes('Position', plot_row2_D);
%
% f = VBR_sols(1).VBR.ISV.f ;
% % AndradePsP, two frequencies
% for j = 1:lines
%         %LineW = LineW_vec(j);
%         clr = colorscale(j,:) ;
%         %f_M = (VBR_sols(j).VBR.Gu_vec(1))/(VBR_sols(j).VBR.eta_total(1)) ;
%         %I_fM = find(f>=f_M,1);
%         Ma = VBR_sols(j).VBR.AndradePsP.Ma ;
%         Mcomp = VBR_sols(j).VBR.AndradePsP.Mcomp ;
%         %MAT(:,j) = Qs ;
%
%         plot(log10(f),Ma./1e9,'k--','LineWidth', LineW-1, 'Color', clr); hold on;
%         plot(log10(f),Mcomp./1e9,'k-','LineWidth', LineW, 'Color', clr); hold on;
%         %plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;
%
%         % PLOT DATA:
%         %plot(data(j).logf,data(j).G,'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
%
% end
%
% axis tight
% %xlim([1.8e1 3e2])
% %ylim([1e-6 5e-4])
%
% xlabel('log frequency', 'fontname','Times New Roman','fontsize', LBLFNT)
% ylabel('Modulus, M (GPa)', 'fontname','Times New Roman','fontsize', LBLFNT)
% set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
% set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);
% %set(gca,'XTickLabel', [])
