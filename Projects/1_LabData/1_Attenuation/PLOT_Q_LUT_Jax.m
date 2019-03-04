
addpath('./datasets') ;
%VBR_LUT = load('VBR_LUT_labdata') ;
VBR_LUT = load('VBR_LUT_labdata_y190304');
VBR = VBR_LUT.VBR ;
VBR.in.SV_vectors;

Make_DATA ;
data = Data;

%find_index_f ;

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


% VBR = VBR_sols(1).VBR ;
% d = VBR.SSV.dg_um(1) ;
% P = VBR.ISV.P_GPa(1) ;

%% PLOT =======================================================
% MAKE THE DATA LISTS
% [BH 2019: What was this for? ]
% for j = 1:nlines
%     data(j).T = VBR_sols(j).T_params ;
%     data(j).logf = Data.SundCoop(j).logf ;
%     data(j).logQ = log10(1./Data.SundCoop(j).Qinv) ;
%     data(j).G = Data.SundCoop(j).G ;
% end

%colorscale = [ 1 0.5 0.2 ; 1 0 0 ] ;
%cool to warm:
nlines = 1
colorscale(:,1) = linspace(0.5,1,nlines) ;
colorscale(:,2) = linspace(0,0,nlines) ;
colorscale(:,3) = linspace(1,0,nlines) ;

%%  function to find the VBR solutions at each expt condition =================



%%  ==================================================
%%  finding and PLOTTING !
%%  ==================================================

clf;
%%  Q vs FREQUENCY (BURGERS) ==================================================
axes('Position', plot_row1_A);

nlines = 1 ; %length(VBR_sols(1).T_params) ;
for j = 1:nlines
        %LineW = LineW_vec(j);
        clr = colorscale(j,:) ;
        % %f_M = (VBR_sols(j).VBR.Gu_vec(1))/(VBR_sols(j).VBR.eta_total(1)) ;
        % %I_fM = find(f>=f_M,1);

        state = data.TanJax.exptCond ;
        [i_T_d1, i_g_d2, i_P_d3] = find_index_f(VBR,state) ;

        Qs = VBR.out.anelastic.eBurgers.Q(i_T_d1, i_g_d2, i_P_d3,:) ;
        Q = squeeze(Qs) ;
        f_vec = VBR.in.SV.f ;
        plot(log10(f_vec),log10(1./Q),'k-','LineWidth', LineW, 'Color', clr); hold on;
        % plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;

        % PLOT DATA:
        plot(state.logf,log10((squeeze(data.TanJax.Results.Qinv))),'k.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
        %plot(data.TanJax.exptCond.logf,1./(data.TanJax.Results.Qinv),'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;

end

axis tight
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])
title(['Extended Burgers'],'fontname','Times New Roman','fontsize',LBLFNT);
xlabel('log_{10} frequency', 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('log_{10} Q^{-1}, attenuation', 'fontname','Times New Roman','fontsize', LBLFNT)
%ylabel('log_{10} Q^{-1}, (J_1/J_2)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);






%
%
% %%  G vs FREQUENCY (BURGERS) ==================================================
% axes('Position', plot_row2_C);
%
% f = VBR_sols(1).VBR.ISV.f ;
% % AndradePsP, two frequencies
% for j = 1:nlines
%         %LineW = LineW_vec(j);
%         clr = colorscale(j,:) ;
%         %f_M = (VBR_sols(j).VBR.Gu_vec(1))/(VBR_sols(j).VBR.eta_total(1)) ;
%         %I_fM = find(f>=f_M,1);
%         %M = VBR_sols(j).VBR.eBurgers.M(:,1) ;
%
%         %plot(log10(f),M./1e9,'k-','LineWidth', LineW, 'Color', clr); hold on;
%         %plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;
%
%         % PLOT DATA:
%         plot(data.TanJax.exptCond.logf,data.TanJax.Results.G,'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;
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
%
%
%
% return
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
