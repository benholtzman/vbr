clf;
addpath('./datasets') ;
load('VBR_sols')
nlines = length(VBR_sols(1).T_params) ;

%COLORS

% left bottom width height
W = 0.33 ;
H = 0.36 ;
plot_row1_A = [0.1 0.60 W H ] ;
plot_row1_B = [0.55 0.60 W H ] ;
plot_row2_C = [0.1 0.12 W H ] ;
plot_row2_D = [0.55 0.12 W H ] ;


LBLFNT = 18 ;
LineW = 2 ;
LineW_vec = linspace(1,3,nlines);
dotsize = 12;
dotsize_D = 20 ;

VBR = VBR_sols(1).VBR ;
d = VBR.SSV.dg_um(1) ;
P = VBR.ISV.P_GPa(1) ;

%% PLOT =======================================================
% MAKE THE DATA LISTS
Make_DATA ;

for j = 1:nlines
    data(j).T = VBR_sols(j).T_params ;
    data(j).logf = Data.SundCoop(j).exptCond.logf ;
    data(j).logQ = log10(1./Data.SundCoop(j).Results.Qinv) ;
    data(j).G_GPa = Data.SundCoop(j).Results.G ;
end

%colorscale = [ 1 0.5 0.2 ; 1 0 0 ] ;
%cool to warm:
colorscale(:,1) = linspace(0.5,1,nlines) ;
colorscale(:,2) = linspace(0,0,nlines) ;
colorscale(:,3) = linspace(1,0,nlines) ;

%%  Q vs FREQUENCY (Andrade) ==================================================
axes('Position', plot_row1_A);

f = VBR_sols(1).VBR.ISV.f ;
% Andrade, two frequencies
for j = 1:nlines
        %LineW = LineW_vec(j);
        clr = colorscale(j,:) ;
       % f_M = (VBR_sols(j).VBR.Gu_vec(1))/(VBR_sols(j).VBR.eta_total(1)) ;
       % I_fM = find(f>=f_M,1);
       Q = VBR_sols(j).VBR.AndradePsP.Qa(:) ;
       Qcomp = VBR_sols(j).VBR.AndradePsP.Q_comp(:) ;
       % MAT(:,j) = Qs ;

       plot(log10(f),log10(Q),'k--','LineWidth', LineW-1, 'Color', clr); hold on;
       plot(log10(f),log10(Qcomp),'k-','LineWidth', LineW, 'Color', clr); hold on;
        %plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;

        % PLOT DATA:
        plot(data(j).logf,data(j).logQ,'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;

end

axis tight
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])
title(['Andrade, P=' num2str(P) ' GPa, d=' num2str(d) ' \mum'],'fontname','Times New Roman','fontsize',LBLFNT);
xlabel('log frequency', 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('log Q, (J_1/J_2)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);



%%  G vs FREQUENCY (Andrade) ==================================================
axes('Position', plot_row2_C);

% Andrade, two frequencies
for j = 1:nlines
        %LineW = LineW_vec(j);
        clr = colorscale(j,:) ;
        %f_M = (VBR_sols(j).VBR.Gu_vec(1))/(VBR_sols(j).VBR.eta_total(1)) ;
        %I_fM = find(f>=f_M,1);
        Ma = VBR_sols(j).VBR.AndradePsP.Ma(:) ;
        Mcomp = VBR_sols(j).VBR.AndradePsP.M_comp(:) ;
        %MAT(:,j) = Qs ;

        plot(log10(f),Ma./1e9,'k--','LineWidth', LineW-1, 'Color', clr); hold on;
        plot(log10(f),Mcomp./1e9,'k-','LineWidth', LineW, 'Color', clr); hold on;
        %plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;

        % PLOT DATA:
        plot(data(j).logf,data(j).G_GPa,'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;

end

axis tight
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])

xlabel('log frequency', 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('Modulus, M (GPa)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);
%set(gca,'XTickLabel', [])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  Q vs FREQUENCY (BURGERS) ==================================================
axes('Position', plot_row1_B);

f = VBR_sols(1).VBR.ISV.f ;
% Andrade, two frequencies
for j = 1:n_temp
        %LineW = LineW_vec(j);
        clr = colorscale(j,:) ;
        %f_M = (VBR_sols(j).VBR.Gu_vec(1))/(VBR_sols(j).VBR.eta_total(1)) ;
        %I_fM = find(f>=f_M,1);
        Qs = VBR_sols(j).VBR.eBurgers.Q(:,1) ;
        %MAT(:,j) = Qs ;

        plot(log10(f),log10(Qs),'k-','LineWidth', LineW, 'Color', clr); hold on;
        %plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;

        % PLOT DATA:
        plot(data(j).logf,data(j).logQ,'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;

end

axis tight
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])
title(['Extended Burgers'],'fontname','Times New Roman','fontsize',LBLFNT);
xlabel('log frequency', 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('log Q, (J_1/J_2)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);




%%  G vs FREQUENCY (BURGERS) ==================================================
axes('Position', plot_row2_D);

f = VBR_sols(1).VBR.ISV.f ;
% Andrade, two frequencies
for j = 1:n_temp
        %LineW = LineW_vec(j);
        clr = colorscale(j,:) ;
        %f_M = (VBR_sols(j).VBR.Gu_vec(1))/(VBR_sols(j).VBR.eta_total(1)) ;
        %I_fM = find(f>=f_M,1);
        M = VBR_sols(j).VBR.eBurgers.M(:,1) ;
        %MAT(:,j) = Qs ;

        plot(log10(f),M./1e9,'k-','LineWidth', LineW, 'Color', clr); hold on;
        %plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;

        % PLOT DATA:
        plot(data(j).logf,data(j).G,'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;

end

axis tight
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])

xlabel('log frequency', 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('relaxed Modulus, G (GPa)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);
%set(gca,'XTickLabel', [])
