clf; 
addpath('./datasets') ;
load('VBR_sols')

%COLORS

% left bottom width height
W = 0.33 ;
H = 0.36 ;
plot_row1_A = [0.1 0.60 W H ] ;
plot_row1_B = [0.55 0.60 W H ] ;
plot_row2_C = [0.1 0.12 W H ] ;
plot_row2_D = [0.55 0.12 W H ] ;


LBLFNT = 18 ;
TITLFNT = 16 ; 
LineW = 2 ; 
dotsize = 12; 
dotsize_D = 20 ;


%% PLOT =======================================================
% MAKE THE DATA LISTS
Make_DATA ;


%%  Q vs FREQUENCY (TanJax - eBurgers/AndradePsP fits) ==================================================
axes('Position', plot_row1_A); 

ExptSet = 1 ;

f = VBR_sols(1,1).VBR.ISV.f ;  
sz_dataset = size(Data.TanJax);
nT = sz_dataset(2);

%colorscale = [ 1 0.5 0.2 ; 1 0 0 ] ;
%cool to warm: 
colorscale(:,1) = linspace(0.5,1,nT) ;
colorscale(:,2) = linspace(0,0,nT) ;
colorscale(:,3) = linspace(1,0,nT) ;

for jT = 1:nT        
       %LineW = LineW_vec(j);
       clr = colorscale(jT,:) ;
       Qa = VBR_sols(ExptSet,jT).VBR.AndradePsP.Qa(1,:) ; % these axes flipped from 0.93 to 0.94
       %Qcomp = VBR_sols(ExptSet,jT).VBR.AndradePsP.Q_comp(:) ;  
       Qeb = VBR_sols(ExptSet,jT).VBR.eBurgers.Q(1,:) ;
  
        
       plot(log10(f),log10(1./Qa),'k--','LineWidth', LineW-1, 'Color', clr); hold on;  
       plot(log10(f),log10(1./Qeb),'k-','LineWidth', LineW, 'Color', clr); hold on; 
        %plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;
        
        % PLOT DATA: 
        f_data = Data.TanJax(jT).exptCond.f ;
        Qinv_data = Data.TanJax(jT).Results.Qinv ; 
        plot(log10(f_data), log10(Qinv_data),'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;

end

axis tight
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])
%title(['Andrade, P=' num2str(P) ' GPa, d=' num2str(d) ' \mum'],'fontname','Times New Roman','fontsize',LBLFNT);
title(['Tan/Jackson data, And.+eBurg. fit'],'fontname','Times New Roman','fontsize',TITLFNT);
xlabel('log frequency', 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('log Q^{-1}, (J_1/J_2)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);


%%  GribbCoop/SundCoop data, Andrade fit ==================================================
axes('Position', plot_row2_C); 

ExptSet = 2 ;

f = VBR_sols(1,1).VBR.ISV.f ;  
sz_dataset = size(Data.GribCoop);
nT = sz_dataset(2);

%colorscale = [ 1 0.5 0.2 ; 1 0 0 ] ;
%cool to warm: 
clear colorscale
colorscale(:,1) = linspace(0.5,1,nT) ;
colorscale(:,2) = linspace(0,0,nT) ;
colorscale(:,3) = linspace(1,0,nT) ;

for jT = 1:nT        
       %LineW = LineW_vec(j);
       clr = colorscale(jT,:) ;
       Qa = VBR_sols(ExptSet,jT).VBR.AndradePsP.Qa(1,:) ;
       %Qcomp = VBR_sols(ExptSet,jT).VBR.AndradePsP.Q_comp(:) ;  
       Qeb = VBR_sols(ExptSet,jT).VBR.eBurgers.Q(1,:) ;
  
        
       plot(log10(f),log10(1./Qa),'k--','LineWidth', LineW-1, 'Color', clr); hold on;  
       plot(log10(f),log10(1./Qeb),'k-','LineWidth', LineW, 'Color', clr); hold on; 
        %plot(log10(f(I_fM)),log10(Qs(I_fM)),'r.', 'MarkerSize',dotsize); hold on;
        
        % PLOT DATA: 
        f_data = Data.GribCoop(jT).exptCond.f ;
        Qinv_data = Data.GribCoop(jT).Results.Qinv ; 
        plot(log10(f_data), log10(Qinv_data),'g.', 'MarkerSize',dotsize_D, 'Color', clr); hold on;

end

axis tight
%xlim([1.8e1 3e2])
%ylim([1e-6 5e-4])
%title(['Andrade, P=' num2str(P) ' GPa, d=' num2str(d) ' \mum'],'fontname','Times New Roman','fontsize',LBLFNT);
title(['Gribb/Cooper data, And.+eBurg PsP fit'],'fontname','Times New Roman','fontsize',TITLFNT);
xlabel('log frequency', 'fontname','Times New Roman','fontsize', LBLFNT)
ylabel('log Q^{-1}, (J_1/J_2)', 'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'fontname','Times New Roman','fontsize', LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on','ticklength',[0.03 0.03],'linewidth',1);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%





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




