function F2=PLOT_wholeBox_ch(Box,Obs,Fit_Params,ylimits)

i_Var1=Fit_Params.best_Var1;
j_Var2=Fit_Params.best_Var2;

sz = size(Box);
nj = sz(1) ; % variable DEPTH
nk = sz(2) ; % variable dT 


frame_with_VBR=Fit_Params.Frame_Selection;

SP(1).xlab='T [^oC]';
SP(1).xlims=[0 1500];
SP(2).xlab='M [GPa]';
SP(2).xlims=[40 90];
SP(3).xlab='log Q';
SP(3).xlims=[0 4];
SP(4).xlab='V_s [km/s]';
SP(4).xlims=[3 8];

%% PLOTTING ====================================
F2=figure('color',[1 1 1]);
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.1 0.5 8 4]);

LBLFNT = 14 ;

colorsR = linspace(0,1,nj) ;
colorsG = linspace(0,0,nj) ;
colorsB = linspace(1,0,nj) ;
bestColor = [0 0.6 0] ;

colors =  [colorsR ; colorsG ; colorsB] ;
LineW_Standard = 0.5; 
LineW_Emphasize = 3;

for k = 1:nk
    for j = 1:nj
        
        VBRan=Box(j,k).Frames(frame_with_VBR(j,k)).VBR.out.anelastic;
        
        T_C = Box(j,k).Frames(frame_with_VBR(j,k)).T-273 ;		       
        M = VBRan.(Fit_Params.VBR_anelastic_method).Ma./1e9;		
        Q = log10(VBRan.(Fit_Params.VBR_anelastic_method).Qa);
        Vs = (VBRan.(Fit_Params.VBR_anelastic_method).Vave)./1e3;              
       
        M=mean(M,2);
        Q=mean(Q,2);
        
        Z_km = Box(j,k).run_info.Z_km ; 
                
        if j == i_Var1 && k == j_Var2
            Best(1).var=T_C;
            Best(2).var=M;
            Best(3).var=Q;
            Best(4).var=Vs; 
            Best(5).var=Z_km;
        end
        
        LineW =LineW_Standard;
        clr = colors(:,j) ;
        
        subplot(1,4,1)
        hold on 
        plot(T_C,Z_km, 'LineWidth', LineW, 'Color', clr)
        ylabel('Z [km]','fontname','Times New Roman','fontsize',LBLFNT)
        
        subplot(1,4,2)
        hold on
        plot(M,Z_km, 'LineWidth', LineW, 'Color', clr)
        
        subplot(1,4,3)
        hold on 
        plot(Q,Z_km, 'LineWidth', LineW, 'Color', clr)
        
        subplot(1,4,4)
        hold on
        plot(Vs,Z_km, 'LineWidth', LineW, 'Color', clr); 
    end
end


for ip = 1:4
    subplot(1,4,ip)
    hold on 
    LineW=LineW_Emphasize;
    clr=bestColor;
    plot(Best(ip).var,Best(5).var, 'LineWidth', LineW, 'Color', clr)
    
    set(gca,'box','on','xminortick','on','yminortick','on','Ydir','rev',...
            'fontname','Times New Roman','fontsize', LBLFNT);
    ylim(ylimits)
    xlabel(SP(ip).xlab,'fontname','Times New Roman','fontsize',LBLFNT); 
    xlim(SP(ip).xlims)
    axis tight
end   

subplot(1,4,4)
hold on
plot(Obs.Vs/1e3,Obs.depth,'--k','LineWidth',1.5)



end
    
