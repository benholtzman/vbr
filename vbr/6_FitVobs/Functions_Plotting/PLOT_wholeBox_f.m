function PLOT_wholeBox_f(Box,VarInfo,i_Var1,j_Var2,obsVs,obsdepth,frame_with_VBR,ylimits,figName)

sz = size(Box);
nj = sz(1) ; % variable DEPTH
nk = sz(2) ; % variable dT 

sz_f = length(Box(1,1).Movie.info.freq_vec) ;
ni = sz_f ; %  number of frequencies
i = 1 ; % the frequency i wanna plot


%% PLOTTING ====================================
figure('color',[1 1 1])
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.1 0.5 8 4]);
Z_km = Box(1,1).Movie.info.Z_km ; 
LBLFNT = 14 ;

colorsR = linspace(0,1,nj) ;
colorsG = linspace(0,0,nj) ;
colorsB = linspace(1,0,nj) ;
bestColor = [0 0.6 0] ;

colors =  [colorsR ; colorsG ; colorsB] ;
LineW_Standard = 0.5; 
LineW_Emphasize = 3;
 %linspace(4,2,nj) ; 

% ======================
subplot(1,4,1)
for k = 1:nk
    for j = 1:nj
%         [val,frame_with_VBR]=min(abs(Box(j,k).Movie.info.timesteps_myrs - target_age)); 
        T_C = Box(j,k).Movie.Frames(frame_with_VBR).T-273 ;
		%T_C = Box(j,k).Movie.Frames(frame_with_VBR).VBR.ISV.T_K-273 ;
        LineW =LineW_Standard;
		linecolor = colors(:,j) ;
        plot(T_C,Z_km, 'LineWidth', LineW, 'Color', linecolor); hold on;          
    end
end
T_C = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).T-273 ;
plot(T_C,Z_km, 'LineWidth', LineW_Emphasize, 'Color', bestColor); hold on;

    set(gca,'YDir','reverse')
    xlabel('T (C)','fontname','Times New Roman','fontsize',LBLFNT); 
    ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
    ylim(ylimits)
    set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
    axis tight


% ======================
subplot(1,4,2)


for k = 1:nk
    for j = 1:nj
        M = (Box(j,k).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Ma(i,:))./1e9;
		%M = (Box(j,k).Movie.Frames(frame_with_VBR).VBR.Gu)./1e9;
                LineW =LineW_Standard;
		linecolor = colors(:,j) ;
        plot(M,Z_km, 'LineWidth', LineW, 'Color', linecolor); hold on;  
    end
end
M = (Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Ma(i,:))./1e9 ;
plot(M,Z_km, 'LineWidth', LineW_Emphasize, 'Color', bestColor); hold on;

set(gca,'YDir','reverse')
xlabel('M (GPa)','fontname','Times New Roman','fontsize',LBLFNT); 
%ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
ylim(ylimits)
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
axis tight


% ======================
subplot(1,4,3)

for k = 1:nk
    for j = 1:nj
        Q = log10(Box(j,k).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Qa(i,:));
        LineW =LineW_Standard;
		linecolor = colors(:,j) ;
        plot(Q,Z_km, 'LineWidth', LineW, 'Color', linecolor); hold on;  
    end
end
logQ = log10(Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Qa(i,:)) ;
plot(logQ,Z_km, 'LineWidth', LineW_Emphasize, 'Color', bestColor); hold on;

set(gca,'YDir','reverse')
ylim(ylimits)
xlabel('log Q','fontname','Times New Roman','fontsize',LBLFNT); 
%ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
xlim([0 4])


    
% ======================
subplot(1,4,4)
minVs = 8;
maxVs = 5; 
for k = 1:nk
    for j = 1:nj
        Vs = (Box(j,k).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Va(i,:))./1e3;
        minVs = min([min(Vs) minVs]);
        maxVs = max([max(Vs) maxVs]);
        LineW =LineW_Standard;
		linecolor = colors(:,j) ;       
        plot(Vs,Z_km, 'LineWidth', LineW, 'Color', linecolor); hold on;  
    end
end
Vs = (Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Va(i,:))./1e3 ;
plot(Vs,Z_km, 'LineWidth', LineW_Emphasize, 'Color', bestColor); hold on;

plot(obsVs,obsdepth, 'k-', 'LineWidth', 3 ); hold on;
ylim(ylimits)
set(gca,'YDir','reverse')
xlabel('Vs (km/s)','fontname','Times New Roman','fontsize',LBLFNT); 
%ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
xlim([minVs 4.8])    

print(gcf,'-depsc',figName)
end
    
