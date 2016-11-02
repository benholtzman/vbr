function PLOT_GIA_beta(Box,VarInfo,i_Var1,j_Var2,obsVs,obsdepth,frame_with_VBR,ylimits)

figname = 'FIG_VBR_beta'

% load Box data
  Box_dir = '../../../0_BOXES/';
  Box_name ='Box_120Ma_plate_shallow_PTdep_BnR_y150520_GIA_VBR';
  Box_name = [Box_dir Box_name]
  load(Box_name)

Z_km = Box(1,1).Movie.info.Z_km ; 

if VarInfo.Var2_n > 1
title_name = ['Best fit: ' VarInfo.Var1_name '=' num2str(VarInfo.Var1_range(i_Var1)) ...
              ' ' VarInfo.Var1_units ', ' VarInfo.Var2_name '=' num2str(VarInfo.Var2_range(j_Var2)) ...
              ' ' VarInfo.Var2_units];
else
title_name = ['Best fit: ' VarInfo.Var1_name '=' num2str(VarInfo.Var1_range(i_Var1)) ...
              ' ' VarInfo.Var1_units];
end

VBR_color = [0 0.6 0.1];
VBR_color_edge = VBR_color;%[0 0 0];

%% number of frequencies:
sz_vbr = size(Box(1,8).Movie.Frames(10).VBR.AndradePsP.Qa); 
nf = sz_vbr(1) ; 

i_f_sw_min = 12 ; 
i_f_sw_max = 18 ;

i_f_GIA_min = 1 ; 
i_f_GIA_max = 10 ;

%% PLOTTING ====================================
figure; %('color',[1 1 1])
%set(gcf, 'Units', 'inches');
%set(gcf, 'Position', [0.1 0.5 8 4]);


H = 0.6 ; 
W = 0.22 ; 
dx = 0.03 ;
lmarg = 0.07 ;
bot = 0.3 ;
% figure size
  figsize_inches = [0.1 0.6 8 6];
% left bottom width height
plot_A = [lmarg bot W/2 H ] ; 
plot_B = [lmarg+dx+0.5*W bot W H ] ;
plot_C = [lmarg+2*dx+1.5*W bot W H];
plot_D = [lmarg+3*dx+2.5*W bot W H ] ;  
plot_E = [lmarg+4*dx+3.5*W bot W H ] ; 
  
LBLFNT = 14 ;

colorsR = linspace(0,1,nf) ;
colorsG = linspace(0,0,nf) ;
colorsB = linspace(1,0,nf) ;

bestColor = [0 0.5 0] ;

colors =  [colorsR ; colorsG ; colorsB] ;
LineW_Std = 1.0; 
LineW_Emphasize = 3;
 %linspace(4,2,nj) ; 

% ======================================
% TEMPERATURE 
axes('position', plot_A); hold on; 


T_C = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).T-273 ;
linecolor = [0.8 0 0] ;
plot(T_C,Z_km, 'LineWidth', LineW_Std, 'Color', linecolor); hold on;   

    set(gca,'YDir','reverse')
    xlabel('T (C)','fontname','Times New Roman','fontsize',LBLFNT); 
    ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
    ylim(ylimits)
    set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
    axis tight
    
% ======================================
% Vs (for surface wave frequencies) 
axes('position', plot_B); hold on; 


for i_freq = i_f_sw_min:i_f_sw_max 
	Vs_rel = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Va(i_freq,:) ;
	linecolor = (colors(:,i_freq))' ; 
	plot(Vs_rel./1e3,Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end

Vsu = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Vsu ; 
plot(Vs_rel./1e3,Z_km, 'k--'); hold on; 

titletext = ['log_{10} f = ' num2str(i_f_sw_min) ' - ' num2str(i_f_sw_max) ] ;
set(gca, 'YTick', []);
    set(gca,'YDir','reverse')
    title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
    xlabel('V_s [km/s]','fontname','Times New Roman','fontsize',LBLFNT); 
    %ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
    ylim(ylimits)
    set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
    xlim([3.9 5.1])
    
% ======================================
% M (for GIA frequencies) 
axes('position', plot_C); hold on; 


for i_freq = i_f_GIA_min:i_f_GIA_max 
	M_rel = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Ma(i_freq,:) ;
	linecolor = colors(:,i_freq)' ; 
	plot(M_rel./1e9,Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end

titletext = ['log_{10} f = ' num2str(i_f_GIA_min) ' - ' num2str(i_f_GIA_max) ]  ;

set(gca, 'YTick', []);
    set(gca,'YDir','reverse')
    title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
    xlabel('M [GPa]','fontname','Times New Roman','fontsize',LBLFNT); 
    %ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
    ylim(ylimits)
    set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
    axis tight

% ======================================
% M (for GIA frequencies) 
axes('position', plot_D); hold on; 

% anharmonic shear modulus (T,P dependent) ; 
GuTP = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Gu ; 

for i_freq = i_f_GIA_min:i_f_GIA_max 
	M_rel = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Ma(i_freq,:)' ;
	Mrel_ManhTP = M_rel./GuTP ; 
	linecolor = colors(:,i_freq)' ; 
	plot(Mrel_ManhTP,Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end

set(gca, 'YTick', []);
    set(gca,'YDir','reverse')
    xlabel('M/M(T,P,f=inf.)','fontname','Times New Roman','fontsize',LBLFNT); 
    %ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
    ylim(ylimits)
    set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
    axis tight
% ======================================

print(gcf,'-depsc',figname)

end

