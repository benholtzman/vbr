function PLOT_GIA_veryVBR_beta(Box,VarInfo,i_Var1,j_Var2,obsVs,obsdepth,frame_with_VBR,ylimits)
clf; % delete this when making multiple plots.. 
figname = 'FIG_GIA_veryVBR_beta'

% load Box data
%  Box_dir = '../../../0_BOXES/';
%  Box_name ='Box_120Ma_plate_shallow_PTdep_BnR_y150520_GIA_VBR';
%  Box_name = [Box_dir Box_name]
%  load(Box_name)

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
sz_vbr = size(Box(i_Var1,j_Var2).Movie.Frames(1).VBR.ISV.f); 
nf = sz_vbr(2) ; 

% need a way to bring these in automatically.. store in VBR.info ?  doesnt exist yet... 
i_f_GIA_min = 1 ; 
i_f_GIA_max = 10 ;

i_f_int_min = 11 ; 
i_f_int_max = 20 ;

i_f_sw_min = 22 ; 
i_f_sw_max = 28 ;


%% PLOTTING ====================================
% figure size
figsize_inches = [0.1 0.1 6 8];
%figure; %('color',[1 1 1])
%set(gcf, 'Units', 'inches');
%set(gcf, 'Position', [0.1 0.5 8 4]);
%set(gcf, 'Position', figsize_inches);

H = 0.17 ; 
W = 0.18 ; 
dx = 0.02 ;
lmarg = 0.1 ;
row1_bot = 0.8 ;
row2_bot = 0.54 ;
row3_bot = 0.3 ;
row4_bot = 0.06 ;

% left bottom width height
row1_A = [lmarg row1_bot 0.85*W H ] ; 
row1_B = [lmarg+dx+W row1_bot 0.85*W H ] ;
row1_C = [lmarg+2*dx+2*W row1_bot 0.85*W H];
row1_D = [lmarg+3*dx+3*W row1_bot 0.55*W H ] ;  

row2_A = [lmarg row2_bot W H ] ; 
row2_B = [lmarg+1*dx+1*W row2_bot W H ] ; 
row2_C = [lmarg+2*dx+2*W row2_bot W H ] ; 

row3_A = [lmarg row3_bot W H ] ; 
row3_B = [lmarg+1*dx+1*W row3_bot W H ] ; 
row3_C = [lmarg+2*dx+2*W row3_bot W H ] ; 

row4_A = [lmarg row4_bot W H ] ; 
row4_B = [lmarg+1*dx+1*W row4_bot W H ] ; 
row4_C = [lmarg+2*dx+2*W row4_bot W H ] ;  
  
LBLFNT = 9 ;

colorsR = linspace(0,1,nf) ;
colorsG = linspace(0,0,nf) ;
colorsB = linspace(1,0,nf) ;

bestColor = [0 0.5 0] ;

colors =  [colorsR ; colorsG ; colorsB] ;
LineW_Std = 1.0; 
LineW_Emphasize = 3;
 %linspace(4,2,nj) ; 

% ======================================================================
% TOP ROW...  frequency independent things... 
% ======================================================================
% ==========================
% TEMPERATURE 
axes('position', row1_A); hold on; 

T_C = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).T-273 ;
linecolor = [0.8 0 0] ;
plot(T_C,Z_km, 'LineWidth', LineW_Std, 'Color', linecolor); hold on;   

    set(gca,'YDir','reverse')
    xlabel('T (C)','fontname','Times New Roman','fontsize',LBLFNT); 
    ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
    ylim(ylimits)
    set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
    axis tight

% ==========================
axes('position', row1_B); hold on; 

Gu = (Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Gu)./1e9 ;
plot(Gu,Z_km, 'k-','LineWidth', LineW_Std); hold on;   

    set(gca,'YDir','reverse')
    xlabel('Mu [GPa]','fontname','Times New Roman','fontsize',LBLFNT); 
    %ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
    ylim(ylimits)
    set(gca, 'YTick', []);
    set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
    axis tight

% ==========================
axes('position', row1_C); hold on; 
linecolor = [0.8 0 0] ;

eta = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.eta_total ;
eta_cap = 1e24 ;
eta(find(eta>=eta_cap)) = eta_cap ; 
log10_eta = log10(eta) ; 
plot(log10_eta,Z_km, 'LineWidth', LineW_Std, 'Color', linecolor); hold on;   

    set(gca,'YDir','reverse')
    xlabel('log_{10} \eta [Pa.s]','fontname','Times New Roman','fontsize',LBLFNT); 
    %ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
    ylim(ylimits)
    set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
    axis tight
    set(gca, 'YTick', []);

% ======================================================================
% ROW 2.. J1 
% ======================================================================  
% ======================================
% J1 @ GIA low frequencies... 
axes('position', row2_A); hold on; 

for i_freq = i_f_GIA_min:i_f_GIA_max
	J1 = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.J1(i_freq,:) ;
	M1 = (1./J1)./1e9 ; 
	linecolor = (colors(:,i_freq))' ; 
	plot(M1,Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end

%Vsu = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Vsu ; 
%plot(Vs_rel./1e3,Z_km, 'k--'); hold on; 

titletext = 'GIA (low) freqs' ;
%titletext = ['log_{10} f = ' num2str(i_f_sw_min) ' - ' num2str(i_f_sw_max) ] ;
%set(gca, 'YTick', []);
set(gca,'YDir','reverse')
ylim(ylimits)
% ====
title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
xlabel('M1 [GPa]','fontname','Times New Roman','fontsize',LBLFNT); 
ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT);
% ====
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
axis tight
xlim([0 82])
  
% ======================================
% J1 @ Seasonal (intermediate) frequencies... 
axes('position', row2_B); hold on;  

for i_freq = i_f_int_min:i_f_int_max 
	J1 = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.J1(i_freq,:) ;
	M1 = (1./J1)./1e9 ; 
	linecolor = (colors(:,i_freq))' ; 
	plot(M1,Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end

%Vsu = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Vsu ; 
%plot(Vs_rel./1e3,Z_km, 'k--'); hold on; 

titletext = 'Seasonal (int) freqs' ;
%titletext = ['log_{10} f = ' num2str(i_f_sw_min) ' - ' num2str(i_f_sw_max) ] ;
set(gca, 'YTick', []);
set(gca,'YDir','reverse')
ylim(ylimits)
% ====
title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
xlabel('M1 [GPa]','fontname','Times New Roman','fontsize',LBLFNT); 
%ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT);
% ====
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
axis tight
xlim([0 82])


% ======================================
% J1 @ Seismic (high) frequencies... 
axes('position', row2_C); hold on; 

for i_freq = i_f_sw_min:i_f_sw_max 
	J1 = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.J1(i_freq,:) ;
	M1 = (1./J1)./1e9 ; 
	linecolor = (colors(:,i_freq))' ; 
	plot(M1,Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end

%Vsu = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Vsu ; 
%plot(Vs_rel./1e3,Z_km, 'k--'); hold on; 

titletext = 'Seismic (high) freqs' ;
%titletext = ['log_{10} f = ' num2str(i_f_sw_min) ' - ' num2str(i_f_sw_max) ] ;
set(gca, 'YTick', []);
set(gca,'YDir','reverse')
ylim(ylimits)
% ====
title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
xlabel('M1 [GPa]','fontname','Times New Roman','fontsize',LBLFNT); 
%ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT);
% ====
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
axis tight
xlim([0 82])
% ======================================

% ======================================================================
% ROW 3.. J2 
% ====================================================================== 
low_f = 1e-12 ; 
M2_max  = 1e6/low_f ;
% ======================================
% J2 @ GIA low frequencies... 
axes('position', row3_A); hold on; 

for i_freq = i_f_GIA_min:i_f_GIA_max
	J2 = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.J2(i_freq,:) ;
	M2 = (1./J2) ;
	M2(find(M2>=M2_max)) = M2_max ;  
	linecolor = (colors(:,i_freq))' ; 
	plot(log10(M2),Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end

%Vsu = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Vsu ; 
%plot(Vs_rel./1e3,Z_km, 'k--'); hold on; 

titletext = 'GIA (low) freqs' ;
%titletext = ['log_{10} f = ' num2str(i_f_sw_min) ' - ' num2str(i_f_sw_max) ] ;
%set(gca, 'YTick', []);
set(gca,'YDir','reverse')
ylim(ylimits)
% ====
%title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
xlabel('log_{10}M2 [Pa]','fontname','Times New Roman','fontsize',LBLFNT); 
ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT);
% ====
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
axis tight
xlim([5 log10(M2_max)])
  
% ======================================
% J2 @ Seasonal (intermediate) frequencies... 
axes('position', row3_B); hold on;  

for i_freq = i_f_int_min:i_f_int_max 
	J2 = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.J2(i_freq,:) ;
	M2 = (1./J2) ;
	M2(find(M2>=M2_max)) = M2_max ; 
	linecolor = (colors(:,i_freq))' ; 
	plot(log10(M2),Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end

%Vsu = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Vsu ; 
%plot(Vs_rel./1e3,Z_km, 'k--'); hold on; 

%titletext = 'Seasonal (int) freqs' ;
%titletext = ['log_{10} f = ' num2str(i_f_sw_min) ' - ' num2str(i_f_sw_max) ] ;
set(gca, 'YTick', []);
set(gca,'YDir','reverse')
ylim(ylimits)
% ====
%title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
xlabel('log_{10}M2 [Pa]','fontname','Times New Roman','fontsize',LBLFNT); 
%ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT);
% ====
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
axis tight
xlim([5 log10(M2_max)])


% ======================================
% J1 @ Seismic (high) frequencies... 
axes('position', row3_C); hold on; 

for i_freq = i_f_sw_min:i_f_sw_max 
	J2 = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.J2(i_freq,:) ;
	M2 = (1./J2) ;
	M2(find(M2>=M2_max)) = M2_max ;  
	linecolor = (colors(:,i_freq))' ; 
	plot(log10(M2),Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end

%Vsu = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Vsu ; 
%plot(Vs_rel./1e3,Z_km, 'k--'); hold on; 

%titletext = 'Seismic (high) freqs' ;
%titletext = ['log_{10} f = ' num2str(i_f_sw_min) ' - ' num2str(i_f_sw_max) ] ;
set(gca, 'YTick', []);
set(gca,'YDir','reverse')
ylim(ylimits)
% ====
%title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
xlabel('log_{10}M2 [Pa]','fontname','Times New Roman','fontsize',LBLFNT); 
%ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT);
% ====
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
axis tight
xlim([5 log10(M2_max)])

% ======================================

% ======================================================================
% ROW 4.. M or other measureable...  
% ======================================================================  
% ======================================
% M @ GIA low frequencies... 
axes('position', row4_A); hold on; 

for i_freq = i_f_GIA_min:i_f_GIA_max
	Ma = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Ma(i_freq,:) ;
	linecolor = (colors(:,i_freq))' ; 
	plot(Ma./1e9,Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end
plot(Gu,Z_km, 'k-','LineWidth', LineW_Std); hold on;  
%Vsu = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Vsu ; 
%plot(Vs_rel./1e3,Z_km, 'k--'); hold on; 

titletext = 'GIA (low) freqs' ;
%titletext = ['log_{10} f = ' num2str(i_f_sw_min) ' - ' num2str(i_f_sw_max) ] ;
%set(gca, 'YTick', []);
set(gca,'YDir','reverse')
ylim(ylimits)
% ====
%title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
xlabel('M_{rlx} [GPa]','fontname','Times New Roman','fontsize',LBLFNT); 
ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT);
% ====
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
axis tight
xlim([0 82])
  
% ======================================
% M @ Seasonal (intermediate) frequencies... 
axes('position', row4_B); hold on;  

for i_freq = i_f_int_min:i_f_int_max 
	Ma = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Ma(i_freq,:) ;
	linecolor = (colors(:,i_freq))' ; 
	plot(Ma./1e9,Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end
plot(Gu,Z_km, 'k-','LineWidth', LineW_Std); hold on;  
%Vsu = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Vsu ; 
%plot(Vs_rel./1e3,Z_km, 'k--'); hold on; 

%titletext = 'Seasonal (int) freqs' ;
%titletext = ['log_{10} f = ' num2str(i_f_sw_min) ' - ' num2str(i_f_sw_max) ] ;
set(gca, 'YTick', []);
set(gca,'YDir','reverse')
ylim(ylimits)
% ====
%title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
xlabel('M_{rlx} [GPa]','fontname','Times New Roman','fontsize',LBLFNT); 
%ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT);
% ====
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
axis tight
xlim([0 82])


% ======================================
% M @ Seismic (high) frequencies... 
axes('position', row4_C); hold on; 

for i_freq = i_f_sw_min:i_f_sw_max 
	Ma = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Ma(i_freq,:) ;
	linecolor = (colors(:,i_freq))' ; 
	plot(Ma./1e9,Z_km, 'Color', linecolor, 'LineWidth',LineW_Std); hold on; 
end
plot(Gu,Z_km, 'k-','LineWidth', LineW_Std); hold on;  

%Vsu = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.Vsu ; 
%plot(Vs_rel./1e3,Z_km, 'k--'); hold on; 

%titletext = 'Seismic (high) freqs' ;
%titletext = ['log_{10} f = ' num2str(i_f_sw_min) ' - ' num2str(i_f_sw_max) ] ;
set(gca, 'YTick', []);
set(gca,'YDir','reverse')
ylim(ylimits)
% ====
%title(titletext,'fontname','Times New Roman','fontsize',LBLFNT); 
xlabel('M_{rlx} [GPa]','fontname','Times New Roman','fontsize',LBLFNT); 
%ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT);
% ====
set(gca,'box','on','xminortick','on','yminortick','on', 'fontname','Times New Roman','fontsize', LBLFNT);
axis tight
xlim([0 82])
% ======================================
%	M_rlx = Box(i_Var1,j_Var2).Movie.Frames(frame_with_VBR).VBR.AndradePsP.Ma(i_freq,:) ;

print(gcf,'-depsc',figname)

end

