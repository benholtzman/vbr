function PLOT_bestProfiles_f(Box,VarInfo,i_Var1,j_Var2,tsnap,obsVs,obsdepth,ylimits,figName)
figname = 'FIG3_BestProfs_CHANGENAME'

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
%% Patch the frequencies



%% ================================
% PLOTTING  setting
%  ================================
% figure size
%figsize_inches = [0.1 0.6 8 6];

top = 0.3 ; 
%bot = 0.1 ; 
Lmarg = 0.1 ;
dx = 0.017 ; 
W = 0.16 ; 
H = 0.60 ;

% left bottom width height
topA = [Lmarg top W H ] ;
topB = [Lmarg+dx+W top W H ] ;
topC = [Lmarg+2*dx+2*W top W H ] ;
topD = [Lmarg+3*dx+3*W top W H ] ;
topE = [Lmarg+4*dx+4*W top W H ] ;
%botA = [Lmarg bot W H ] ;
%botB = [Lmarg+dx+W bot W H ] ;
%botC = [Lmarg+2*dx+2*W bot W H ] ;

figure; %('color',[1 1 1])
%set(gcf, 'Units', 'inches');
%set(gcf, 'Position', figsize_inches);


%set(gcf, 'Units', 'inches');
%set(gcf, 'Position', [0.1 0.1 11 4]);
 
% patchline(V_VBR, Z_km,'linestyle','-','edgecolor',heat,'linewidth',1,'edgealpha',1); hold on;
%% =================
% Plot T and stuff
% =============

% ===========================================
% TEMPERATURE
axes('position', topA,'box','on'); hold on; 
%subplot(1,5,1)   
   plot(Box(i_Var1,j_Var2).Movie.Frames(tsnap).T-273,Z_km,'k')
   hold on
   iLAB=Box(i_Var1,j_Var2).Movie.info.i_LAB_vec(tsnap);
   plot([0 2000],[Z_km(iLAB) Z_km(iLAB)],'--k')
   hold off
   xlabel('T (C)'); 
   ylabel('z (km)'); 
   set(gca,'ydir','rev')  
   
   
% =========================================== 
% Modulus      
axes('position', topB,'box','on'); hold on;   
   Ma1 = Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.AndradePsP.Ma(1,:)*1e-9;
   Ma2 = Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.AndradePsP.Ma(end,:)*1e-9;   
   h=fill([Ma1 fliplr(Ma2)],[Z_km' fliplr(Z_km')],VBR_color,'FaceAlpha', 1);
   set(h,'EdgeColor',VBR_color_edge);
   
   hold on
   plot(Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.Gu*1e-9,Z_km,'k','linestyle','--') 
   plot([20 115],[Z_km(iLAB) Z_km(iLAB)],'--k')
%    plot(Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.AndradePsP.Ma*1e-9,Z_km)
   hold off
   xlabel('m (GPa)'); 
   %ylabel('z (km)');
   xlim([45 95])
   set(gca,'ydir','rev')  
   set(gca, 'YTick', []);
   
% ===========================================  
% ATTENUATION    
axes('position', topC,'box','on'); hold on;
%    Q = log10(Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.AndradePsP.Qa(:,:));
   Q1 = log10(Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.AndradePsP.Qa(1,:));
   Q2 = log10(Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.AndradePsP.Qa(end,:));
   h=fill([Q1 fliplr(Q2)],[Z_km' fliplr(Z_km')],VBR_color,'FaceAlpha', 1);
   set(h,'EdgeColor',VBR_color_edge);
   hold on
   plot([0 4],[Z_km(iLAB) Z_km(iLAB)],'--k')
   hold off
%    plot(Q,Z_km)
   xlabel('log_1_0 Q'); 
   %ylabel('z (km)'); 
   set(gca,'ydir','rev')   
   set(gca, 'YTick', []);    
   xlim([0 4])
   title(title_name);
   
   
% ===========================================
% Velocity
axes('position', topD,'box','on'); hold on; 
   V1 = Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.AndradePsP.Va(1,:)*1e-3;
   V2 = Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.AndradePsP.Va(end,:)*1e-3;   
   h=fill([V1 fliplr(V2)],[Z_km' fliplr(Z_km')],VBR_color,'FaceAlpha', 1);
   set(h,'EdgeColor',VBR_color_edge);

%    plot(Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.AndradePsP.Va/1e3,Z_km)
   hold on 
   plot(obsVs,obsdepth,'k','linewidth',2)
   plot([3 6],[Z_km(iLAB) Z_km(iLAB)],'--k')
   hold off
   xlabel('V_s (km s^{-1})'); 
   %ylabel('z (km)'); 
   set(gca,'ydir','rev') 
   xlim([4 5])
   set(gca, 'YTick', []);
   
% ===========================================
% VISCOSITY 
axes('position', topE,'box','on'); hold on; 
   eta_log = log10(Box(i_Var1,j_Var2).Movie.Frames(tsnap).VBR.eta_total(:,:));
   plot(eta_log,Z_km,'k')
   hold on
   plot([14 24],[Z_km(iLAB) Z_km(iLAB)],'--k')
   hold off
   xlabel('log_1_0 eta (Pa s)'); 
   %ylabel('z (km)'); 
   set(gca,'ydir','rev')    
   xlim([17 22])
   set(gca, 'YTick', []);
   
   
%for ii = 1:5
%    subplot(1,5,ii)
%    ylim(ylimits)
%end


% ===========================================
print(gcf,'-depsc',figName)

end

