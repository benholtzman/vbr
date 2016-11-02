function F3=PLOT_bestProfiles_ch(Box,VarInfo,Fit_Params,Obs)

obsVs=Obs.Vs/1e3; 
obsdepth=Obs.depth; 

i_Var1=Fit_Params.best_Var1;
j_Var2=Fit_Params.best_Var2;

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

%% ================================
% PLOTTING  setting
%  ================================

top = 0.3 ; 
%bot = 0.1 ; 
Lmarg = 0.1 ;
dx = 0.017 ; 
W = 0.16 ; 
H = 0.60 ;

% location=[left bottom width height]
topA = [Lmarg top W H ] ;
topB = [Lmarg+dx+W top W H ] ;
topC = [Lmarg+2*dx+2*W top W H ] ;
topD = [Lmarg+3*dx+3*W top W H ] ;
topE = [Lmarg+4*dx+4*W top W H ] ;

F3=figure('color',[1 1 1]);

%% =================
% Plot T and stuff
% =============

Frm = Fit_Params.Frame_Selection(i_Var1,j_Var2);
Z_km = Box(i_Var1,j_Var2).run_info.Z_km;
% ===========================================
% TEMPERATURE
axes('position', topA,'box','on'); hold on;    
   plot(Box(i_Var1,j_Var2).Frames(Frm).T-273,Z_km,'k')
   xlabel('T [C]'); 
   ylabel('z [km]'); 
   set(gca,'ydir','rev')  
   
   
% =========================================== 
% Modulus      
axes('position', topB,'box','on'); hold on; 
   VBR=Box(i_Var1,j_Var2).Frames(Frm).VBR.out.anelastic.(Fit_Params.VBR_anelastic_method);
   Ma1 = VBR.Ma(:,1)*1e-9;
   Ma2 = VBR.Ma(:,end)*1e-9;   
   h=fill([Ma1' fliplr(Ma2')],[Z_km' fliplr(Z_km')],VBR_color,'FaceAlpha', 1);
   set(h,'EdgeColor',VBR_color_edge);
   

   xlabel('M [GPa]'); 
   xlim([45 95])
   set(gca,'ydir','rev')  
   set(gca, 'YTick', []);
   
% ===========================================  
% ATTENUATION    
axes('position', topC,'box','on'); hold on;
   Q1 = log10(VBR.Qa(:,1));
   Q2 = log10(VBR.Qa(:,end));
   h=fill([Q1' fliplr(Q2')],[Z_km' fliplr(Z_km')],VBR_color,'FaceAlpha', 1);
   set(h,'EdgeColor',VBR_color_edge);
   xlabel('log_1_0 Q'); 
   set(gca,'ydir','rev')   
   set(gca, 'YTick', []);    
   xlim([0 4])
   title(title_name);
   
   
% ===========================================
% Velocity
axes('position', topD,'box','on'); hold on; 
   V1 = VBR.Va(:,1)*1e-3;
   V2 = VBR.Va(:,end)*1e-3;   
   h=fill([V1' fliplr(V2')],[Z_km' fliplr(Z_km')],VBR_color,'FaceAlpha', 1);
   set(h,'EdgeColor',VBR_color_edge);
   hold on 
   plot(obsVs,obsdepth,'k','linewidth',2)

   hold off
   xlabel('V_s (km s^{-1})');     
   set(gca,'ydir','rev') 
   xlim([4 5])
   set(gca, 'YTick', []);
   
% ===========================================
% VISCOSITY 
axes('position', topE,'box','on'); hold on; 
   visc=Box(i_Var1,j_Var2).Frames(Frm).VBR.out.viscous.(Fit_Params.VBR_visc_method);   
   semilogx(visc.eta_total,Z_km,'k')
   xlabel('\eta [Pa s]'); 
   
   set(gca,'ydir','rev','xscale','log')    
   xlim(10.^[19 22])
   set(gca, 'YTick', []);
  
end

