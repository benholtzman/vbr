

function PLOT_bestProfiles_f(Box,iTemp,iZLAB,tsnap,obsVs,obsdepth,ylimits)
Z_km = Box(1,1).Movie.info.Z_km ; 
%% =================
% Plot T and stuff
% =============
figure('color',[1 1 1])
subplot(1,5,1)   
   plot(Box(iTemp,iZLAB).Movie.Frames(tsnap).T-273,Z_km)
   xlabel('T (C)'); ylabel('z (km)'); title(['best fit, ' ...
       num2str(Box(iTemp,iZLAB).info.T_pot_C) '^oC, ' ...
       num2str(Box(iTemp,iZLAB).info.ZLAB/1e3) ' km'])
   set(gca,'ydir','rev')
% subplot(1,4,2)   
%    plot(Box(iTemp,iZLAB).Movie.Frames(tsnap).rho,Z_km)
%    xlabel('\rho [kg m^{-3}]'); ylabel('z (km)'); 
%    set(gca,'ydir','rev')   
subplot(1,5,2)  
   plot(Box(iTemp,iZLAB).Movie.Frames(tsnap).VBR.Gu*1e-9,Z_km,'linestyle','--')
   hold on
   plot(Box(iTemp,iZLAB).Movie.Frames(tsnap).VBR.AndradePsP.Ma*1e-9,Z_km)
   hold off
   xlabel('m (GPa)'); ylabel('z (km)'); 
   set(gca,'ydir','rev')        
subplot(1,5,3)  
   Q = log10(Box(iTemp,iZLAB).Movie.Frames(tsnap).VBR.AndradePsP.Qa(:,:));
   plot(Q,Z_km)
   xlabel('log Q'); ylabel('z (km)'); 
   set(gca,'ydir','rev')       
   xlim([0 4])
subplot(1,5,4)  
   plot(Box(iTemp,iZLAB).Movie.Frames(tsnap).VBR.AndradePsP.Va/1e3,Z_km)
   hold on 
   plot(obsVs,obsdepth,'k','linewidth',2)
   hold off
   xlabel('V_s (km s^{-1})'); ylabel('z (km)'); 
   set(gca,'ydir','rev') 
   xlim([4 5])
subplot(1,5,5)     
   xlabel('eta (Pa s)'); ylabel('z (km)'); 
   set(gca,'ydir','rev')    

for ii = 1:5
    subplot(1,5,ii)
    ylim(ylimits)
end

end

