%% Plot it!    
clear
dir = '../../../../0_BOXES/individual_runs/';
fbase = 'Vup_con_nomelting_';
yrange = [0 160];
figure('color',[1 1 1])
   
for ifile = 1:5
   fname = [dir fbase num2str(ifile) '.mat']
   load(fname)   
   zplt = Info.z_km; 
   dnm = [num2str(settings.Vbg) ' cm/yr']; 
   
   subplot(1,2,1)
      if ifile == 1; 
      plot(Vars.T(:,1),zplt,'--b','displayname','T(0)','linewidth',2)
      hold on
      plot(Vars.Tsol(:,1),zplt,'--k','displayname','T_{sol}','linewidth',2)      
      hold off
      end
          hold all
      
      
      plot(Vars.T(:,end),zplt,'displayname',dnm,'linewidth',2)
      set(gca,'ydir','rev'); ylim(yrange)
      xlabel('[^oC]'); ylabel('depth [km]')
      
   subplot(1,2,2) 
      if ifile > 1; hold all; end
      Xvar=Vars.Gdot(:,end);%.*(Vars.phi(:,it)>settings.phimin*10);
      plot(Xvar,zplt,'displayname',dnm,'linewidth',2)            
      set(gca,'ydir','rev'); ylim(yrange)
      xlabel('\Gamma [s^{-1}]'); ylabel('depth [km]')          
      if ifile > 1; hold off; end
      
end