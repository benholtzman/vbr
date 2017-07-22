%wtf%% Plot it!    
   zplt = settings.Zinfo.z_km; 
   if exist('yrange','var')==0; yrange = [min(zplt) max(zplt)]; end
   figure('color',[1 1 1])
  
   subplot(1,4,1)
      for it = 1:numel(Info.t)
          Xvar=Vars.T(:,it);
          plot(Xvar,zplt,'marker','none')
          hold all 
      end    
      plot(Vars.Tsol(:,1),zplt,'--k'); 
      set(gca,'ydir','rev'); ylim(yrange)
      xlabel('T [C]'); ylabel('depth [km]')
      
   subplot(1,4,2)
      for it = 1:numel(Info.t)
          Xvar=Vars.phi(:,it);
          plot(Xvar,zplt,'marker','none')
          hold all
      end
      set(gca,'ydir','rev'); ylim(yrange)
      xlabel('\phi'); ylabel('depth [km]')
         
   subplot(1,4,3)
      for it = 1:numel(Info.t)
          Xvar=Vars.eta_nonad(:,it)./Vars.eta_nonad(end,1);
          semilogx(Xvar,zplt,'marker','none')
          hold all 
      end    
      semilogx([10 10],[min(zplt) max(zplt)],'--k')
      xlim([0.1 1000])
      set(gca,'ydir','rev'); ylim(yrange)
      xlabel('non-adiabatic \eta / \eta_o'); ylabel('depth [km]')
      
   subplot(1,4,4)
      for it = 1:numel(Info.t)
          Xvar=Vars.eta(:,it)./Vars.eta(end,1);
          semilogx(Xvar,zplt,'marker','none')
          hold all 
      end    
      semilogx([10 10],[min(zplt) max(zplt)],'--k')
      xlim([0.1 1000])
      set(gca,'ydir','rev'); ylim(yrange)
      xlabel('\eta / \eta_o'); ylabel('depth [km]')
      
