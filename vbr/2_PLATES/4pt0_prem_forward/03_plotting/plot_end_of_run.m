%% Plot it!    
   zplt = Info.z_km; 
   if exist('yrange','var')==0; yrange = [0 max(zplt)];end
   
   
if strcmp(settings.Flags.problem,'One_Phase_T')==0   
    figure('color',[1 1 1])
   subplot(1,4,2)
      for it = 1:3:numel(Info.t)
          Xvar=Vars.phi(:,it);
          plot(Xvar,zplt,'marker','none')
          hold all
      end
      plot(Xvar,zplt,'--k','marker','none','linewidth',2)      
      set(gca,'ydir','rev'); ylim(yrange)
      xlabel('\phi'); ylabel('depth [km]')
      
   subplot(1,4,1)
      
       for it = 1:3:numel(Info.t)
          Xvar=Vars.T(:,it);
          plot(Xvar,zplt,'marker','none')
          hold all
%            Xvar=Vars.Tsol(:,it);
%           plot(Xvar,zplt,'marker','none','linestyle','--')
      end
      plot(Vars.Tsol(:,1),zplt,'--k','displayname','T_{sol}','linewidth',2)
      plot(Vars.Tsol(:,end),zplt,'--k','displayname','T_{sol}','linewidth',2)
      plot(Vars.T(:,end),zplt,'b','displayname','final','linewidth',2)
      set(gca,'ydir','rev'); ylim(yrange)
      xlabel('[^oC]'); ylabel('depth [km]')

      subplot(1,4,3)       
      for it = 1:2:numel(Info.t)
          Xvar=Vars.eta(:,it);
          Xvar(Xvar>1e25)=1e25; 
          semilogx(Xvar,zplt,'marker','none','linewidth',2)
          hold all 
      end    
      semilogx(Xvar,zplt,'--k','marker','none','linewidth',2)
      semilogx([min(Xvar) max(Xvar)],[Info.zLABeta(end) Info.zLABeta(end)]/1e3,'--k')
      semilogx([min(Xvar) max(Xvar)],[Info.zLAB(end) Info.zLAB(end)]/1e3,'--r')
      set(gca,'ydir','rev'); ylim(yrange); xlim([1e19 1e25])
      set(gca,'xtick',[10.^(19:24)])
      xlabel('\eta [Pa s]'); ylabel('depth [km]')
     
     
    subplot(1,4,4) 
      for it = 1:2:numel(Info.t)
          Xvar=Vars.Cs(:,it); 
          plot(Xvar,zplt,'marker','none')
          hold all 
      end      
      plot(Xvar,zplt,'--k','marker','none','linewidth',2)
      set(gca,'ydir','rev'); ylim(yrange)
      xlabel('C_s [wt %]'); ylabel('depth [km]')
      xlim([0 max([10*1e-4 max(max(Vars.Cs))])])
      
      for ip = 1:2
              subplot(1,4,ip)
              hold on
              plot(get(gca,'xlim'),[Info.zLAB(end) Info.zLAB(end)]*1e-3,'--k')
              plot(get(gca,'xlim'),[Info.zLAB(1) Info.zLAB(1)]*1e-3,'--r')
      end
      
       x_width=8 ;y_width=3;
else
   figure('color',[1 1 1])
      
    subplot(1,2,1)
    
    for it = 1:3:numel(Info.t)
        Xvar=Vars.T(:,it);
        plot(Xvar,zplt,'marker','none')
        hold all
    end
    plot(Vars.Tsol(:,1),zplt,'--k','displayname','T_{sol}','linewidth',2)
    plot(Vars.Tsol(:,end),zplt,'--k','displayname','T_{sol}','linewidth',2)
    plot(Vars.T(:,end),zplt,'b','displayname','final','linewidth',2)
    set(gca,'ydir','rev'); ylim(yrange)
    xlabel('[^oC]'); ylabel('depth [km]')
    
   
      
    subplot(1,2,2)
    for it = 1:2:numel(Info.t)
        Xvar=Vars.eta(:,it);
        Xvar(Xvar>1e25)=1e25;
        semilogx(Xvar,zplt,'marker','none','linewidth',2)
        hold all
    end
    semilogx(Xvar,zplt,'--k','marker','none','linewidth',2)
    semilogx([min(Xvar) max(Xvar)],[Info.zLABeta(end) Info.zLABeta(end)]/1e3,'--k')
    semilogx([min(Xvar) max(Xvar)],[Info.zLAB(end) Info.zLAB(end)]/1e3,'--r')
    set(gca,'ydir','rev'); ylim(yrange); xlim([1e19 1e25])
    set(gca,'xtick',[10.^(19:24)])
    
    xlabel('\eta [Pa s]'); ylabel('depth [km]')
    x_width=4 ;y_width=3;

end

   set(gcf, 'PaperUnits', 'inches');

 set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
 
 
 
     if exist([savedir '/individual_runs'],'dir')~=7
         mkdir([savedir '/individual_runs']);
     end
     if exist([savedir '/individual_runs/' savebase],'dir')~=7
         mkdir([savedir '/individual_runs/' savebase]);
     end
     indifile=[savebase '_' num2str((ivar2+(ivar1-1)*nvar2))];
     savename = [savedir '/individual_runs/' savebase '/' indifile];
     saveas(gcf,[savename '.eps'],'epsc')
 
 
%  
%  figure('color',[1 1 1])
%  plot(Info.tMyrs,Info.zLABeta/1e3)
%  xlabel('t [Myrs]')
%  ylabel('z_{LAB} [km]')
%  set(gca,'ydir','rev')

