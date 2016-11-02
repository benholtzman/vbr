yrange = [25 70];
%% Save settings
   savedir='../../../../0_BOXES/individual_runs';
   savename = [savedir '/TestH2O_1.mat'];
   
   videoname='TestH2O.avi'; 
   videoname = [pwd '/' videoname];
   vidObj = VideoWriter(videoname);
   open(vidObj);

%% Load it!
   disp('Loading:')
   disp(savename)
   
   if exist(savename,'file')==0
       disp('file does not exist')
       
   else
       disp(savename)
       load(savename)
       
       %%     Plot it!
       t_vec=Info.tMyrs;
       disp('Plotting')
       zplt = Info.z_km;
%        yrange = [min(zplt) max(zplt)];
       yrange = [30 150];
       [val,yrangeid(1)]=min(abs(zplt-yrange(1)));       
       %%     Watch it!
       if exist('MOVVV','var')==1; clear MOVVV; end
       x_width=12 ;y_width=5;
       f2=figure('color',[1 1 1],'units','inches','position',[.1 .1 x_width y_width]);
       for it = 1:numel(t_vec)
           subplot(1,4,1)
            Xvar=Vars.T(:,it);
            plot(Xvar,zplt,'r','marker','none','linewidth',2)
            hold on
            plot(Vars.Tsol(:,it),zplt,'--k','linewidth',2)
            plot(Vars.Tsol(:,1),zplt,'--k','linewidth',1)
            hold off
            set(gca,'ydir','rev'); ylim(yrange)
            xlim([round(Vars.T(yrangeid(1))/100)*100 max(max(Vars.T))+100])
            xlabel('[^oC]'); ylabel('dedpth [km]')
           
           subplot(1,4,2)      
            Xvar=Vars.phi(:,it);
            plot(Xvar,zplt,'marker','none','linewidth',2)
            set(gca,'ydir','rev'); ylim(yrange); 
            xlim([0 max(max(round(100*Vars.phi)/100))]);
            xlabel('\phi'); ylabel('depth [km]')
            
          
      
          subplot(1,4,3)
            Xvar=Vars.eta(:,it);
            Xvar(Xvar>1e25)=1e25;
            semilogx(Xvar,zplt,'marker','none','linewidth',2)
            hold on       
            semilogx([min(Xvar) max(Xvar)],[Info.zLABeta(it) Info.zLABeta(it)]/1e3,'--k')
            semilogx([min(Xvar) max(Xvar)],[Info.zLAB(it) Info.zLAB(it)]/1e3,'--r')
            hold off
            set(gca,'ydir','rev'); ylim(yrange); xlim([1e19 1e25])
            set(gca,'xtick',[10.^(19:24)])
            xlabel('\eta [Pa s]'); ylabel('depth [km]')     
      
         subplot(1,4,4)
            Xvar=Vars.Cs(:,it); 
            plot(Xvar,zplt,'marker','none','linewidth',2)            
            set(gca,'ydir','rev'); ylim(yrange)
            xlabel('C_s [wt %]'); ylabel('depth [km]')
            xlim([0 max(max(Vars.Cs))])
           
            
         set(gcf, 'PaperUnits', 'inches');        
         set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
         MOVVV(it)=getframe(f2);
         writeVideo(vidObj,MOVVV(it));
       end        
       
       pause(0.2)
       MOVVV(it)=getframe(f2);
       writeVideo(vidObj,MOVVV(it));
       
       close(vidObj);
   end