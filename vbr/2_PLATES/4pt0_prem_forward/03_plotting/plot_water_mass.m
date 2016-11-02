%% Plot it!    
   zplt = Info.z_km; 
   oplt = 'no';
     
   if exist('yrange','var')==0; yrange = [0 160];end
   if strcmp(oplt,'yes')==0
     f000=figure('color',[1 1 1]);
   else
     figure(f000)
   end
  
   subplot(2,1,1)
      H2O = Vars.Cs.*(1-Vars.phi) + Vars.Cf.*(Vars.phi); 
      for it = 1:5:numel(Info.t)
          Xvar=H2O(:,it); 
          plot(Xvar,zplt,'marker','none')
          hold all 
      end      
      plot(Xvar,zplt,'--k','marker','none','linewidth',2)
      plot(H2O(:,1),zplt,'--k','marker','none','linewidth',2)
      set(gca,'ydir','rev'); ylim(yrange)
      xlabel('Bulk H2O [wt %]'); ylabel('depth [km]')
      
    subplot(2,1,2)
      Deficit=zeros(size(Info.t)); 
      for it = 1:numel(Info.t); 
         Intgr8 = cumtrapz(zplt,H2O(:,it));
         Deficit(it)=Intgr8(end);
      end
      if strcmp(oplt,'yes')==1; hold all; end
      plot(Info.tMyrs,(Deficit-Deficit(1))./Deficit(1)); 
      if strcmp(oplt,'yes')==1; hold off; end
%       
%       H_Incoming = abs(Info.t(end)*Vars.Vbgz(end,end))*Vars.Cs(end,end) 
