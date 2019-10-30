
  
% thermal_states,tau_points,freq_points

% % store the output  
%   Time.TimeFast=TimeFast;
%   Time.TimeSlow=TimeSlow;
%   Err.ErrJ1=ErrJ1;
%   Err.ErrJ2=ErrJ2;
%   Err.ErrM=ErrM;
%   Err.ErrVs=ErrVs;
%   Err.ErrQinv=ErrQinv;  
  
  tau_points=Ranges.tau_points; 
  thermal_states=Ranges.thermal_states;
  freq_points=Ranges.freq_points;
  
  figure('color',[1 1 1])
  
  abscissa='thermal_states';
  
  if strcmp(abscissa,'tau_points')
      subplot(1,2,1)
      plot(tau_points,Time.TimeFast(1,:,1),'r','displayname','FastBurger','marker','.');
      hold on
      plot(tau_points,Time.TimeSlow(1,:,1),'k','displayname','SlowBurger','marker','.');
      xlabel('n_{\tau}')
      ylabel('elapsed time [s]')
      legend('location','northwest')
      title(['n_{thermal states}=' num2str(thermal_states) ', n_{freq}=' num2str(freq_points)])
      
      subplot(1,2,2)
      semilogy(tau_points,Err.ErrVs(1,:,1)+1e-50,'k','displayname','V_s','marker','.');
      hold all
      semilogy(tau_points,Err.ErrM(1,:,1)+1e-50,'displayname','M','marker','.');
      semilogy(tau_points,Err.ErrJ1(1,:,1)+1e-50,'displayname','J_1','marker','.');
      semilogy(tau_points,Err.ErrJ2(1,:,1)+1e-50,'displayname','J_2','marker','.');
      semilogy(tau_points,Err.ErrQinv(1,:,1)+1e-50,'displayname','Q^{-1}','marker','.');
      xlabel('n_{\tau}')
      ylabel('err')
      legend('location','northwest')
  elseif strcmp(abscissa,'freq_points')
      subplot(1,2,1)
      plot(freq_points,squeeze(Time.TimeFast(1,1,:)),'r','displayname','FastBurger','marker','.');
      hold on
      plot(freq_points,squeeze(Time.TimeSlow(1,1,:)),'k','displayname','SlowBurger','marker','.');
      xlabel('n_{freq}')
      ylabel('elapsed time [s]')
      legend('location','northwest')
      title(['n_{thermal states}=' num2str(thermal_states) ', n_{tau}=' num2str(tau_points)])
      
      subplot(1,2,2)
      semilogy(freq_points,squeeze(Err.ErrVs(1,1,:))+1e-50,'k','displayname','V_s','marker','.');
      hold all
      semilogy(freq_points,squeeze(Err.ErrM(1,1,:))+1e-50,'displayname','M','marker','.');
      semilogy(freq_points,squeeze(Err.ErrJ1(1,1,:))+1e-50,'displayname','J_1','marker','.');
      semilogy(freq_points,squeeze(Err.ErrJ2(1,1,:))+1e-50,'displayname','J_2','marker','.');
      semilogy(freq_points,squeeze(Err.ErrQinv(1,1,:))+1e-50,'displayname','Q^{-1}','marker','.');
      xlabel('n_{freq}')
      ylabel('err')
      legend('location','northwest')
  elseif strcmp(abscissa,'thermal_states')
      
      ntau=numel(tau_points); 
      nfreq= numel(freq_points);       
      nclr = nfreq; 
      
%       titlename=['n_{freq}=' num2str(freq_points) ', (min,max) freq: (' ...
%           num2str(10^Ranges.freqmin) ',' num2str(10^Ranges.freqmax) ')'];
             
      titlename=['n_{tau}=' num2str(tau_points) ', (min,max) freq: (' ...
          num2str(10^Ranges.freqmin) ',' num2str(10^Ranges.freqmax) ')'];
%       for ip = 1:ntau
      for ifreq = 1:nfreq
          ip = 1;           
          iclr = ifreq; 
          
          clr = [(iclr-1)/(nclr-1) 0 1-(iclr-1)/(nclr-1)]; 
          
%           Dn=num2str(tau_points(ip)); 
          Dn=num2str(freq_points(ifreq)); 
          
          subplot(2,3,1)
          if iclr > 1; hold on; end
          semilogy(thermal_states,squeeze(Time.TimeFast(:,ip,ifreq)),'color',clr,'displayname','FastBurger','marker','o');
          hold on
          semilogy(thermal_states,squeeze(Time.TimeSlow(:,ip,ifreq)),'color',clr,'displayname','SlowBurger','marker','.');
          ylabel('elapsed time [s]')          
          title(titlename)
          
          subplot(2,3,2)
          if iclr > 1; hold on; end
          semilogy(thermal_states,squeeze(Err.ErrJ1(:,ip,ifreq))+1e-50,'color',clr,'displayname',Dn,'marker','.');
          ylabel('err(J_1)')
          
          subplot(2,3,3)
          if iclr > 1; hold on; end
          semilogy(thermal_states,squeeze(Err.ErrJ2(:,ip,ifreq))+1e-50,'color',clr,'displayname',Dn,'marker','.');
          ylabel('err(J_2)')
          
          subplot(2,3,4)
          if iclr > 1; hold on; end
          semilogy(thermal_states,squeeze(Err.ErrM(:,ip,ifreq))+1e-50,'color',clr,'displayname',Dn,'marker','.');
          ylabel('err(M)')
          
          subplot(2,3,5)
          if iclr > 1; hold on; end
          semilogy(thermal_states,squeeze(Err.ErrQinv(:,ip,ifreq))+1e-50,'color',clr,'displayname',Dn,'marker','.');
          ylabel('err(Q^{-1})')
          
          subplot(2,3,6)
          if iclr > 1; hold on; end
          semilogy(thermal_states,squeeze(Err.ErrVs(:,ip,ifreq))+1e-50,'color',clr,'displayname',Dn,'marker','.');                    
          ylabel('err(V_s)')
          legend('location','northwest')
          
      end
      
      for ip = 1:6
          subplot(2,3,ip)
          xlabel('N_{thermal states}')
          
      end
  end
  
  set(gcf, 'Units', 'inches');
  x_width=12 ;y_width=6;
  set(gcf, 'Position', [0 0 x_width y_width]); %
  
  
