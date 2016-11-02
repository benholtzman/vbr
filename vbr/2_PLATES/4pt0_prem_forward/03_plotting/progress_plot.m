function fprog = progress_plot(Vark,z,fprog)

    z= z/1e3; 

    if isempty(fprog.fighand)
      fprog.fighand=figure('color',[1 1 1]);
    else
      figure(fprog.fighand);
    end
    
    subplot(2,2,2);
      plot(Vark.Cbulk,z,'k','displayname','Cbulk');
      hold on; 
      plot(Vark.Strg,z,'--r','displayname','Strg');
      plot(Vark.Cparg,z,'displayname','Cparg');
      hold off; 
      xlabel('H_2O [wt %]')
      legend('location','southeast')
    subplot(2,2,4);
      plot(Vark.Cbulk,z);
      xlim([0 0.1])
      xlabel('Bulk H_2O [wt %]')
    subplot(2,2,1);
      plot(Vark.T,z,'displayname','T'); 
      hold on; 
      plot(Vark.Tsol,z,'--r','displayname','T_{sol}'); 
      hold off; 
      xlim([800 max([max(Vark.T) max(Vark.Tsol)])])
      xlabel('[^oC]')
      legend('location','southwest')
    subplot(2,2,3);
      plot(Vark.phi,z); 
      maxphi = max([max(Vark.phi) 0.02]);
      xlim([0 maxphi])
      xlabel('\phi')

    for ip = 1:4
     subplot(2,2,ip)
     set(gca,'ydir','rev')
     ylim([0 max(z)])
    end

   
   set(gcf, 'PaperUnits', 'inches');
   set(gcf, 'PaperPosition', [0 0 fprog.x_width fprog.y_width]); %
   
   pause(fprog.prog_pause)
   saveas(gcf,'Progress.eps','epsc')

end
