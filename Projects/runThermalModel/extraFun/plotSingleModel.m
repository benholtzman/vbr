function fig=plotSingleModel(Vars,Info)
  fig=figure();
  subplot(1,2,1)
  plot(Vars(:,end).Tsol,Info.z_km,'--k')
  hold on
  for it = 1:2: numel(Info.tMyrs)
    cf=(it - 1) / (numel(Info.tMyrs)-1);
    plot(Vars.T(:,it),Info.z_km,'color',[0,0,cf])
  end
  box on
  set(gca,'ydir','reverse')
  xlabel('geotherms and solidii [C]')
  ylabel('depth [km]')

  subplot(1,2,2)
  plot(Info.tMyrs,Info.zSOL/1000,'k','linewidth',2)
  hold on
  plot(Info.tMyrs,Info.zLAB/1000,'--k','linewidth',2)
  for it = 1: numel(Info.tMyrs)
    cf=(it - 1) / (numel(Info.tMyrs)-1);
    plot(Info.tMyrs(it),Info.zSOL(it)/1000,'color',[0,0,cf])
    plot(Info.tMyrs(it),Info.zLAB(it)/1000,'color',[0,0,cf])
  end
  box on
  xlabel('model time [Myr]')
  ylabel('zSOL, zLAB [km]')
end
