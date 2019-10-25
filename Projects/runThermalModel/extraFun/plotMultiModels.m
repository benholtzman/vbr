function fig = plotMultiModels(Box,iBox,settings,Work)

  fig=figure();

  subplot(2,2,1)
  hold on
  plot(Box(iBox).Frames(end).Tsol,Box(iBox).run_info.Z_km,'--k')
  for it = 1:2: numel(Box(iBox).run_info.tMyrs)
    cf=(it - 1) / (numel(Box(iBox).run_info.tMyrs)-1);
    plot(Box(iBox).Frames(it).T,Box(iBox).run_info.Z_km,'color',[0,0,cf])
  end
  box on
  set(gca,'ydir','reverse')
  xlabel('geotherms and solidii [C]')
  ylabel('depth [km]')

  subplot(2,2,2)
  zIso=extractIsotherm(Box(iBox),1200);

  plot(Box(iBox).run_info.tMyrs,Box(iBox).run_info.zSOL/1000,'k','linewidth',2)
  hold on
  plot(Box(iBox).run_info.tMyrs,Box(iBox).run_info.zLAB/1000,'--k','linewidth',2)
  plot(Box(iBox).run_info.tMyrs,zIso,'--r','linewidth',2)

  for it = 1: numel(Box(iBox).run_info.tMyrs)
    cf=(it - 1) / (numel(Box(iBox).run_info.tMyrs)-1);
    plot(Box(iBox).run_info.tMyrs(it),Box(iBox).run_info.zSOL(it)/1000,'color',[0,0,cf])
    plot(Box(iBox).run_info.tMyrs(it),Box(iBox).run_info.zLAB(it)/1000,'color',[0,0,cf])
    plot(Box(iBox).run_info.tMyrs(it),zIso(it),'color',[0,0,cf])
  end
  box on
  xlabel('model time [Myr]')
  ylabel('zSOL, zLAB [km]')


  subplot(2,2,3)
  hold on
  for iv2 = 1: settings.Box.nvar2
    cf=(iv2 - 1) / (settings.Box.nvar2 -1);
    plot(Box(1,iv2).Frames(end).Tsol,Box(1,iv2).run_info.Z_km,'color',[0,cf,0],'linestyle','--')
  end
  for iv1 = 1:settings.Box.nvar1
    cf=(iv1 - 1) / (settings.Box.nvar1 -1);
    plot(Box(iv1,1).Frames(end).T,Box(iv1,1).run_info.Z_km,'color',[cf,0,0])
  end
  box on
  set(gca,'ydir','reverse')
  xlabel('final geotherms and solidii [C]')
  ylabel('depth [km]')

  subplot(2,2,4)
  hold on
  for iBox=1:Work.nBox
    cf=(iBox - 1) / (Work.nBox  -1);
    plot(Box(iBox).run_info.tMyrs,Box(iBox).run_info.zSOL/1000,'color',[cf,0,0],'marker','.')
  end
  box on
  ylim([0,settings.Zinfo.asthenosphere_max_depth+10])
  xlabel('model time [Myr]')
  ylabel('solidus-geotherm intersection [km]')

end
