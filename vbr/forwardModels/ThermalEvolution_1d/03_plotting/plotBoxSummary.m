function  plotBoxSummary(Box,settings,Work,varargin)


  % default options
  ValidOpts=struct('plot_every_k',[],'depth_range',[],'plot_every_dt',[],...
                   'isotherm_Vals',[],'savenames',[],'iBoxes',[]);
  Options=struct('plot_every_k',10,...
    'depth_range',[0,max(Box(1).run_info.Z_km)],...
    'plot_every_dt',0,'isotherm_Vals',[1200],'savenames',{'none'},'iBoxes',[0]);
  Options=validateStructOpts('plotSummary',varargin,Options,ValidOpts);

  % plot single boxes
  for iBox_indx=1:numel(Options.iBoxes)
    iBox=Options.iBoxes(iBox_indx);
    if iBox > 0
      [Vars,Info,settings] = pullFromBox(Box,iBox);
      plotSummary(Vars,Info,'plot_every_k',Options.plot_every_k,...
        'plot_every_dt',Options.plot_every_dt,'depth_range',Options.depth_range);
    end
  end

  % plot summary of all boxes
  sumFig=figure();
  subplot(1,3,1)
  hold on
  if isfield(settings.Box,'nvar2')
    for iv2 = 1: settings.Box.nvar2
      cf=(iv2 - 1) / (settings.Box.nvar2 -1);
      plot(Box(1,iv2).Frames(end).Tsol,Box(1,iv2).run_info.Z_km,'color',[0,cf,0],'linestyle','--')
    end
  end
  for iv1 = 1:settings.Box.nvar1
    cf=(iv1 - 1) / (settings.Box.nvar1 -1);
    plot(Box(iv1,1).Frames(end).T,Box(iv1,1).run_info.Z_km,'color',[cf,0,0])
  end
  box on
  set(gca,'ydir','reverse')
  xlabel('temperature [C]')
  ylabel('depth [km]')
  title('Final profiles')


  subplot(1,3,2)
  hold on
  for iBox=1:Work.nBox
    cf=(iBox - 1) / (Work.nBox  -1);
    plot(Box(iBox).run_info.tMyrs,Box(iBox).run_info.zSOL/1000,'color',[cf,0,0],'marker','.')
  end
  box on
  ylim([0,settings.Zinfo.asthenosphere_max_depth+10])
  xlabel('model time [Myr]')
  ylabel('z_{SOL} [km]')

  subplot(1,3,3)
  hold on
  for iBox=1:Work.nBox
    cf=(iBox - 1) / (Work.nBox  -1);
    plot(sqrt(Box(iBox).run_info.tMyrs),Box(iBox).run_info.zSOL/1000,'color',[cf,0,0],'marker','.')
  end
  box on
  ylim([0,settings.Zinfo.asthenosphere_max_depth+10])
  xlabel('sqrt(model time [Myr])')
  ylabel('z_{SOL} [km]')

  % if ~strcmp(Options.savename,'none')
  %   disp(['saving fig to ',Options.savename])
  % end

end

function zIso=extractIsothermDepth(Box0,isoVal)

  tsteps=Box0.run_info.tMyrs;
  zIso=zeros(size(tsteps));
  z=Box0.run_info.Z_km;
  for it=1:numel(tsteps)
    T=Box0.Frames(it).T;
    zIso(it)=min(z(T>=isoVal));
  end
end
