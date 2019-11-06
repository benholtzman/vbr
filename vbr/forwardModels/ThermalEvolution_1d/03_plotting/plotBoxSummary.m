function  plotBoxSummary(Box,settings,varargin)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % fig = plotBoxSummary(Box,settings,varargin)
  %
  % builds summary plot for parameter sweep, plots single run results for
  % specified runs.
  %
  % Parameters
  % ----------
  %   Box       the box structure array
  %   settings  the box settings for the sweep
  %
  %
  %   optional keyword parameters:
  %   'plot_every_k',N    integer value to plot every N timesteps (default 10)
  %   'plot_every_dt',dt  scalar value for time in Myrs, will plot every dt.
  %   'depth_range',[zmin,zmax]   the depth range to plot, 2 element array
  %                               (default is whole domain)
  %   'savename',savename.png     file to save to. will split the filename and
  %                               save savename_summary.png and
  %                               savename_box_N.png for each box plotted. No
  %                               figures saved by default.
  %   'iBoxes',[1,3,4]            array of box numbers to plot, none by default
  %   'gridtime',5                the time at which to build the zSOL grid (5
  %                               Myrs is default)
  %
  % Output
  % ------
  %   none
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % default options
  ValidOpts=struct('plot_every_k',[],'depth_range',[],'plot_every_dt',[],...
                  'savename',[],'iBoxes',[],'gridtime',[]);
  Options=struct('plot_every_k',10,...
    'depth_range',[0,max(Box(1).run_info.Z_km)],...
    'plot_every_dt',0,'savename','none','iBoxes',[0],'gridtime',10);
  Options=validateStructOpts('plotSummary',varargin,Options,ValidOpts);

  % plot single boxes
  for iBox_indx=1:numel(Options.iBoxes)
    iBox=Options.iBoxes(iBox_indx);
    if iBox > 0
      [Vars,Info,settings] = pullFromBox(Box,iBox);
      if strcmp(Options.savename,'none')
        sname='none';
      else
        b=Options.savename;
        prfx=b(1:strfind(b,'.')-1);
        sffx=b(strfind(b,'.'):end);
        sname=[prfx,'_box_',num2str(iBox),sffx];
      end
      plotSummary(Vars,Info,'plot_every_k',Options.plot_every_k,...
        'plot_every_dt',Options.plot_every_dt,'depth_range',Options.depth_range,...
        'savename',sname);
    end
  end

  % plot summary of all boxes
  sumFig=figure();

  % temperature profiles
  subplot(2,2,1)
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

  % zSOl vs time
  subplot(2,2,2)
  Work.nBox=settings.Box.nvar1 * settings.Box.nvar2;
  hold on
  for iBox=1:Work.nBox
    cf=(iBox - 1) / (Work.nBox  -1);
    plot(Box(iBox).run_info.tMyrs,Box(iBox).run_info.zSOL/1000,'color',[cf,0,0],'marker','.')
  end
  box on
  ylim([0,settings.Zinfo.asthenosphere_max_depth+10])
  xlabel('model time [Myr]')
  ylabel('z_{SOL} [km]')

  % zSOL vs sqrt(time)
  subplot(2,2,4)
  hold on
  for iBox=1:Work.nBox
    cf=(iBox - 1) / (Work.nBox  -1);
    plot(sqrt(Box(iBox).run_info.tMyrs),Box(iBox).run_info.zSOL/1000,'color',[cf,0,0],'marker','.')
  end
  box on
  ylim([0,settings.Zinfo.asthenosphere_max_depth+10])
  xlabel('sqrt(model time [Myr])')
  ylabel('z_{SOL} [km]')


  % grid of zSOL at a given time
  zSOLgrid=buildGrid(Box,settings,Options.gridtime);
  subplot(2,2,3)
  if settings.Box.nvar2>1
    imagesc(zSOLgrid);
    set(gca,'xtick',1:1:settings.Box.nvar1,'xticklabel',settings.Box.var1range)
    set(gca,'ytick',1:1:settings.Box.nvar2,'yticklabel',settings.Box.var2range)
    colorbar
    title(['z_{SOL} at ',num2str(Options.gridtime),' Myrs'])
    xlabel([strrep(settings.Box.var1name,'_','\_'),settings.Box.var1units])
    ylabel([strrep(settings.Box.var2name,'_','\_'),settings.Box.var2units])
  else
    plot(settings.Box.var1range,zSOLgrid,'.k')
    xlabel([strrep(settings.Box.var1name,'_','\_'),settings.Box.var1units])
    ylabel(['z_{SOL} at ',num2str(Options.gridtime),' Myrs'])
  end

  % save it
  if ~strcmp(Options.savename,'none')
    b=Options.savename;
    prfx=b(1:strfind(b,'.')-1);
    sffx=b(strfind(b,'.'):end);
    sname=[prfx,'_summary',sffx];
    saveas(sumFig,sname)
  end

end

function zSOLgrid = buildGrid(Box,settings,ttarg)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % zSOLgrid = buildGrid(Box,settings,ttarg)
  %
  % builds the zSOL grid for a given time
  %
  % Parameters
  % ----------
  %   Box       the box structure array
  %   settings  the box settings for the sweep
  %   ttarg     the target time in Myrs to build the grid for
  %
  % Output
  % ------
  %   zSOLgrid  the matrix of zSOL at ttarg. If nvar2 == 1, will be 1d.
  %             otherwise, size(zSOLgrid)=(nvar2,nvar1)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % initialize the grid
  if settings.Box.nvar2>1
    zSOLgrid=zeros(settings.Box.nvar2,settings.Box.nvar1);
  else
    zSOLgrid=zeros(settings.Box.nvar1,1);
  end

  % pull out value at each time
  for ivar1=1:settings.Box.nvar1
    if settings.Box.nvar2>1
      for ivar2=1:settings.Box.nvar2
        [val,indx]=min(abs(Box(ivar1,ivar2).run_info.tMyrs-ttarg));
        zSOLgrid(ivar2,ivar1)=Box(ivar1,ivar2).run_info.zSOL(indx)/1000;
      end
    else
      [val,indx]=min(abs(Box(ivar1).run_info.tMyrs-ttarg));
      zSOLgrid(ivar1)=Box(ivar1).run_info.zSOL(indx)/1000;
    end
  end
end
