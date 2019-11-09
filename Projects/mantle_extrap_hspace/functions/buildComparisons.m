function buildComparisons(VBR,HS,figDir)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % buildComparisons(Box,figDir)
  %
  % compares anelastic methods, generates some figures
  %
  % Parameters
  % ----------
  % VBR          the VBR structure
  % HS           halfspace cooling structure
  % figDir       directory to save figures to
  %
  % Output
  % ------
  % none         (figures written to file)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % close all
  % plotProfiles(VBR,HS,figDir); % depth profiles
  % plotQprofiles(VBR,HS,figDir); % Q depth profiles
  plotPreMelt(VBR,HS,figDir); % premelting analysis
  % plotVsBoxes(VBR,HS,figDir); % Vs depth profiles
  % plotComparisons(VBR,HS,figDir); % comparisons of anelastic methods
  % close all
end

function plotPreMelt(VBR,HS,figDir)
  ifreq=1;
  istep=1;
  z=HS.z_km;

  figure()
  subplot(1,4,1)
  plot(VBR.in.SV.T_K(:,istep)-273,z,'k','Displayname','T')
  hold on
  plot(VBR.in.SV.Tsolidus_K(:,istep)-273,z,'r','Displayname','Tsol')
  xlabel('T [C]')

  subplot(1,4,2)
  plot(VBR.in.SV.T_K(:,istep)./VBR.in.SV.Tsolidus_K(:,istep),z,'k')
  hold on
  plot([0.9,0.9],[min(z),max(z)],'--k')
  plot([0.95,0.95],[min(z),max(z)],'--k')
  plot([1,1],[min(z),max(z)],'--k','linewidth',1.5)
  plot([1.1,1.1],[min(z),max(z)],'--k')
  xlim([0.6,1.2])
  xlabel('T/Tsol')

  subplot(1,4,3)
  plot(log10(VBR.out.anelastic.YT2016_solidus.J1(:,istep,ifreq)),z,'k','Displayname','J1')
  hold on
  plot(log10(VBR.out.anelastic.YT2016_solidus.J2(:,istep,ifreq)),z,'b','Displayname','J2')
  xlabel('log10(J)')
  legend('location','southeast')
  xlim([-16,-10])

  subplot(1,4,4)
  plot(log10(VBR.out.anelastic.YT2016_solidus.Q(:,istep,ifreq)),z,'k','Displayname','Q')
  xlabel('log10(Q)')
  xlim([0,6])

  for iplt=1:4
    subplot(1,4,iplt)
    set(gca,'ydir','reverse')
    ylim([25,150])
  end

end


function plotQprofiles(VBR,HS,figDir)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plotQprofiles(VBR,HS,figDir)
  %
  % plots Qinv profiles for each box
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  z=HS.z_km;

  fig=figure('Position', [10 10 700 400],'PaperPosition',[0,0,6,3],'PaperPositionMode','manual','DefaultAxesFontSize',8);

  meths=meths=fieldnames(VBR.out.anelastic);
  Nmeths=numel(meths);
  Nfreqs=numel(VBR.in.SV.f);
  for ifreq=1:Nfreqs
    for imeth=1:Nmeths
      iplt=imeth+Nmeths*(ifreq-1);
      ax_container(imeth,ifreq)=subplot(Nfreqs,Nmeths,iplt);
    end
  end

  tMyr=HS.t_Myr;
  for indx=1:numel(tMyr)

    R=tMyr(indx) / max(tMyr);
    R=R*(R<=1)+1*(R>1);
    RGB=[R,0,1-R];

    for ifreq=1:Nfreqs
      for imeth=1:Nmeths
        meth=meths{imeth};
        iplt=imeth+Nmeths*(ifreq-1);
        set(fig,'currentaxes',ax_container(imeth,ifreq))
        hold on
        plot(log10(VBR.out.anelastic.(meth).Q(:,indx,ifreq)),z,'color',RGB)
      end
    end
  end

  for ifreq=1:Nfreqs
    for imeth=1:Nmeths
      meth=meths{imeth};
      iplt=imeth+Nmeths*(ifreq-1);
      set(fig,'currentaxes',ax_container(imeth,ifreq))
      xlabel('log10(Q)')
      ylabel('depth [km]')
      title([strrep(meth,'_','\_'),' ',num2str(VBR.in.SV.f(ifreq))])
      set(gca,'ydir','reverse','box','on')
      xlim([1,6])
      ylim([0,200])
    end
  end

  saveas(gcf,[figDir,'Qprofiles.png'])

end

function plotVsBoxes(Box,figDir)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plotVsBoxes(Box,figDir)
  %
  % plots Vs profiles for each box, method, freq
  %
  % Parameters
  % ----------
  % VBR          the VBR structure
  % figDir       directory to save figures to
  %
  % Output
  % ------
  % none         (figures to screen, written to file)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % close all
  ts=0:1:40;

  for iBox=1:numel(Box)
    B=Box(iBox);
    z=run_info.Z_km;

    fig=figure('Position', [10 10 700 400],'PaperPosition',[0,0,6,3],'PaperPositionMode','manual','DefaultAxesFontSize',8);

    meths=meths=fieldnames(VBR.out.anelastic);
    Nmeths=numel(meths);
    Nfreqs=numel(VBR.in.SV.f);
    for ifreq=1:Nfreqs
      for imeth=1:Nmeths
        iplt=imeth+Nmeths*(ifreq-1);
        ax_container(imeth,ifreq)=subplot(Nfreqs,Nmeths,iplt);
      end
    end

    tMyr=run_info.tMyrs;
    for t_indx=1:numel(ts)
      [val,indx]=min(abs(ts(t_indx)-tMyr));

      R=tMyr(indx) / max(ts);
      R=R*(R<=1)+1*(R>1);
      RGB=[R,0,1-R];

      for ifreq=1:Nfreqs
        for imeth=1:Nmeths
          meth=meths{imeth};
          iplt=imeth+Nmeths*(ifreq-1);
          set(fig,'currentaxes',ax_container(imeth,ifreq))
          hold on
          plot((VBR.out.anelastic.(meth).V(:,indx,ifreq))/1e3,z,'color',RGB)
        end
      end
    end

    for ifreq=1:Nfreqs
      for imeth=1:Nmeths
        meth=meths{imeth};
        iplt=imeth+Nmeths*(ifreq-1);
        set(fig,'currentaxes',ax_container(imeth,ifreq))
        xlabel('Vs [km/s]')
        ylabel('depth [km]')
        title([strrep(meth,'_','\_'),' ',num2str(VBR.in.SV.f(ifreq))])
        set(gca,'ydir','reverse','box','on')
        xlim([4,5])
        ylim([0,200])
      end
    end

    saveas(gcf,[figDir,'/Box_',num2str(iBox),'_Vsprofiles.png'])
  end
end

function plotProfiles(VBR,HS,figDir)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plotProfiles(VBR,figDir)
  %
  % plots profiles for each box
  %
  % Parameters
  % ----------
  % VBR          the VBR structure
  % figDir       directory to save figures to
  %
  % Output
  % ------
  % none         (figures to screen, written to file)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % close all


  z=HS.z_km;
  tMyr=HS.t_Myr;

  fig=figure('Position', [10 10 700 400],'PaperPosition',[0,0,6,3],'PaperPositionMode','manual','DefaultAxesFontSize',8);
  ax_T=subplot(1,4,1);
  ax_phi=subplot(1,4,2);
  ax_eta=subplot(1,4,3);
  ax_V=subplot(1,4,4);

  set(fig,'currentaxes',ax_T)
  plot(VBR.in.SV.Tsolidus_K(:,1)-273,z,'k')


  for indx=1:numel(tMyr)

    R=tMyr(indx) / max(tMyr);
    R=R*(R<=1)+1*(R>1);
    RGB=[R,0,1-R];

    % temperature v depth
    set(fig,'currentaxes',ax_T)
    hold on
    plot(VBR.in.SV.T_K(:,indx)-273,z,'color',RGB)

    % phi v depth
    set(fig,'currentaxes',ax_phi)
    hold on
    plot(VBR.in.SV.phi(:,indx),z,'color',RGB)

    % eta v depth
    set(fig,'currentaxes',ax_eta)
    hold on
    plot(log10(VBR.out.viscous.HK2003.eta_total(:,indx)),z,'color',RGB)

    % Vs averaged over methods and frequencies
    meths=fieldnames(VBR.out.anelastic);
    Vave=zeros(size(VBR.in.SV.phi(:,indx)));
    for imeth =1:numel(meths)
      meth=meths{imeth};
      Vave=Vave+VBR.out.anelastic.(meth).Vave(:,indx);
    end
    Vave=Vave / numel(meths);

    set(fig,'currentaxes',ax_V)
    hold on
    plot(Vave/1e3,z,'color',RGB)

  end

  set(fig,'currentaxes',ax_phi)
  xlim([0,0.015])
  phiticks=[0,0.005,0.01];
  set(ax_phi,'xtick',phiticks,'xticklabel',phiticks,'yticklabel',{})
  xlabel('\phi')
  set(ax_phi,'ydir','reverse')
  box on

  set(fig,'currentaxes',ax_eta)
  set(ax_eta,'ydir','reverse','yticklabel',{})
  xlabel('log10(\eta) [Pa s]')
  xlim([17,26])
  box on

  set(fig,'currentaxes',ax_V)
  set(ax_V,'ydir','reverse','yticklabel',{})
  xlabel('V_s [km/s]')
  xlim([3.8,4.8])
  title('freq-method average')
  box on

  set(fig,'currentaxes',ax_T)
  xlabel('T [C]')
  ylabel('Depth [km]')
  Tmajor=0:500:1500;
  set(ax_T,'xtick',Tmajor,'xticklabel',Tmajor)
  set(ax_T,'ydir','reverse')
  xlim([0,1700])
  box on


  saveas(gcf,[figDir,'/HS_profiles.png'])

end

function plotComparisons(Box,figDir)
  ts=0:.5:40;
  % close all
  meths=fieldnames(Box(1).VBR.out.anelastic);
  ref_meth=meths{1};
  BadLines(numel(Box(1).VBR.in.SV.f)*numel(meths))=struct();

  for iBox=1:numel(Box)
    figure('Position', [10 10 700 400],'PaperPosition',[0,0,6,3],'PaperPositionMode','manual','DefaultAxesFontSize',8)
    B=Box(iBox);
    t=run_info.tMyrs;
    z=run_info.Z_km;
    f=VBR.in.SV.f;
    % diff between methods at every depth, time, frequency
    dV=struct();
    dQinv=struct();
    for imeth =1:numel(meths)
      meth=meths{imeth};
      dV.(meth)=VBR.out.anelastic.(meth).V/1e3;
      dQinv.(meth)=VBR.out.anelastic.(meth).Q;
    end

    % loop over time, plot differences by method, frequency
    z_mask=(z >= 45)&(z <= 55);
    lnsty='-';
    NCs=numel(meths);
    clrs={'k','r','b','m','c'};
    ax_V=subplot(3,3,[1,2]);
    ax_Q=subplot(3,3,[4,5]);
    ax_leg=subplot(3,3,[3,6]);
    ax_T=subplot(3,3,[7,8]);
    iBad=1;
    for ifreq=1:numel(f)

      for imeth =1:numel(meths)
        meth=meths{imeth};
        dV_t=mean(dV.(meth)(z_mask,:,ifreq),1);
        dQinv_t=mean(dQinv.(meth)(z_mask,:,ifreq),1);

        lab=strrep([meth,', ',num2str(f(ifreq))],'_','\_');

        set(gcf,'currentaxes',ax_V);
        hold all
        plot(t,dV_t,'DisplayName',lab,'LineStyle',lnsty,'color',clrs{imeth},'linewidth',2)

        set(gcf,'currentaxes',ax_Q);
        hold all
        plot(t,log10(dQinv_t),'DisplayName',lab,'LineStyle',lnsty,'color',clrs{imeth},'linewidth',2)

        set(gcf,'currentaxes',ax_leg);
        hold all
        BadLines(iBad).ln=plot(t,log10(dQinv_t),'DisplayName',lab,'LineStyle',lnsty,'color',clrs{imeth},'linewidth',2);
        iBad=iBad+1;

      end
      lnsty='--';
    end

    % add the legend
    set(gcf,'currentaxes',ax_V);
    title([info.var1name,'=',num2str(info.var1val),info.var1units])
    set(ax_V,'xticklabel',{});
    box on
    ylabel(['Vs [km/s]'])

    set(gcf,'currentaxes',ax_Q);
    ylabel('log10(Q)')
    box on

    set(gcf,'currentaxes',ax_leg);
    pos=get(ax_leg,'position');
    L=legend('location','north');
    set(L,'linewidth',0,'position',pos,'edgecolor',[1,1,1])
    set(ax_leg,'visible','off')
    for iBad=1:numel(BadLines)
      set(BadLines(iBad).ln,'visible','off')
    end

    set(gcf,'currentaxes',ax_T)
    T=mean(VBR.in.SV.T_K(z_mask,:)-273,1);
    Tsol=mean(VBR.in.SV.Tsolidus_K(z_mask,:)-273,1);
    plot(t,T./Tsol,'k')
    solparams = Params_Anelastic('YT2016_solidus');
    hold on
    for iAp=1:numel(solparams.Ap_Tn_pts)
      plot([t(1),t(end)],[solparams.Ap_Tn_pts(iAp),solparams.Ap_Tn_pts(iAp)],'--k')
    end

    % hold on
    % plot(t,Tsol,'r','DisplayName','Tsol')
    xlabel('t [Myrs]')
    box on
    ylabel('T / Tsol')


    saveas(gcf,[figDir,'/Box_',num2str(iBox),'_comparison.png'])
  end
end
