function buildComparisons(Box,figDir)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % buildComparisons(Box,figDir)
  %
  % compares anelastic methods, generates some figures
  %
  % Parameters
  % ----------
  % Box          the VBR box
  % figDir       directory to save figures to
  %
  % Output
  % ------
  % none         (figures to screen, written to file)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  close all

  plotBoxes(Box,figDir); % depth profiles
  plotComparisons(Box,figDir); % comparisons of anelastic methods

end


function plotBoxes(Box,figDir)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plotBoxes(Box,figDir)
  %
  % plots profiles for each box
  %
  % Parameters
  % ----------
  % Box          the VBR box
  % figDir       directory to save figures to
  %
  % Output
  % ------
  % none         (figures to screen, written to file)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ts=0:1:40;

  for iBox=1:numel(Box)
    B=Box(iBox);
    z=B.run_info.Z_km;

    fig=figure('Position', [10 10 700 400],'PaperPosition',[0,0,6,3],'PaperPositionMode','manual','DefaultAxesFontSize',8);
    ax_T=subplot(1,4,1);
    ax_phi=subplot(1,4,2);
    ax_eta=subplot(1,4,3);
    ax_V=subplot(1,4,4);

    set(fig,'currentaxes',ax_T)
    plot(B.VBR.in.SV.Tsolidus_K(:,1)-273,z,'k')

    tMyr=B.run_info.tMyrs;

    for t_indx=1:numel(ts)
      [val,indx]=min(abs(ts(t_indx)-tMyr));

      R=tMyr(indx) / max(ts);
      R=R*(R<=1)+1*(R>1);
      RGB=[R,0,1-R];

      % temperature v depth
      set(fig,'currentaxes',ax_T)
      hold on
      plot(B.VBR.in.SV.T_K(:,indx)-273,z,'color',RGB)

      % phi v depth
      set(fig,'currentaxes',ax_phi)
      hold on
      plot(B.VBR.in.SV.phi(:,indx),z,'color',RGB)

      % eta v depth
      set(fig,'currentaxes',ax_eta)
      hold on
      plot(log10(B.VBR.out.viscous.HK2003.eta_total(:,indx)),z,'color',RGB)

      % Vs averaged over methods and frequencies
      meths=fieldnames(B.VBR.out.anelastic);
      Vave=zeros(size(B.VBR.in.SV.phi(:,indx)));
      for imeth =1:numel(meths)
        meth=meths{imeth};
        Vave=Vave+B.VBR.out.anelastic.(meth).Vave(:,indx);
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
    title([B.info.var1name,'=',num2str(B.info.var1val),B.info.var1units])
    Tmajor=0:500:1500;
    set(ax_T,'xtick',Tmajor,'xticklabel',Tmajor)
    set(ax_T,'ydir','reverse')
    xlim([0,1700])
    box on


    saveas(gcf,[figDir,'/Box_',num2str(iBox),'_profiles.png'])
  end
end

function plotComparisons(Box,figDir)
  ts=0:.5:40;

  meths=fieldnames(Box(1).VBR.out.anelastic);
  ref_meth=meths{1};

  for iBox=1:numel(Box)
    figure()
    B=Box(iBox);
    t=B.run_info.tMyrs;
    z=B.run_info.Z_km;
    f=B.VBR.in.SV.f;
    % diff between methods at every depth, time, frequency
    dV=struct();
    dQinv=struct();
    for imeth =1:numel(meths)
      meth=meths{imeth};
      dV.(meth)=abs(B.VBR.out.anelastic.(meth).V-B.VBR.out.anelastic.(ref_meth).V)./B.VBR.out.anelastic.(ref_meth).V;
      dQinv.(meth)=abs(B.VBR.out.anelastic.(meth).Qinv-B.VBR.out.anelastic.(ref_meth).Qinv)./B.VBR.out.anelastic.(ref_meth).Qinv;
    end

    % loop over time, plot differences by method, frequency
    z_mask=(z > 30)&(z < 150);
    lnsty='-';
    NCs=numel(meths);
    clrs={'k','r','b','m','c'};
    ax_V=subplot(2,1,1);
    ax_Q=subplot(2,1,2);
    for ifreq=1:numel(f)

      for imeth =1:numel(meths)
        meth=meths{imeth};
        if ~strcmp(meth,ref_meth)
          dV_t=mean(dV.(meth)(z_mask,:,ifreq),1);
          dQinv_t=mean(dQinv.(meth)(z_mask,:,ifreq),1);

          lab=[meth,', ',num2str(f(ifreq))];

          set(gcf,'currentaxes',ax_V);
          hold all
          plot(t,dV_t,'DisplayName',lab,'LineStyle',lnsty,'color',clrs{imeth},'linewidth',2)
          box on
          ylabel(['dV (rel. to',ref_meth,')'])

          set(gcf,'currentaxes',ax_Q);
          hold all
          plot(t,dQinv_t,'DisplayName',lab,'LineStyle',lnsty,'color',clrs{imeth},'linewidth',2)

          box on
          ylabel(['dQ^{-1} (rel. to',ref_meth,')'])
          xlabel('t [Myrs]')

        end
      end
      lnsty='--';
    end

    % add the legend
    set(gcf,'currentaxes',ax_V);
    legend('location','eastoutside')
    title([B.info.var1name,'=',num2str(B.info.var1val),B.info.var1units])
    set(gcf,'currentaxes',ax_Q);
    legend('location','eastoutside')

    saveas(gcf,[figDir,'/Box_',num2str(iBox),'_comparison.png'])
  end
end
