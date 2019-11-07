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
  % VBRsettings  structure of some VBR settings
  %
  % Output
  % ------
  % none         (figures to screen, written to file)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  close all

  % for each Box, plot a thermal evolution and comparison panel
  plotBoxes(Box,figDir);

end


function plotBoxes(Box,figDir)


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


    saveas(gcf,[figDir,'/Box_',num2str(iBox),'.png'])
  end
end
