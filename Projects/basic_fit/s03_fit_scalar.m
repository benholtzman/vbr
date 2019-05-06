
  % (YS, TF) at ~80 km:(4.0598,3.9878)
  Obs.Vs=[4.0598,3.9878];
  Obs.Vs_name={'Yellowstone';'Twin Falls VF'};
  Obs.N=numel(Obs.Vs);

  plot_i_freq=1;
  plot_meth='eBurgers';

  VBR_file='./VBR_output.mat';
  load(VBR_file)

  meths={'eBurgers';'AndradePsP'};
  Fit=struct();


  for iObs =1:Obs.N
    Vs_Obs=Obs.Vs(iObs);
    for i_freq=1:numel(VBR.in.SV.f)

      for i_meth=1:numel(meths)
          meth=meths{i_meth};
          Vs=squeeze(VBR.out.anelastic.(meth).V(:,:,:,i_freq))/1e3;
          dV=abs(Vs-Vs_Obs)./Vs_Obs;
          Fit.(meth).dV(iObs,i_freq).dVgrid=dV;
          lin_el=find(dV == min(dV(:)));
          Fit.(meth).dV(iObs,i_freq).best.phi=VBR.in.SV.phi(lin_el);
          Fit.(meth).dV(iObs,i_freq).best.T_C=VBR.in.SV.T_K(lin_el)-273;
          Fit.(meth).dV(iObs,i_freq).best.dg_m=VBR.in.SV.dg_um(lin_el)/1e6;
          Fit.(meth).dV(iObs,i_freq).best.Vs=Vs(lin_el);
          Fit.(meth).dV(iObs,i_freq).best.dV=dV(lin_el);
          Fit.(meth).dV(iObs,i_freq).best.lin_el=lin_el;
      end
    end
  end


  meth=plot_meth;
  i_freq=plot_i_freq;
  phi=log10(squeeze(VBR.in.SV.phi(:,1,1)));
  T=(squeeze(VBR.in.SV.T_K(1,:,1))-273);
  dg=log10(squeeze(VBR.in.SV.dg_um(1,1,:))/1e6);

  for iObs =1:Obs.N
    Vs_Obs=Obs.Vs(iObs);

    figure
    disp(Obs.Vs_name{iObs})
    disp(Fit.(meth).dV(iObs,i_freq).best)
    dV=Fit.(meth).dV(iObs,i_freq).dVgrid;
    lin_el=Fit.(meth).dV(iObs,i_freq).best.lin_el;
    [best_phi,best_T,best_dg] = ind2sub(size(dV),lin_el);
    maxdV=-1;
    mindV=-6;

    subplot(1,3,1)
    period=1./VBR.in.SV.f(i_freq);
    imagesc(T,phi,log10(squeeze(dV(:,:,best_dg))));
    title_txt=[Obs.Vs_name{iObs},':',meth,', ',num2str(period),' s'];
    title(title_txt)
    xlabel('T [C]')
    ylabel('log10(phi)')
    colormap(hot)
    caxis([mindV,maxdV])

    subplot(1,3,2)
    imagesc(T,dg,log10(squeeze(dV(best_phi,:,:)))');
    xlabel('T [C]')
    ylabel('log10(dg) [m]')
    colormap(hot)
    caxis([mindV,maxdV])

    subplot(1,3,3)
    imagesc(phi,dg,log10(squeeze(dV(:,best_T,:))));
    xlabel('log10(phi)')
    ylabel('log10(dg) [m]')
    colormap(hot)
    caxis([mindV,maxdV])
    colorbar

  end
