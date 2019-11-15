function [OutVBR,data]=FitData_YT16_visc()

  % put VBR in the path
  path_to_top_level_vbr='../../../';
  addpath(path_to_top_level_vbr)
  vbr_init
  addpath('./functions')

  data = loadYT2016visc();
  figure('DefaultAxesFontSize',12)
  OutVBR=struct();
  if isfield(data,'visc')
    N=numel(data.visc.sample_list);
    dg_range=max(data.visc.dg_um)-min(data.visc.dg_um);

    clrs={'k','r','b','c','m','g','p'}
    for isamp=1:N

      VBR.in=struct();
      VBR.in.viscous.methods_list={'xfit_premelt'};

      % pull this sample's data
      samp=data.visc.sample_list(isamp);
      dg=data.visc.dg_um(data.visc.sample==samp)(1);
      T_C=data.visc.T_C(data.visc.sample==samp);
      eta=data.visc.eta(data.visc.sample==samp);

      [T_C,I]=sort(T_C); eta=eta(I);
      VBR.in.viscous.xfit_premelt=setBorneolParams();


      VBR.in.SV.T_K = T_C+273 ;
      sz=size(VBR.in.SV.T_K);
      VBR.in.SV.dg_um= dg * ones(sz);
      VBR.in.SV.P_GPa = zeros(sz); % pressure [GPa]
      VBR.in.SV.phi = zeros(sz); % melt fraction
      VBR.in.SV.Tsolidus_K = 43.0 + 273 ;

      samp_field=['sample_',num2str(samp)];
      disp(['Calculating ',samp_field])
      VBR.in.viscous.xfit_premelt.dg_um_r=dg;
      VBR.in.viscous.xfit_premelt.Tr_K=T_C(1)+273;
      VBR.in.viscous.xfit_premelt.eta_r=eta(1);
      % fprintf([num2str(log10(params.eta_r)),',',num2str(params.Tr_K-273),',',num2str(params.dg_um_r),'\n'])

      % The pre-melting scaling takes into account the change in activation volume.
      % only want to use the lt 23 value
      VBR.in.viscous.xfit_premelt.H=data.table3_H.(samp_field).lt23.H*1e3;
      disp([samp_field,' lt: ',num2str(VBR.in.viscous.xfit_premelt.H)])
      [VBR_bysamp] = VBR_spine(VBR);
      VBReta=VBR_bysamp.out.viscous.xfit_premelt.diff.eta;
      OutVBR(isamp).sampleVBR=VBR_bysamp;

      VBR.in.viscous.xfit_premelt=setBorneolParams();
      [VBR] = VBR_spine(VBR);
      OutVBR(isamp).fullVBR=VBR;
      VBRetaGlob=VBR.out.viscous.xfit_premelt.diff.eta;

      clr=clrs{isamp};

      subplot(1,2,1)
      hold on
      semilogy(T_C,eta,'.','color',clr,'displayname',[num2str(dg),',',num2str(samp)],'MarkerSize',12)
      semilogy(T_C,VBReta,'color',clr,'displayname',[num2str(dg),',',num2str(samp)],'LineWidth',1.5)


      subplot(1,2,2)
      hold on
      semilogy(T_C,eta,'.','color',clr,'displayname',[num2str(dg),',',num2str(samp)],'MarkerSize',12)
      semilogy(T_C,VBRetaGlob,'color',clr,'displayname',[num2str(dg),',',num2str(samp)],'LineWidth',1.5)
    end
    subplot(1,2,1)
    title('H, dg\_ref, T\_ref, eta\_r set by sample')
    xlabel('T [C]'); ylabel('eta [Pa s]')
    legend('location','southwest')
    box on

    subplot(1,2,2)
    title('H=147 kJ/mol, dg\_ref=34.2 um, T\_ref=23 C, eta\_r=7e13 Pas')
    xlabel('T [C]'); ylabel('eta [Pa s]')
    box on
  else
    disp('This function requires data!')
  end

end
