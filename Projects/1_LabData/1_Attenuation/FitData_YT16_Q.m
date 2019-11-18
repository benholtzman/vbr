function [OutVBR,data]=FitData_YT16_Q()

  % put VBR in the path
  path_to_top_level_vbr='../../../';
  addpath(path_to_top_level_vbr)
  vbr_init
  addpath('./functions')

  % load data if it exists
  viscData = loadYT2016visc();
  Qdata=loadYT2016Q();


  figure('DefaultAxesFontSize',12)
  OutVBR=struct();
  if isfield(viscData,'visc') && isfield(Qdata,'Qinv')

    experimental_Ts=unique(Qdata.Qinv.T_C); % the temp conditions at which f was varied
    N=numel(experimental_Ts); % number of experimental P/T conditions
    clrs={'k','r','b','c','m','g','y'};
    samp=41;

    % loop over exp. conditions (T is varied only), calculate Q
    for iexp=1:N

      This_T_C=experimental_Ts(iexp);
      VBR.in=struct();
      VBR.in.elastic.methods_list={'anharmonic';'anh_poro'};
      VBR.in.viscous.methods_list={'xfit_premelt'};
      VBR.in.anelastic.methods_list={'xfit_premelt'};

      % pull this sample's viscData
      dg=viscData.visc.dg_um(viscData.visc.sample==samp)(1);
      T_Cvisc=viscData.visc.T_C(viscData.visc.sample==samp);
      eta=viscData.visc.eta(viscData.visc.sample==samp);
      [T_Cvisc,I]=sort(T_Cvisc); eta=eta(I);

      % set this sample's viscosity parameters
      VBR.in.viscous.xfit_premelt=setBorneolParams();
      VBR.in.viscous.xfit_premelt.dg_um_r=dg;
      VBR.in.viscous.xfit_premelt.Tr_K=T_Cvisc(1)+273;
      VBR.in.viscous.xfit_premelt.eta_r=eta(1);

      % set anharmonic conditions
      VBR.in.elastic.anharmonic=Params_Elastic('anharmonic');
      [Gu_o,dGdT,dGdT_ave]= YT16_E(This_T_C);
      Gu_o=Gu_o-0.05;

      % Gu_o is for a given T, set anharmonic derives to 0
      VBR.in.elastic.anharmonic.Gu_0_ol = Gu_o;
      VBR.in.elastic.anharmonic.dG_dT = 0;
      VBR.in.elastic.anharmonic.dG_dP = 0;
      disp(['Calculating for T=',num2str(This_T_C)])
      disp(sprintf('  Gu(T=%0.1f C)=%0.3f GPa',This_T_C,Gu_o))

      % adjust some anelastic settings
      VBR.in.anelastic.xfit_premelt.tau_pp=2*1e-5;

      % set experimental conditions
      VBR.in.SV.T_K = This_T_C+273 ;
      sz=size(VBR.in.SV.T_K);
      VBR.in.SV.dg_um= dg.* ones(sz);
      VBR.in.SV.P_GPa = 1.0132e-04 .* ones(sz); % pressure [GPa]
      VBR.in.SV.rho =1011 .* ones(sz); % density [kg m^-3]
      VBR.in.SV.sig_MPa =1000 .* ones(sz)./1e6; % differential stress [MPa]
      VBR.in.SV.phi = zeros(sz); % melt fraction
      VBR.in.SV.Tsolidus_K = 43.0 + 273 ;
      VBR.in.SV.Ch2o_0=zeros(sz);

      VBR.in.SV.f=logspace(-4,2,50);

      samp_field=['sample_',num2str(samp)];

      % The pre-melting scaling takes into account the change in activation volume.
      % only want to use the lt 23 value
      VBR.in.viscous.xfit_premelt.H=viscData.table3_H.(samp_field).lt23.H*1e3;

      [VBR_bysamp] = VBR_spine(VBR);
      VBR_Q_samp=VBR_bysamp.out.anelastic.xfit_premelt.Qinv;
      VBR_G_samp=VBR_bysamp.out.anelastic.xfit_premelt.M/1e9;
      OutVBR(iexp).sampleVBR=VBR_bysamp;

      Q_obs=Qdata.Qinv.Qinv(Qdata.Qinv.T_C==This_T_C);
      Q_obs_f=Qdata.Qinv.f(Qdata.Qinv.T_C==This_T_C);
      E_obs=Qdata.E.E(Qdata.E.T_C==This_T_C);
      E_obs_f=Qdata.E.f(Qdata.E.T_C==This_T_C);

      if iexp > numel(clrs)
        icolor=iexp-numel(clrs);
        lnsty='--';
      else
        icolor=iexp;
        lnsty='';
      end
      clr=clrs{icolor};

      subplot(2,1,2)
      hold on
      loglog(Q_obs_f,Q_obs,'.','color',clr,'displayname',[num2str(dg),',',num2str(samp)],'MarkerSize',12)
      loglog(VBR.in.SV.f,VBR_Q_samp,[lnsty,clr],'displayname',[num2str(dg),',',num2str(samp)],'LineWidth',1.5)

      subplot(2,1,1)
      hold on
      semilogx(E_obs_f,E_obs,'.','color',clr,'displayname',[num2str(dg),',',num2str(samp)],'MarkerSize',12)
      semilogx(VBR.in.SV.f,VBR_G_samp,[lnsty,clr],'displayname',[num2str(dg),',',num2str(samp)],'LineWidth',1.5)

    end
    subplot(2,1,2)
    xlabel('f [Hz]'); ylabel('Q^{-1}')
    box on

    subplot(2,1,1)
    xlabel('f [Hz]'); ylabel('M [GPa]')
    box on
  else
    disp('This function requires data!')
  end

end
