function VBRs = FitData_McCT11_all()
  % comparison of maxwell times, viscosities, unrelaxed moduli

  % put VBR in the path
  path_to_top_level_vbr='../../../';
  addpath(path_to_top_level_vbr)
  vbr_init
  addpath('./functions')

  ViscData = tryDataLoadVisc();
  if isfield(ViscData,'T_C')
    % Fig1=compare_viscosity(ViscData);
    Fig2=plot_relaxSpectrum();

    [VBRs] = calculate_Fits();
    Fig3 = plot_J1J2(VBRs);
  end

end

function [VBRs] = calculate_Fits()
  % sample 15, constant grain size 22 um
  Gu_40=2.4;
  Gu_20=2.5;
  dGdT=(Gu_40-Gu_20)/20*1e9;
  Tref=20;

  % anharmonic parameters
  VBR.in.elastic.methods_list={'anharmonic'};
  VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticit
  VBR.in.elastic.anharmonic.Gu_0_ol = Gu_40;
  VBR.in.elastic.anharmonic.dG_dT = dGdT;
  VBR.in.elastic.anharmonic.T_K_ref = 20+273;
  VBR.in.elastic.anharmonic.dG_dP = 0;

  % viscous parameters
  VBR.in.viscous.methods_list={'xfit_premelt'}; % far enough from solidus, A terms will go to 1
  VBR.in.viscous.xfit_premelt=SetBorneolParamsMcCT_lt7();

  % anelastic parameters
  VBR.in.anelastic.methods_list={'xfit_mxw'};
  % use fit 2! only diff is tau_n_p
  VBR.in.anelastic.xfit_mxw=Params_Anelastic('xfit_mxw');

  % set state variables
  VBR.in.SV.T_K = 23 +273 ; % temperature [K]
  VBR.in.SV.dg_um = 5; % grain size
  VBR.in.SV.P_GPa = 1000 / 1e9 ; % pressure [GPa] (does not affect this calc)
  VBR.in.SV.rho = 1000; % density [kg m^-3] (does not affect this calc)
  VBR.in.SV.sig_MPa = .1; % differential stress [MPa] (does not affect this calc)
  VBR.in.SV.phi = 0; % melt fraction
  VBR.in.SV.Tsolidus_K = (204.5 + 273) ; % solidus of pure borneol
  VBR.in.SV.f=logspace(-6,10,50);
  [VBR] = VBR_spine(VBR);
  VBRs.fit1.VBR=VBR;

  VBR.in.anelastic.xfit_mxw.fit='fit2';
  [VBR] = VBR_spine(VBR);
  VBRs.fit2.VBR=VBR;

end

function fig = plot_relaxSpectrum()

  fig = figure()

  tau_norm=logspace(-14,1,100);
  % plot the relexation spectrum
  params=Params_Anelastic('xfit_mxw');
  [X_tau] = Q_xfit_mxw_xfunc(tau_norm,params);
  hold all
  loglog(tau_norm,X_tau,'LineWidth',1.5,'r');

  params.fit='fit2';
  [X_tau] = Q_xfit_mxw_xfunc(tau_norm,params);
  loglog(tau_norm,X_tau,'LineWidth',1.5,'b');

  relaxData=tryDataLoadRelax();
  if isfield(relaxData,'relax_data_dg3to8')
    loglog(relaxData.relax_fit1_pts.tau_norm,relaxData.relax_fit1_pts.relax_fit1_pts,'.r')
    loglog(relaxData.relax_fit2_pts.tau_norm,relaxData.relax_fit2_pts.relax_fit2_pts,'.b')
    p1=relaxData.relax_PREM_pts.relax_PREM_pts(1);
    loglog(relaxData.relax_PREM_pts.tau_norm,[p1,p1],'k','linewidth',5)
    loglog(relaxData.relax_data_dg3to8.tau_norm,relaxData.relax_data_dg3to8.relax_data_dg3to8,'.k','MarkerSize',10)
  end

  xlabel('normalized time scale')
  ylabel('normalized relaxation spectrum')
  ylim([1e-4,2])
  xlim([1e-14,1e1])
  set(gca,'Xdir','reverse','XMinorTick','on','YMinorTick','on')
end


function fig = plot_J1J2(VBRs)


  fig = figure();
  ax_j1=subplot(1,2,1);
  ax_j2=subplot(1,2,2);

  set(gcf,'CurrentAxes',ax_j1)
  JU=1./VBRs.fit1.VBR.out.elastic.anharmonic.Gu;
  Fn=VBRs.fit1.VBR.out.anelastic.xfit_mxw.f_norm;
  semilogx(Fn,VBRs.fit1.VBR.out.anelastic.xfit_mxw.J1/JU,'r','linewidth',1.5)

  hold on;
  semilogx(Fn,VBRs.fit2.VBR.out.anelastic.xfit_mxw.J1/JU,'--r','linewidth',1.5)

  set(gcf,'CurrentAxes',ax_j2)
  loglog(Fn,VBRs.fit1.VBR.out.anelastic.xfit_mxw.J2/JU,'b','linewidth',1.5)
  hold on
  loglog(Fn,VBRs.fit2.VBR.out.anelastic.xfit_mxw.J2/JU,'--b','linewidth',1.5)

  data=tryDataLoadRelax();
  if isfield(data,'j1_data')

    set(gcf,'CurrentAxes',ax_j1)
    hold on;
    semilogx(data.j1_data.fnorm,data.j1_data.j1_data,'.k','MarkerSize',10)
    semilogx(data.j1_fit1.fnorm,data.j1_fit1.j1_fit1,'.r')
    semilogx(data.j1_fit2.fnorm,data.j1_fit2.j1_fit2,'.r')


    set(gcf,'CurrentAxes',ax_j2)
    hold on
    loglog(data.j2_data.fnorm,data.j2_data.j2_data,'.k','MarkerSize',10)
    loglog(data.j2_fit1.fnorm,data.j2_fit1.j2_fit1,'.b')
    loglog(data.j2_fit2.fnorm,data.j2_fit2.j2_fit2,'.b')
    % loglog(relaxData.relax_fit2_pts.tau_norm,relaxData.relax_fit2_pts.relax_fit2_pts,'.b')
    % p1=relaxData.relax_PREM_pts.relax_PREM_pts(1);
    % loglog(relaxData.relax_PREM_pts.tau_norm,[p1,p1],'k','linewidth',5)
    % loglog(relaxData.relax_data_dg3to8.tau_norm,relaxData.relax_data_dg3to8.relax_data_dg3to8,'.k','MarkerSize',10)
  end

  set(gcf,'CurrentAxes',ax_j1)
  xlabel('normalized frequency')
  ylabel('J1 / Ju')
  xlim([1e-1,1e11])
  ylim([1,3])
  set(gca,'XMinorTick','on','YMinorTick','on')

  set(gcf,'CurrentAxes',ax_j2)
  xlabel('normalized frequency')
  ylabel('J2 / Ju')
  xlim([1e-1,1e11])
  ylim([1e-4,2])
  set(gca,'XMinorTick','on','YMinorTick','on')


end


function fig = plot_EQ_fig9(data,ViscData,VBRs)

  include_E=0;

  fig=figure();

  if include_E==1
    ax_E=subplot(2,2,1);
    ax_E_norm=subplot(2,2,2);
    ax_Qinv=subplot(2,2,3);
    ax_Qinv_norm=subplot(2,2,4);
  else
    ax_Qinv=subplot(1,2,1);
    ax_Qinv_norm=subplot(1,2,2);
  end

  clrs={'k';'r';'b';'c';'m';'y';'g'};

  % plot the VBR
  for iT=1:numel(VBRs)

    VBR=VBRs(iT).VBR;
    This_T_C=VBR.in.SV.T_K-273;
    fld=['T_',num2str(round(This_T_C))];
    clr=clrs{iT};
    ClrStruct.(fld)=clr;

    fd=VBR.in.SV.f;
    M=VBR.out.anelastic.xfit_mxw.M/1e9;
    Mnorm=VBR.out.anelastic.xfit_mxw.M./VBR.out.elastic.anharmonic.Gu;
    Qinv=VBR.out.anelastic.xfit_mxw.Qinv;
    fnorm=VBR.out.anelastic.xfit_mxw.f_norm;

    if include_E==1
      set(fig,'CurrentAxes',ax_E)
      hold on
      semilogx(fd,M,clr);

      set(fig,'CurrentAxes',ax_E_norm)
      hold on
      semilogx(fnorm,Mnorm,clr);
    end

    set(fig,'CurrentAxes',ax_Qinv)
    hold on
    loglog(fd,Qinv,clr);

    set(fig,'CurrentAxes',ax_Qinv_norm)
    hold on
    loglog(fnorm,Qinv,clr);
  end

  % plot the data
  if isfield(data,'Qinv')
    for iT=1:numel(data.T_list)
      This_T_C=data.T_list(iT);

      fld=['T_',num2str(round(This_T_C))];
      clr=ClrStruct.(fld);


      fd=data.f_Hz(data.T_C==This_T_C);
      fd(fd<1e-4)=1e-4; % correction for the data grab
      M=data.E_GPa(data.T_C==This_T_C);
      Qinv=data.Qinv(data.T_C==This_T_C);

      iVisc=find(ViscData.T_C==This_T_C);
      mxwll=ViscData.tau_m_s(iVisc);
      fnorm=1./mxwll;
      Gfac=1./ViscData.GU_at_T_GPa(iVisc);

      if include_E==1
        set(fig,'CurrentAxes',ax_E)
        hold on
        semilogx(fd,M,['.',clr],'displayname','none','MarkerSize',10);

        set(fig,'CurrentAxes',ax_E_norm)
        hold on
        semilogx(fd/fnorm,M*Gfac,['.',clr],'displayname','none','MarkerSize',10);
      end

      set(fig,'CurrentAxes',ax_Qinv)
      hold on
      loglog(fd,Qinv,['.',clr],'MarkerSize',10);

      set(fig,'CurrentAxes',ax_Qinv_norm)
      hold on
      loglog(fd/fnorm,Qinv,['.',clr],'MarkerSize',10);
    end
  end

  if include_E==1
    set(fig,'CurrentAxes',ax_E)
    box on
    xlabel('f [Hz]'); ylabel('E [GPa]')
    ylim([.5,3])
    xlim([1e-4,10])
    set(gca,'XMinorTick','on')

    set(fig,'CurrentAxes',ax_E_norm)
    box on
    xlabel('f_N'); ylabel('E / E_U')
    ylim([0,1])
    xlim([1e-1,1e5])
    set(gca,'XMinorTick','on')
  end
  set(fig,'CurrentAxes',ax_Qinv)
  box on
  xlabel('f [Hz]'); ylabel('Qinv')
  xlim([1e-4,10])
  ylim([1e-2,2])
  set(gca,'XMinorTick','on','YMinorTick','on')

  set(fig,'CurrentAxes',ax_Qinv_norm)
  box on
  xlabel('f_N'); ylabel('Qinv')
  xlim([1e-1,1e5])
  ylim([1e-2,2])
  set(gca,'XMinorTick','on','YMinorTick','on')

end

function fig = compare_viscosity(data)
  % figure 8 of McCT11: viscosity vs grain size
  % pull the data
  eta=data.eta_a((data.T_C>=22.3)&(data.T_C<=23.7));
  dg=data.dg_um((data.T_C>=22.3)&(data.T_C<=23.7));
  dg=dg(eta>0);
  eta=eta(eta>0);

  etab=data.eta_b((data.T_C>=22.3)&(data.T_C<=23.7));
  dgb=data.dg_um((data.T_C>=22.3)&(data.T_C<=23.7));
  dgb=dgb(etab>0);
  etab=etab(etab>0);

  % VBR calc for dg_um > 7
  T_C=mean([22.3,23.7]);
  % viscous parameters
  VBR.in.viscous.methods_list={'xfit_premelt'}; %VBR.in.viscous.methods_list={'xfit_premelt'};
  VBR.in.viscous.xfit_premelt=SetBorneolParamsMcCT_gt7();
  VBR.in.SV.dg_um = logspace(0,2,100); % grain size
  VBR.in.SV.T_K = (T_C +273 ); % temperature [K]
  VBR.in.SV.P_GPa = 1000 / 1e9 ; % pressure [GPa] (does not affect this calc)
  VBR.in.SV.phi = 0; % melt fraction
  VBR.in.SV.Tsolidus_K = (204.5 + 273) ; % solidus of pure borneol
  [VBR_gt7] = VBR_spine(VBR) ;

  % VBR calc for dg_um < 7
  VBR.in.viscous.xfit_premelt=SetBorneolParamsMcCT_lt7();
  [VBR] = VBR_spine(VBR) ;
  fig=figure()
  loglog(dgb(dgb>=7),etab(dgb>=7),'.b','MarkerSize',10);
  hold on
  loglog(dgb(dgb<7),etab(dgb<7),'.r','MarkerSize',10);
  loglog(VBR_gt7.in.SV.dg_um,VBR_gt7.out.viscous.xfit_premelt.diff.eta,'b','LineWidth',1.5)
  loglog(VBR.in.SV.dg_um,VBR.out.viscous.xfit_premelt.diff.eta,'r','LineWidth',1.5)
  ylim([1e12,1e15])
  xlim([1e0,1e2])
  xlabel('Grain Size [um]')
  ylabel('Viscosity [Pa s]')
  title('T between 22.3, 23.7 C')
end


function data = tryDataLoadVisc()

  dataDir='../../../../vbrWork/expt_data/3_attenuation/McCT11/McCT11_new/';
  data=struct();
  if exist([dataDir,'McCT11_table1.csv'],'file')
    disp('loading')
    d=csvread([dataDir,'McCT11_table1.csv']);
    d=d(2:end,:);
    % sample	dg_um	T_C	eta_Pas_a	eta_Pas_b	eta_ave_Pas	tau_m_s	strain	GU_at_T_Gpa
    flds={'sample';'dg_um';'T_C';'eta_a';'eta_b';'eta_Pas';'tau_m_s';'strain';'GU_at_T_GPa'};
    for ifld=1:numel(flds)
      if ~strcmp(flds{ifld},'skip')
        data.(flds{ifld})=d(:,ifld);
      end
    end
    data.T_list=unique(data.T_C);
    data.dg_list=unique(data.dg_um);
    data.sample_list=unique(data.sample);
  else
    disp('not loading')
  end

end

function data = tryDataLoadRelax();
  data=struct();

  dataDir='../../../../vbrWork/expt_data/3_attenuation/McCT11/McCT11_new/';
  fi_list={'j1_data';'j1_fit1';'j1_fit2';'j2_data';'j2_fit1';'j2_fit2';...
           'relax_fit1_pts';'relax_fit2_pts';'relax_PREM_pts';'relax_data_dg3to8'};

  for ifi=1:numel(fi_list)
    fl=fi_list{ifi};
    if exist([dataDir,fl,'.csv'])
      d=csvread([dataDir,fl,'.csv']);
      if numel(strfind(fl,'relax')>0)
        data.(fl).tau_norm=d(:,1);
        data.(fl).(fl)=d(:,2);
      else
        data.(fl).fnorm=d(:,1);
        data.(fl).(fl)=d(:,2);
      end
    end
  end

end

function data = tryDataLoadFig9()
  % loads data from fig 9: Qinv, E vs freq at constant grain size, varying T

  dataDir='../../../../vbrWork/expt_data/3_attenuation/McCT11/McCT11_new/';
  data=struct();
  if exist([dataDir,'sample_15_Tdependence_fig9.csv'],'file')
    disp('loading')
    d=csvread([dataDir,'sample_15_Tdependence_fig9.csv']);
    d=d(2:end,:);
    data.T_C=d(:,1);
    data.f_Hz=d(:,4);
    data.Qinv=d(:,5);
    data.E_GPa=d(:,6);
    data.T_list=unique(data.T_C);
  else
    disp('not loading')
  end

end

function params = SetBorneolParamsMcCT_gt7()
  % set the viscous parameters for borneol
  % near-solidus and melt effects
  params.alpha=0;
  params.T_eta=0.94; % eqn 17,18- T at which homologous T for premelting.
  params.gamma=5;
  % flow law constants for YT2016
  params.Tr_K=22.5+273; %
  params.Pr_Pa=1000; % p7817 of YT2016, second paragraph
  params.eta_r=5.25*1e13; %
  params.H=85.4*1e3; % activation energy [J/mol], figure 20 of YT2016
  params.V=0; % activation vol [m3/mol], figure 20 of YT2016
  params.R=8.314; % gas constant [J/mol/K]
  params.m=1; % grain size exponent
  params.dg_um_r=21.4 ; % caption of Fig 9. % 24.4; % reference grain size [um]
end

function params = SetBorneolParamsMcCT_lt7()
  % set the viscous parameters for borneol
  % near-solidus and melt effects
  params.alpha=0;
  params.T_eta=0.94; % eqn 17,18- T at which homologous T for premelting.
  params.gamma=5;
  % flow law constants for YT2016
  params.Tr_K=23.6+273; %
  params.Pr_Pa=1000; % p7817 of YT2016, second paragraph
  params.eta_r=2.04*1e12; %
  params.H=85.4*1e3; % activation energy [J/mol], figure 20 of YT2016
  params.V=0; % activation vol [m3/mol], figure 20 of YT2016
  params.R=8.314; % gas constant [J/mol/K]
  params.m=3; % grain size exponent
  params.dg_um_r=3.35 ; % caption of Fig 9. % 24.4; % reference grain size [um]
end
