function fit_HK2003()
  close all;

  % put VBR in the path
  path_to_top_level_vbr='../../../';
  addpath(path_to_top_level_vbr)
  vbr_init

  expt_data=tryLoad(); % load data if it exists

  %%%%%%%%%%%%%%%%%%%%%%%
  %% Test the dry case %%
  %%%%%%%%%%%%%%%%%%%%%%%

  %% the thing to vary
  stress=logspace(1,3,55);

  %% construct state variable fields
  VBR.in.SV.T_K = (1250 + 273)*ones(size(stress));
  VBR.in.SV.P_GPa = 300 * 1e6/1e9*ones(size(stress));
  VBR.in.SV.sig_MPa = stress;
  VBR.in.SV.dg_um = 15*ones(size(stress));
  VBR.in.SV.phi = 0.0*ones(size(stress));
  VBR.in.SV.Ch2o = 0*ones(size(stress));

  %% write method list (these are the things to calculate)
  VBR.in.viscous.methods_list={'HK2003'};
  VBR.in.viscous.HK2003 = Params_Viscous('HK2003');
  VBR.in.viscous.HK2003.possible_mechs={'disl';'diff';};
  VBR.in.GlobalSettings.melt_enhancement=0;

  %% calculate viscosity
  [VBR] = VBR_spine(VBR);

  %% plot
  figure('color',[1 1 1])
  loglog(stress,VBR.out.viscous.HK2003.diff.sr,'--k','displayname','diff','linewidth',1.5)
  hold on
  loglog(stress,VBR.out.viscous.HK2003.disl.sr,'k','displayname','disl','linewidth',1.5)
  if isfield(VBR.out.viscous.HK2003,'gbs')
     loglog(stress,VBR.out.viscous.HK2003.gbs.sr,'--b','displayname','gbs','linewidth',1.5)
  end
  loglog(stress,VBR.out.viscous.HK2003.sr_tot,'k','linewidth',2,...
        'displayname','composite')

  if isfield(expt_data,'fig2a')
    loglog(expt_data.fig2a.sigma_MK2000,expt_data.fig2a.sr_MK2000,'ok','markersize',10)
    loglog(expt_data.fig2a.sigma_Disl,expt_data.fig2a.sr_Disl,'.k','markersize',10)
  end
  ylim([1e-7,1e-4])
  xlabel('\sigma [MPa]')
  ylabel('\epsilon [s^{-1}]')
  title('Dry: see figure 2A of Hirth and Kohlstedt 2003')
  legend('location','northwest')

  %%%%%%%%%%%%%%%%%%%%%%%
  %% Test the wet case %%
  %%%%%%%%%%%%%%%%%%%%%%%

  %% the thing to vary
  Ch2o = logspace(1,4,100);

  %% construct state variable fields
  VBR.in.SV.T_K = (1250 + 273)*ones(size(Ch2o));
  VBR.in.SV.P_GPa = 300 * 1e6/1e9 *ones(size(Ch2o));
  VBR.in.SV.sig_MPa = 150*ones(size(Ch2o));
  VBR.in.SV.dg_um = 15*ones(size(Ch2o));
  VBR.in.SV.phi = 0.0011*ones(size(Ch2o));
  VBR.in.SV.Ch2o = Ch2o;

  %% calculate viscosity
  [VBR] = VBR_spine(VBR);

  %% plot
  figure('color',[1 1 1])
  loglog(Ch2o,VBR.out.viscous.HK2003.sr_tot,'k','linewidth',2,...
        'displayname','composite')
  hold on
  loglog(Ch2o,VBR.out.viscous.HK2003.diff.sr,'--k','displayname','diff','linewidth',1.5)
  loglog(Ch2o,VBR.out.viscous.HK2003.disl.sr,'k','displayname','disl','linewidth',1.5)
  if isfield(VBR.out.viscous.HK2003,'gbs')
     loglog(Ch2o,VBR.out.viscous.HK2003.gbs.sr,'--b','displayname','gbs','linewidth',1.5)
  end
  xlabel('C_{H_2O} [PPM]')
  ylabel('\epsilon [s^{-1}]')
  title('Wet: see figure 5B of Hirth and Kohlstedt 2003')
  legend('location','northwest')
end

function data=tryLoad()
  dataDir='../../../../vbrWork/expt_data/2_viscosity/';

  data=struct();
  if exist([dataDir,'HK2003_fig2a_data_grab.mat'],'file')
    load([dataDir,'HK2003_fig2a_data_grab.mat'])
    data.fig2a.sigma_Disl=sigma_Disl;
    data.fig2a.sigma_MK2000=sigma_MK2000;
    data.fig2a.sr_MK2000=sr_MK2000;
    data.fig2a.sr_Disl=sr_Disl;
  end

  if exist([dataDir,'HK2003_fig5b_data_grab.mat'],'file')
    load([dataDir,'HK2003_fig5b_data_grab.mat'])
    data.fig5b.sr_HK=sr_HK;
    data.fig5b.Ch2o_HK=Ch2o_HK;
  end


end
