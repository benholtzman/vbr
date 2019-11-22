%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  close all; clear

% put VBR in the path
  path_to_top_level_vbr='../../../';
  addpath(path_to_top_level_vbr)
  vbr_init

% experimental conditions
  stress=logspace(-1,log10(3e3),99); % [MPa]
  grains=logspace(0,4,100); % [um]
  [dg,sig]=meshgrid(grains,stress);
  P_MPa=[300 3*1e9/1e6];
  T_K = [1073 1673];

% construct state variable fields
  VBR.in.SV.sig_MPa = sig;
  VBR.in.SV.dg_um = dg;
  VBR.in.SV.phi = 0.001* ones(size(sig));
  VBR.in.SV.rho = 3300* ones(size(sig));

% write method list
  VBR.in.viscous.methods_list={'HK2003';'HZK2011'};
  VBR.in.GlobalSettings.melt_enhancement=0;

% calculate for the two experimental conditions
  VBR.in.SV.T_K = T_K(1) * ones(size(sig));
  VBR.in.SV.P_GPa = P_MPa(1) * 1e6/1e9 * ones(size(sig));
  [VBR] = VBR_spine(VBR);

  VBR.in.SV.T_K = T_K(2) * ones(size(sig));
  VBR.in.SV.P_GPa = P_MPa(2) * 1e6/1e9 * ones(size(sig));
  [VBR2] = VBR_spine(VBR);

% set up the contours / colors
  vals = [10 20 30];
  for iv = 1:2
    viscstud=VBR.in.viscous.methods_list{iv};

    visc=VBR.out.viscous;
    sr_tot = visc.(viscstud).sr_tot;
    sr_1 = visc.(viscstud).diff.sr;
    sr_2 = visc.(viscstud).disl.sr;
    sr_3 = visc.(viscstud).gbs.sr;

    clr = (sr_1./sr_tot > 0.5) * vals(1) ...
        + (sr_2./sr_tot > 0.5) * vals(2) ...
        + (sr_3./sr_tot > 0.5) * vals(3);

    VBR.out.viscous.(viscstud).clr = clr;

    visc=VBR2.out.viscous;
    sr_tot = visc.(viscstud).sr_tot;
    sr_1 = visc.(viscstud).diff.sr;
    sr_2 = visc.(viscstud).disl.sr;
    sr_3 = visc.(viscstud).gbs.sr;

    clr = (sr_1./sr_tot > 0.5) * vals(1) ...
        + (sr_2./sr_tot > 0.5) * vals(2) ...
        + (sr_3./sr_tot > 0.5) * vals(3);

    VBR2.out.viscous.(viscstud).clr = clr;

  end

% contour the results
  contlevs=3:-1:-20;

  figure('color',[1 1 1]')
  subplot(2,2,1)
  visc=VBR.out.viscous; viscstud='HK2003';
  cf=contourf((dg),(sig),(visc.(viscstud).clr),30,'linestyle','none');
  clims=caxis; caxis(clims)
  set(gca,'xscale','log','yscale','log','yminortick','on','xminortick','on')
  xlabel('grain size [{\mu}m]'); ylabel('deviatoric stress [MPa]')
  title(['strain rate [s^{-1}]'])% at ' num2str(P_MPa(1)) ' [MPa], ' num2str(T_K(1)) ' [K]'])
  hold on
  contour((dg),(sig),log10(visc.(viscstud).sr_tot),contlevs,'linewidth',1','linecolor','k')
  contour((dg),(sig),log10(visc.(viscstud).sr_tot),[-12 -12],'linewidth',2','linecolor','k','displayname','10^{-12} s^{-1}')

  subplot(2,2,3)
  visc=VBR.out.viscous; viscstud='HZK2011';
  cf=contourf((dg),(sig),(visc.(viscstud).clr),30,'linestyle','none');
  clims=caxis; caxis(clims)
  set(gca,'xscale','log','yscale','log','yminortick','on','xminortick','on')
  xlabel('grain size [{\mu}m]'); ylabel('deviatoric stress [MPa]')
  title(['strain rate [s^{-1}]'])% at ' num2str(P_MPa(1)) ' [MPa], ' num2str(T_K(1)) ' [K]'])
  hold on
  contour((dg),(sig),log10(visc.(viscstud).sr_tot),contlevs,'linewidth',1','linecolor','k')
  contour((dg),(sig),log10(visc.(viscstud).sr_tot),[-12 -12],'linewidth',2','linecolor','k','displayname','10^{-12} s^{-1}')

  subplot(2,2,2)
  visc=VBR2.out.viscous; viscstud='HK2003';
  cf=contourf((dg),(sig),(visc.(viscstud).clr),30,'linestyle','none');
  clims=caxis; caxis(clims)
  set(gca,'xscale','log','yscale','log','yminortick','on','xminortick','on')
  xlabel('grain size [{\mu}m]'); ylabel('deviatoric stress [MPa]')
  title(['strain rate [s^{-1}]'])% at ' num2str(P_MPa(1)) ' [MPa], ' num2str(T_K(1)) ' [K]'])
  hold on
  contour((dg),(sig),log10(visc.(viscstud).sr_tot),contlevs,'linewidth',1','linecolor','k')
  contour((dg),(sig),log10(visc.(viscstud).sr_tot),[-12 -12],'linewidth',2','linecolor','k','displayname','10^{-12} s^{-1}')

  subplot(2,2,4)
  visc=VBR2.out.viscous; viscstud='HZK2011';
  cf=contourf((dg),(sig),(visc.(viscstud).clr),30,'linestyle','none');
  clims=caxis; caxis(clims)
  set(gca,'xscale','log','yscale','log','yminortick','on','xminortick','on')
  xlabel('grain size [{\mu}m]'); ylabel('deviatoric stress [MPa]')
  title(['strain rate [s^{-1}]'])% at ' num2str(P_MPa(1)) ' [MPa], ' num2str(T_K(1)) ' [K]'])
  hold on
  contour((dg),(sig),log10(visc.(viscstud).sr_tot),contlevs,'linewidth',1','linecolor','k')
  contour((dg),(sig),log10(visc.(viscstud).sr_tot),[-12 -12],'linewidth',2','linecolor','k','displayname','10^{-12} s^{-1}')
