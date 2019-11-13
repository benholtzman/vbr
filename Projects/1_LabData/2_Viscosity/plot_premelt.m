function VBR=plot_premelt()
  close all; clear

  % put VBR in the path
  path_to_top_level_vbr='../../../';
  addpath(path_to_top_level_vbr)
  vbr_init

  VBR.in.viscous.methods_list={'xfit_premelt'};
  VBR.in.viscous.xfit_premelt=SetBorneolParams();

  % construct state variable fields
  VBR.in.SV.T_K= linspace(0,100,5)+273;
  VBR.in.SV.P_GPa = 0* ones(size(VBR.in.SV.T_K));
  VBR.in.SV.sig_MPa = 10* ones(size(VBR.in.SV.T_K));
  VBR.in.SV.dg_um = 4* ones(size(VBR.in.SV.T_K));
  VBR.in.SV.phi = 0.0 * ones(size(VBR.in.SV.T_K));
  VBR.in.SV.rho = 3300* ones(size(VBR.in.SV.T_K));
  VBR.in.SV.Tsolidus_K=80* ones(size(VBR.in.SV.T_K))+273;

  VBR=VBR_spine(VBR);

  % plot
  Tsol=[min(VBR.in.SV.Tsolidus_K(:)),max(VBR.in.SV.Tsolidus_K(:))]-273;
  eta=VBR.out.viscous.xfit_premelt.diff.eta;
  etaMM=[min(eta(:)),max(eta(:))];
  subplot(1,2,1)
  semilogy(VBR.in.SV.T_K-273,VBR.out.viscous.xfit_premelt.diff.eta,'k','linewidth',2,'DisplayName','eta')

  hold on
  semilogy(VBR.in.SV.T_K-273,VBR.out.viscous.xfit_premelt.diff.eta_meltfree,'--b','linewidth',2,'DisplayName','dry eta')
  semilogy(Tsol,etaMM,'--k','linewidth',2,'DisplayName','T_{sol}')
  xlabel('T [C]')
  ylabel('log10(eta)')
  legend('location','northeast')

  subplot(1,2,2)
  semilogy(VBR.in.SV.T_K ./ VBR.in.SV.Tsolidus_K,VBR.out.viscous.xfit_premelt.diff.eta,'k','linewidth',2)
  hold on
  semilogy(VBR.in.SV.T_K ./ VBR.in.SV.Tsolidus_K,VBR.out.viscous.xfit_premelt.diff.eta_meltfree,'--b','linewidth',2)
  semilogy([1 1],etaMM,'--k','linewidth',2)
  xlabel('T / Tsol')
  ylabel('log10(eta)')
end

function params=SetBorneolParams()
  % set the viscous parameters for borneol

  % near-solidus and melt effects
  params.alpha=0;
  params.T_eta=0.94;
  params.gamma=5;

  % flow law constants for YT2016
  params.Tr_K=100+273; % p7817 of YT2016, second paragraph
  params.Pr_Pa=0; % p7817 of YT2016, second paragraph
  params.eta_r=1e6; % figure 20 of reference paper
  params.H=100*1e3; % activation energy [J/mol], figure 20 of YT2016
  params.V=0; % activation vol [m3/mol], figure 20 of YT2016
  params.R=8.314; % gas constant [J/mol/K]
  params.m=3; % grain size exponent
  params.dg_um_r=4; % reference grain size [um]
end
