function VBR = visc_calc_YT2016_solidus(VBR)
  % calculates the viscosity from Hatsuki Yamauchi and Yasuko Takei, JGR 2016,
  % "Polycrystal anelasticity at near-solidus temperatures,"

  % parameter checks
  if ~isfield(VBR.in,'viscous')
    VBR.in.viscous=struct();
  end
  if isfield(VBR.in.viscous,'YT2016_solidus')==0
    params=Params_Viscous('YT2016_solidus');
  else
    params=VBR.in.viscous.YT2016_solidus;
  end

  % calculate solidus-depdendent activation energy factor & melt effects
  Tprime=VBR.in.SV.T_K./VBR.in.SV.Tsolidus_K;
  A_n=calcA_n(Tprime,VBR.in.SV.phi,params);

  % calculate melt-free viscosity
  visc_method=params.eta_dry_method;
  if strcmp(visc_method,'YT2016_solidus')
    % use exactly what is in YT2016
    eta_dry=YT2016_dryViscosity(VBR,params);
  else
    % use a general olivine flow law to get melt-free diffusion-creep visc
    % need to re-run with phi=0 without losing other state variables
    VBRtemp=VBR;
    VBRtemp.in.viscous.methods_list={visc_method}; % only use one method    
    VBRtemp.in.SV.phi=0; % need melt-free viscosity
    VBRtemp=spineGeneralized(VBRtemp,'viscous');
    disp(fieldnames(VBRtemp.out))
    eta_dry = VBRtemp.out.viscous.(visc_method).diff.eta ;
  end

  % calculate full viscosity
  VBR.out.viscous.YT2016_solidus.diff.eta=A_n .* eta_dry;
  VBR.out.viscous.YT2016_solidus.diff.eta_dry=eta_dry;
end

function eta = YT2016_dryViscosity(VBR,params)
  % exact viscosity function in YT2016

  % reference values
    Tr=params.Tr_K; % ref. temp [K]
    Pr=params.Pr_Pa; % ref. pressure [Pa]
    eta_r=params.eta_r; % viscosity at (Tr,Pr,dr) [Pa s]

  % ref. grain size. per p7817 of reference paper (second paragraph). dr=d. i.e.,
  % gain size dependence is in eta_r.
    dr=VBR.in.SV.dg_um; % grain size independent beyond what is in etar, see

  % constants
    H=params.H; % activation energy [J/mol], figure 20 of reference paper
    Vol=params.V; % activation vol [m3/mol], figure 20 of reference paper
    R=params.R; % gas constant [J/mol/K]
    m=params.m; % grain size exponent -- but this does not matter since dr = d.

  % calculate the viscoscity
    P=VBR.in.SV.P_GPa*1e9;
    eta=eta_r.*(VBR.in.SV.dg_um./dr).^m  .* ...
          exp(Vol/R.*(P./VBR.in.SV.T_K-Pr./Tr)) .* ...
          exp(H/R.*(1./VBR.in.SV.T_K-1./Tr));

end

function A_n = calcA_n(Tn,phi,params)
  % calculates near-solidus and melt effects
  T_eta=params.T_eta;
  gamma=params.gamma;
  lambda=params.alpha; % rename to be consistent with YT2016 nomenclature

  A_n=zeros(size(Tn));

  A_n(Tn<T_eta)=1;

  msk=(Tn >= T_eta) & (Tn < 1);
  A_n(msk)=exp(-(Tn(msk)-T_eta)./(Tn(msk)-Tn(msk)*T_eta)*log(gamma));

  msk=(Tn > 1);
  A_n(msk)=exp(-lambda*phi(msk))/gamma;
end
