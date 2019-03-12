function [VBR,telapsed]=spineViscosity(VBR)
  % calculates viscosity with user-selected methods
  telapsed.visc=tic;
  methods_list=VBR.in.viscous.methods_list; % list of methods to use

  if sum(strncmp('HK2003',methods_list,6)) > 0
     if isfield(VBR.in.viscous,'HK2003')==0
         VBR.in.viscous.HK2003=Params_Viscous('HK2003');
     end
     VBR = sr_visc_calc_HK2003(VBR);
  end

  if sum(strncmp('LH2012',methods_list,6)) > 0
     if isfield(VBR.in.viscous,'LH2012')==0
         VBR.in.viscous.LH2012=Params_Viscous('LH2012');
     end
     VBR = sr_visc_calc_LH2012(VBR);
  end

  if sum(strncmp('YT2016_solidus',methods_list,6)) > 0
     if isfield(VBR.in.viscous,'YT2016_solidus')==0
         VBR.in.viscous.YT2016_solidus=Params_Viscous('YT2016_solidus');
     end
     % this method only returns a single field, diffusion-creep viscosity.
     VBR = visc_calc_YT2016_solidus(VBR);
  end

  telapsed.visc=toc(telapsed.visc);
end
