function XtalF=calc_XtalFactor(z,settings)
  if strcmp(settings.Flags.XtalFactor,'constant')
     XtalF = settings.XtalFactor .* ones(size(z)); 
  elseif strcmp(settings.Flags.XtalFactor,'variable')
      z1=settings.XtalFactor_z1;
      z0=settings.XtalFactor_z0;
      Xf0=settings.XtalFactor; 
      XtalF=Xf0*(1-exp(-((z/1e3-z0)/z1).^2));
      XtalF(z/1e3<z0)=0;
      %      XtalF = 0.5*(1+erf((z/1e3 - z0)/z1));
  end
end