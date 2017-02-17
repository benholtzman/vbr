function SurfCoff = SurfaceTension(Vark,settings)
% Stenvenson formulation for driving force from surface tension: 
  SurfCoff = settings.y_sl./(2*Vark.grains*settings.mufo)...
             .*Vark.phis.^(-3/2).*Vark.perms;    
end