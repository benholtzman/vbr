function [Vp,Vs] = el_VpVs_unrelaxed(bulk_mod,shear_mod,rho)

  Vp = sqrt((bulk_mod + 4/3 * shear_mod)./rho);
  Vs = sqrt(shear_mod./rho); 
  
end
