%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Vp,Vs] = el_VpVs_unrelaxed(bulk_mod,shear_mod,rho)
% calculates unrelaxed Vp, Vs from bulk modulus, shear modulus, density 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vp,Vs] = el_VpVs_unrelaxed(bulk_mod,shear_mod,rho)
  Vp = sqrt((bulk_mod + 4/3 * shear_mod)./rho);
  Vs = sqrt(shear_mod./rho);
end
