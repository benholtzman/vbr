%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [VBR] = el_anharmonic(VBR)
% calculate anharomic moduli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [VBR] = el_anharmonic(VBR)
  [VBR] = el_calc_Gu_0(VBR); % set reference shear modulus 
  [VBR] = el_ModUnrlx_dTdP_f(VBR) ; % calculate anharmonic scaling with T, P of interest
end
