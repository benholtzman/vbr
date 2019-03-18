%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate reference shear modulus.
%
% uses a simple mixing model:
%
% VBR.in.SV.chi   composition (1 = olivine, 0 = crustal assemblage)
% VBR.in.elastic.anharmonic.Gu_0_ol  olivine reference shear modulus [GPa]
% VBR.in.elastic.anharmonic.Gu_0_crust crustal reference shear modulus [GPa]
%
% can be arrays of the same size or scalars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [VBR] = el_calc_Gu_0(VBR)
   Gu_0 = VBR.in.elastic.anharmonic.Gu_0_ol;  % reference moduls in GPa;

   chi = VBR.in.SV.chi; % compositional component
   Gc = VBR.in.elastic.anharmonic.Gu_0_crust;
   Gu = Gu_0 .* chi + (1-chi) .* Gc;

   VBR.out.elastic.Gu_0=1e9*Gu; % convert to Pa;
 end
