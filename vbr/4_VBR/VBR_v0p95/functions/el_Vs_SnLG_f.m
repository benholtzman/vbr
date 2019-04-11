%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [VBR] = el_Vs_SnLG_f(VBR)
% calculates Vs from Stixrude & Lithgow-Bertelloni parameterization
%
% paramterization takes T in K, P in GPa,
% anharmonic velocity Vs output in km/s
%
% reference:
% Stixrude and Lithgow‐Bertelloni (2005), "Mineralogy and elasticity of the oceanic
% upper mantle: Origin of the low‐velocity zone." JGR 110.B3,
% https://doi.org/10.1029/2004JB002965
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [VBR] = el_Vs_SnLG_f(VBR)
  dV_P=0.0380.*(VBR.in.SV.P_GPa);
  dV_T=-0.000378.*(VBR.in.SV.T_K-300);
  VBR.out.elastic.SLB2005.Vs = 4.77 +  dV_P + dV_T ;
end
