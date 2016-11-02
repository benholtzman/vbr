%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates Vs from Stixrude & Lithgow-Bertelloni parameterization
%
% paramterization takes T in K, P in GPa,
% anharmonic velocity Vs output in km/sec
%
%
% reference:
% Stixrude and Lithgow‐Bertelloni (2005), "Mineralogy and elasticity of the oceanic
% upper mantle: Origin of the low‐velocity zone." JGR 110.B3.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[VBR] = el_Vs_SnLG_f(VBR)

VBR.out.elastic.SLB2005.Vs = 4.77 + 0.0380.*(VBR.in.SV.P_GPa) ...
                                           - 0.000378.*(VBR.in.SV.T_K-300);  
end




