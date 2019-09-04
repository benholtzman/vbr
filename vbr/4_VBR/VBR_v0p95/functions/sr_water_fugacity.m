%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             sr_water_fugacity.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output:
%        fH2O  water fugacity [MPa]
% input:
%        P_Pa     pressure [Pa]
%        T_K       temperature [K]
%        H2O_PPM   water concentration [PPM]
%
% Kolhstedt:
%  H2O = A_o * exp(-(E + P * V)/(R*T)) * fH2O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fH2O=sr_water_fugacity(H2O_PPM,H2O_o,P_Pa,T_K)
%   define constants
% CHeck which paper this comes from! 
    E = 40 * 1e3; % activation energy [J/mol]
    V = 10 * 1e-6; % activation volum [m3/mol]
    R = 8.314; % gas constant [J/mol/K]
    A_o = 26; % pre-expondential [PPM/MPa]

%   calculate fugacity
    fH2O = (H2O_PPM>=H2O_o).*(H2O_PPM/A_o).*exp((E+P_Pa*V)./(R*T_K));
end
