% PHYSICAL PROPERTIES
% State variables that would / will change
d0_um = 5e3 ; % microns
sig0_MPa = 0.1 ; 
phi_0 = 0.00 ;

% Density kg/m^3
Rho_olv = 3.3e3 ; % 
%Rho_crust = 3.0e3 ; % avg. feldspar density
%Rho_phi = 2.8e3; % melt (-- not needed now)

% thermal properties:
%Kc_olv = 3.3;		% Conductivity W/(m K)
%Kc_crust = 1.5;		% (source?) Conductivity W/(m K)

%cp_olv = 1100;		% Specific heat J /(kg K) PER UNIT MASS
%cp_crust = 800;		% (????) Specific heat J /(kg K)

% thermal diffusivity
%Kd_olv = Kc_olv/(cp_olv*Rho_olv) ;
%Kd_crust = Kc_crust/(cp_crust*Rho_crust) ; 
%Kd_Xu = 1.31e-6;       % Diffusivity at 298 K from Xu et al. 2004 	mm^2 / s

% Elastic properties (when is this used?)  
% unrelaxed elastic modulus
%Gu_crust = 44 ; %GPa
Gu_olv = 60 ; %GPa

% For scaling elastic properties ! 
% this is now in Params_VBR_General
%VBR.ref_val = 0 ; % 0 = Isaak, 1 = Cammarano
% this also affects parameters inside ModUnrlx_dTdP_f.m
%if VBR.ref_val == 0 ;
	% why this value, probably from Faul and Jackson ? 
%	Gu_olv = 65 ; % in GPa (but are converted later)
%elseif VBR.ref_val == 1
	% Cammarano has G_0 = 81 - 31*X_Fe, for X = 0.9 (Mg# 91), G_0 = 78.2
	% for reference temperature of 0 K or 300 K ? 
%	Gu_olv = 78.2 ; % in GPa (but are converted later)
%end
