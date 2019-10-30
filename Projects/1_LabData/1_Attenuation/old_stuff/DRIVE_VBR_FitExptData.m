% DRIVER for simple VBR calculator in T.
% MAKE THIS MODERN !@#!@#!@#

clear all
addpath('./datasets')
VBR_version = 'VBR_v0p95'
addpath(['../../../vbr/4_VBR/',VBR_version ])
addpath(['../../../vbr/4_VBR/', VBR_version, '/functions'],['../../../vbr/4_VBR/', VBR_version, '/params'])
%addpath('./Plotting_VBR')
%global n

n_freq = 30;

% Read in the experimental data structures:
Make_DATA ;

%% MAKE VECTORS =======================================================
% SET UP THE MATRIX HERE...

Params_QMV_scaling; % in 4_VBR
%Params_VBR_general;
% VBR.ExptSet = ExptSet ;
% VBR.ref_val = 0 ; % determines reference Gu in Params_VBR_general

% ===========================
% VARIABLES

%% which set to run :
% 1 =  Jax :  Jackson & Tan et al, 2002 ?
% 2 =  Gribb and Cooper
% 3 =  SundCoop: Sunberg & Cooper 2010 (and Gribb and Cooper ? )
% 4 =  Takei + McCarthy
if VBR_version == 'VBR_v0p93','VBR_v0p91'
    VBR.FLAGS.Anelastic_calc = 'on'
end

for ExptSet = 1:3
    clearvars -except Data VBR ExptSet n_freq VBR_sols

if ExptSet == 1
    sz_dataset = size(Data.TanJax) ;
    n_Temps = sz_dataset(2) ;
    for nT = 1:n_Temps
        T_params(nT) = Data.TanJax(nT).exptCond.T_C ;
    end

	% X-axis = FREQUENCY
	f = (logspace(-3,1,n_freq))' ;
	% state vars etc (constant)
	P_GPa = Data.TanJax(1).exptCond.P_GPa ;
	rho = Data.TanJax(1).exptCond.rho ;
	phi = Data.TanJax(1).exptCond.phi_0 ;
	dg_um = Data.TanJax(1).exptCond.dg_0 ;
	sig_0 = Data.TanJax(1).exptCond.sig_0 ; % in Pa ! placeholder for calculating viscosity?
	% IF DIFFERENT FROM FAUL AND JACKSON !
	Gu_0 = Data.TanJax(1).exptCond.Gu ; % in GPa
    Ch2o_0 = 0.0 ;
    phi_0 = Data.TanJax(1).exptCond.phi_0 ;

elseif ExptSet == 2 % Gribb & Cooper
    sz_dataset = size(Data.GribCoop) ;
    n_Temps = sz_dataset(2) ;
    for nT = 1:n_Temps
        T_params(nT) = Data.GribCoop(nT).exptCond.T_C ;
    end

	% X-axis = FREQUENCY
	f = (logspace(-3,1,n_freq))' ;
	% state vars etc (constant)
	P_GPa = Data.GribCoop(1).exptCond.P_GPa ;
	rho = Data.GribCoop(1).exptCond.rho ;
	phi = Data.GribCoop(1).exptCond.phi_0 ;
	dg_um = Data.GribCoop(1).exptCond.dg_0 ;
	sig_0 = Data.GribCoop(1).exptCond.sig_0 ; % in Pa ! placeholder for calculating viscosity?
	% IF DIFFERENT FROM FAUL AND JACKSON !
	Gu_0 = Data.GribCoop(1).exptCond.Gu ; % in GPa
    Ch2o_0 = 0.0 ;
    phi_0 = Data.GribCoop(1).exptCond.phi_0 ;

elseif ExptSet == 3 % Sundberg & Cooper
    sz_dataset = size(Data.SundCoop) ;
    n_Temps = sz_dataset(2) ;
    for nT = 1:n_Temps
        T_params(nT) = Data.SundCoop(nT).exptCond.T_C ;
    end

	% X-axis = FREQUENCY
	f = (logspace(-3,1,n_freq))' ;
	% state vars etc (constant)
	P_GPa = Data.SundCoop(1).exptCond.P_GPa ;
	rho = Data.SundCoop(1).exptCond.rho ;
	phi = Data.SundCoop(1).exptCond.phi_0 ;
	dg_um = Data.SundCoop(1).exptCond.dg_0 ;
	sig_0 = Data.SundCoop(1).exptCond.sig_0 ; % in Pa ! placeholder for calculating viscosity?
	% IF DIFFERENT FROM FAUL AND JACKSON !
	Gu_0 = Data.SundCoop(1).exptCond.Gu ; % in GPa
    Ch2o_0 = 0.0 ;
    phi_0 = Data.SundCoop(1).exptCond.phi_0 ;
end

%% ====================================================
%% DRIVE THE VBR CALCULATOR ===========================
%% ====================================================


for jT = 1:n_Temps

	T_c = T_params(jT) ;
	T_K = T_c+273 ;
	make_SVs ;
	VBR.ISV.f = f ;

    [VBR] = VBR_spine(VBR) ;

    VBR_sols(ExptSet,jT).VBR = VBR ;
    VBR_sols(ExptSet,jT).T_c = T_c ;
    VBR_sols(ExptSet,jT).T_params = T_params;

end

end


%% SAVE and PLOT =================================================

save VBR_sols VBR_sols
