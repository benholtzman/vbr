% This code will fit an input shear velocity profile (and LAB depth, if
% given) to state variables using VBR

basedir = 'C:\Users\Emily\Documents\VBR\';
addpath([basedir 'vbr/vbr/6_FitVobs/Functions_Emily']);

Work = define_directories(basedir ,'2018-07-02-CHbox');

%%%% GENERATE BOXES %%%%
generate_boxes(Work);
with_melt = 0; run_vbr(Work, with_melt);

%%%% GET SEISMIC OBSERVABLES %%%%
seismic_obs = get_seismic_data(Work);

%%%%% FIT LAB DEPTH %%%%
q_method = 'eBurgers'; % 'AndradePsP'; 'YT_maxwell'; 'eBurgers';
zPlate = fit_LAB_Tp(Work, seismic_obs, q_method);


%%%% MULTIVARIATE SEARCH - TEMP, GRAIN SIZE, PHI %%%%
% Bayesian inversion
