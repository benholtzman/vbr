% This code will fit an input shear velocity profile (and LAB depth, if
% given) to state variables using VBR

%basedir = 'C:\Users\Emily\Documents\VBR\';
basedir = 'D:\vbr\';
addpath([basedir 'vbr/6_FitVobs/Functions_Emily']);

Work = define_directories(basedir ,'2018-11-29-box');

%%%% GENERATE BOXES %%%%
generate_boxes(Work);
with_melt = 0; run_vbr(Work, with_melt);

%%%% GET SEISMIC OBSERVABLES %%%%
seismic_obs = get_seismic_data(Work);

%%%%% FIT PLATE THICKNESS %%%%
seismic_obs.q_method = 'eBurgers'; % 'AndradePsP'; 'YT_maxwell'; 'eBurgers';
zPlate = fit_LAB_Tp(Work, seismic_obs);

%%
%%%% MULTIVARIATE SEARCH - TEMP, GRAIN SIZE, PHI %%%%
% Bayesian inversion
[P_mod] = multivariate_Bayesian_inversion(...
    Work, seismic_obs, zPlate);

