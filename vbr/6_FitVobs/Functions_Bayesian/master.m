% This code will fit an input shear velocity profile (and LAB depth, if
% given) to state variables using VBR

basedir = 'C:\Users\Emily\Documents\VBR\';
%basedir = 'D:\vbr\';
addpath([basedir 'vbr/vbr/6_FitVobs/Functions_Bayesian']);

Work = define_directories(basedir ,'2018-11-29-box');

%%%% GENERATE BOXES %%%%
generate_boxes(Work);
run_vbr(Work);
%
%%%% GET SEISMIC OBSERVABLES %%%%
seismic_obs = get_seismic_data(Work);
%

%%%%% FIT PLATE THICKNESS %%%%
seismic_obs.q_method = 'eBurgers'; % 'AndradePsP'; 'YT_maxwell'; 'eBurgers';
zPlate = fit_LAB_Tp(Work, seismic_obs, 1350);
%%
%%%% MULTIVARIATE SEARCH - TEMP, GRAIN SIZE, PHI %%%%
% Bayesian inversion
seismic_obs.q_method = 'eBurgers'; 
[P_mod] = multivariate_Bayesian_inversion(...
    Work, seismic_obs, zPlate, 1);



% Make it so you cacn run more easily with less keyboard early on
%   (consolidate into one input file)
% save output for playing
% calculate depth profiles for best N models
% play with N models
% different frequency bands
% make work for 1D models
% separate out getting seismic obs bc they SUUUUCCKK
