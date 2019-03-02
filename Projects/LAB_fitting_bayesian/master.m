% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Fit an input shear velocity profile (and LAB depth, if given) to state
% variables using VBR.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% put VBR in the path
  path_to_top_level_vbr='../../';
  addpath(path_to_top_level_vbr)
  vbr_init

% put Project-specific paths in the path
  addpath(genpath('./functions'))
  buildProjectDirectories()
  addpath(genpath('./data'))

% %%%%%%%%%%%%%%%%%%   LAB FITTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% project-specific functions in ./functions
  Files.SV_Box='./data/plate_VBR/thermalEvolution.mat';
  Files.VBR_Box='./data/plate_VBR/thermalEvolution_VBR.mat';
  recalc_SV=0;
  recalc_VBR=0;

% %% Thermal model sweep for thermodynamic state variables
  Trange_C=[1300 1400];
  zPlateRange_km=[80 150];
  if ~exist(Files.SV_Box,'file') || recalc_SV==1
    generate_boxes_ThermalEvolution(Files.SV_Box,Trange_C,zPlateRange_km); % builds T(z,t)
  else
    disp("\nSV Box already exists, delete following file to recalculate")
    disp(['    ',Files.SV_Box,"\n"])
  end

% %% VBR calcuation
  freq = logspace(-2.8,-1,4);
  g_um = 0.01 * 1e6; % 0.01 m grain size to micrometers
  if ~exist(Files.VBR_Box,'file') || recalc_VBR==1
    process_ThermalEvolution_vbr(Files,freq,g_um); % runs T(z,t) through VBR
  else
    disp("\nVBR Box already exists, delete following file to recalculate")
    disp(['    ',Files.VBR_Box,"\n"])
  end


%%%% GET SEISMIC OBSERVABLES %%%%
% seismic_obs = get_seismic_data(Work);
% %
%
% %%%%% FIT PLATE THICKNESS %%%%
% seismic_obs.q_method = 'eBurgers'; % 'AndradePsP'; 'YT_maxwell'; 'eBurgers';
% zPlate = fit_LAB_Tp(Work, seismic_obs, 1350);
% %%
% %%%% MULTIVARIATE SEARCH - TEMP, GRAIN SIZE, PHI %%%%
% % Bayesian inversion
% seismic_obs.q_method = 'eBurgers';
% [P_mod] = multivariate_Bayesian_inversion(...
%     Work, seismic_obs, zPlate, 1);
