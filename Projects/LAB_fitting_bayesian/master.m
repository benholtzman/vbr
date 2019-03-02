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
  addpath(genpath('./'))

% %%%%%%%%%%%%%%%%%%   Generate Thermodynamic State Variables  %%%%%%%%%%%%%%%%%
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

% %%%%%%%%%%%%%%%%%%   Process State Variables with VBR  %%%%%%%%%%%%%%%%%
% %% VBR calcuation
  freq = logspace(-2.8,-1,4);
  g_um = 0.01 * 1e6; % 0.01 m grain size to micrometers
  if ~exist(Files.VBR_Box,'file') || recalc_VBR==1
    process_ThermalEvolution_vbr(Files,freq,g_um); % runs T(z,t) through VBR
  else
    disp("\nVBR Box already exists, delete following file to recalculate")
    disp(['    ',Files.VBR_Box,"\n"])
  end


% %%% GET SEISMIC OBSERVABLES %%%%
  Files.Vs_Model_file='./data/vel_models/Shen_Ritzwoller_2016.mat';
  Files.LAB_Model_file='./data/LAB_models/HopperFischer2018.mat';
  Coords.lat=30; Coords.lon=250; % lat/lon coordinate
  Coords.smooth_rad = 0.5;  % radius over which to average RFs and vel model (degrees)
  Coords.z_min=100; % averaging min depth for asth.
  Coords.z_max=150; % averaging max depth for asth.
  seismic_obs = process_SeismicModels(Files,Coords);


% %%%%% FIT PLATE THICKNESS %%%%
q_method = 'eBurgers'; % 'AndradePsP'; 'YT_maxwell'; 'eBurgers';
zPlate = fit_LAB_Tp(Files.VBR_Box, seismic_obs, 1350,q_method);
% %%
% %%%% MULTIVARIATE SEARCH - TEMP, GRAIN SIZE, PHI %%%%
% % Bayesian inversion
% seismic_obs.q_method = 'eBurgers';
% [P_mod] = multivariate_Bayesian_inversion(...
%     Work, seismic_obs, zPlate, 1);
