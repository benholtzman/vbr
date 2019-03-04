% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Fit an input shear velocity profile (and LAB depth, if given) to state
% variables using VBR.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% put top level of VBR in the path
  addpath('../../'); vbr_init

% put Project-specific paths in the path
  addpath(genpath('./functions'))
  buildProjectDirectories()
  addpath(genpath('./'))

% %%%%%%%%%%%%%%%%%%   Generate Thermodynamic State Variables  %%%%%%%%%%%%%%%%%
% project-specific functions in ./functions
  Files.SV_Box='./data/plate_VBR/thermalEvolution.mat';
  Files.VBR_Box='./data/plate_VBR/thermalEvolution_VBR.mat';
  SVSettings.recalc_SV=0;
  SVSettings.Trange_C=[1300 1400];
  SVSettings.zPlateRange_km=[80 150];

% %% Thermal model sweep for thermodynamic state variables
  if ~exist(Files.SV_Box,'file') || SVSettings.recalc_SV==1
    generate_boxes_ThermalEvolution(Files.SV_Box,SVSettings.Trange_C,...
                                     SVSettings.zPlateRange_km); % builds T(z,t)
  else
    disp("\nSV Box already exists, delete following file to recalculate")
    disp(['    ',Files.SV_Box,"\n"])
  end

% %%%%%%%%%%%%%%%%%%   Process State Variables with VBR  %%%%%%%%%%%%%%%%%%%%%%%
  VBRSettings.freq = logspace(-2.8,-1,4);
  VBRSettings.g_um = 0.01 * 1e6; % 0.01 m grain size to micrometers
  VBRSettings.recalc_VBR=0;
  if ~exist(Files.VBR_Box,'file') || VBRSettings.recalc_VBR==1
    process_ThermalEvolution_vbr(Files,VBRSettings.freq,VBRSettings.g_um);
  else
    disp("\nVBR Box already exists, delete following file to recalculate")
    disp(['    ',Files.VBR_Box,"\n"])
  end

% %%%%%%%%%%%%%%%%%%   GET SEISMIC OBSERVABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Files.Vs_Model_file='./data/vel_models/Shen_Ritzwoller_2016.mat';
  Files.LAB_Model_file='./data/LAB_models/HopperFischer2018.mat';
  Coords.lat=30; Coords.lon=250; % lat/lon coordinate
  Coords.smooth_rad = 0.5;  % radius over which to average RFs and vel model (degrees)
  Coords.z_min=100; % averaging min depth for asth.
  Coords.z_max=150; % averaging max depth for asth.
  seismic_obs = process_SeismicModels(Files,Coords);

% %%%%%%%%%%%%%%%%%%   FIT PLATE THICKNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FitSettings.per_bw_max = 30;  % max period to use for fitting (s)
  FitSettings.per_bw_min = 10; % min period to use for fitting (s)
  FitSettings.set_Tp=1350; % potential temperature to pull best z_plate
  FitSettings.q_method='eBurgers'; % 'AndradePsP'; 'YT_maxwell'; 'eBurgers';
  zPlate = fit_LAB_Tp(Files.VBR_Box, seismic_obs, FitSettings);
  plotFits(Files.VBR_Box,zPlate,seismic_obs,FitSettings);

% %%%%%%%%%%%%%%%%%%   MULTIVARIATE SEARCH - TEMP, GRAIN SIZE, PHI %%%%%%%%%%%%%
% % Bayesian inversion
% seismic_obs.q_method = 'eBurgers';
% [P_mod] = multivariate_Bayesian_inversion(...
%     Work, seismic_obs, zPlate, 1);
