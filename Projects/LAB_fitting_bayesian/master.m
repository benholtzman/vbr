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
  SVSettings.Trange_C=1300:25:1600;
  SVSettings.zPlateRange_km=80:10:150;

  % Thermal model sweep for thermodynamic state variables
  if ~exist(Files.SV_Box,'file') || SVSettings.recalc_SV==1
    generate_boxes_ThermalEvolution(Files.SV_Box,SVSettings.Trange_C,...
                                     SVSettings.zPlateRange_km); % builds T(z,t)
  else
    fprintf(['\nSV Box already exists, delete following file to recalculate\n' ...
        '    %s\n', Files.SV_Box])
  end

% %%%%%%%%%%%%%%%%%%   Process State Variables with VBR  %%%%%%%%%%%%%%%%%%%%%%%
  VBRSettings.freq = logspace(-2.8,-1,4);
  VBRSettings.g_um = 0.01 * 1e6; % 0.01 m grain size to micrometers
  VBRSettings.recalc_VBR=0;
  if ~exist(Files.VBR_Box,'file') || VBRSettings.recalc_VBR==1
    process_ThermalEvolution_vbr(Files,VBRSettings.freq,VBRSettings.g_um);
  else
    fprintf(['\nVBR Box already exists, delete following file to recalculate\n'...
        '    %s\n', Files.VBR_Box])
  end

% %%%%%%%%%%%%%%%%%%   GET SEISMIC OBSERVABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Files.Vs_Model_file='./data/vel_models/Shen_Ritzwoller_2016.mat';
  Files.LAB_Model_file='./data/LAB_models/HopperFischer2018.mat';
  Coords.lat=40; Coords.lon=245; % lat/lon coordinate
  Coords.smooth_rad = 0.5;  % radius over which to average RFs and vel model (degrees)
  Coords.z_min=100; % averaging min depth for asth.
  Coords.z_max=150; % averaging max depth for asth.
  seismic_obs = process_SeismicModels(Files,Coords);

% %%%%%%%%%%%%%%%%%%   FIT PLATE THICKNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FitSettings.per_bw_max = 30;  % max period to use for fitting (s)
  FitSettings.per_bw_min = 10;  % min period to use for fitting (s)
  FitSettings.set_Tp=1350; % potential temperature to pull best z_plate
  FitSettings.q_method='AndradePsP'; % 'AndradePsP'; 'MTH2011'; 'eBurgers';
  zPlate = fit_LAB_Tp(Files.VBR_Box, seismic_obs, FitSettings);
  plotFits(Files.VBR_Box,zPlate,seismic_obs,FitSettings);

% %%%%%%%%%%%%%%%%%%   MULTIVARIATE SEARCH - TEMP, GRAIN SIZE, PHI %%%%%%%%%%%%%
% Bayesian Inference -- maximum probability mapping

% define additional sweep parameters (Tpot taken from above VBR run)
  Files.VBR_bayesian='./data/plate_VBR/BayesianBox.mat';
  sweep_params.phi = (0.0:0.005:0.03); % melt fraction
  sweep_params.gs = linspace(0.001,0.03,10)*1e6; %[0.5 1 5:5:100].*1e3; % grain size
  sweep_params.q_method = FitSettings.q_method;
  sweep_params.per_bw_max = 30; % max period of range of mask (s)
  sweep_params.per_bw_min = 10; % min period of range of mask (s)
  P_mod = run_BayesianInference(Files,sweep_params,seismic_obs,zPlate,1);
  plot_Bayes_surf(P_mod)
