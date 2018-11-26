% This code will fit an input shear velocity profile (and LAB depth, if
% given) to state variables using VBR

Work = define_directories(...
    'C:\Users\Emily\Documents\VBR\','2018-07-02-CHbox');

% GENERATE BOXES
generate_boxes(Work);
with_melt = 0; run_vbr(Work, with_melt);

% GET SEISMIC OBSERVABLES
seismic_obs = get_seismic_data(Work);

% FIT LAB DEPTH 
Work.q_method = 'eBurgers'; % 'AndradePsP'; 'YT_maxwell'; 'eBurgers';
zPlate = fit_LAB_Tp(Work, seismic_obs);


% FIT ASTHENOSPHERIC Vs


%   %   %   %   %   %   %   %   %   %
%  POTENTIAL TEMP/PLATE THICKNESS   %
%   %   %   %   %   %   %   %   %   %


%   %   %   %   %   %   %   %   %   %   %   %   %
%  MULTIVARIATE SEARCH - TEMP, GRAIN SIZE, PHI  %
%   %   %   %   %   %   %   %   %   %   %   %   %
% Bayesian inversion
