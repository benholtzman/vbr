function [posterior_T, best_T_g_phi] = fit_seismic_observations()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fit_seismic_observations
%
% Fit an input shear velocity profile (and LAB depth, if given) to state
% variables using VBR.
%
% Parameters:
% -----------
%       location    structure with the following required fields
%           lat         latitude [degrees North]
%           lon         longitude - assumed to be positive [degrees East]
%           z_min       minimum depth for observation range [km]
%           z_max       maximum depth for observation range [km]
%
%       Vs, Q, and LAB model file names
%           currently hardwired in in the calls to process_SeismicModels()
%
%       Name of a previously calculated parameter sweep (fname), or adjust
%       sweep_params       structure with the following required fields
%               T               vector of temperature values [deg C]
%               phi             vector of melt fractions [vol fraction]
%               gs              vector of grain sizes [micrometres]
%               per_bw_max      maximum period (min. freq.) considered [s]
%               per_bw_min      minimum period (max. freq.) considered [s]
%
% Output:
% -------
%      sc_posterior_S_given_Vs
%               matrix of posterior probabilities for all combinations of
%               sweep_params given the constraints from Vs observations
%
%      sc_posterior_S_given_Qinv
%               matrix of posterior probabilities for all combinations of
%               sweep_params given the constraints from Qinv observations
%
%      posterior_S_given_Vs_and_Qinv
%               matrix of posterior probabilities for all combinations of
%               sweep_params given the constraints from both Vs and Qinv
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup path structure
% put top level of VBR in the path
addpath('../../'); vbr_init

% put Project-specific paths in the path
addpath(genpath('./functions'))
buildProjectDirectories()
addpath(genpath('./'))

%% %%%%%%%%%%%%%%%% Get data for Vs(x, f) and Q(x, f) %%%%%%%%%%%%%%%%% %%
% Vs(x, f) and Q(x, f) are constrained by seismic observations and their %
% uncertainties.                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up location information
location.lat = 40; %lat; % degrees North
location.lon = 240; %lon; % degrees East
location.z_min = 100; % averaging min depth for asth.
location.z_max=150; % averaging max depth for asth.
location.smooth_rad = 5;


ifplot = 1;
vsfile = './data/vel_models/Shen_Ritzwoller_2016.mat';
[obs_Vs, sigma_Vs] = process_SeismicModels('Vs', ...
    location, vsfile, ifplot);

qfile = './data/Qinv_models/Gung_Romanowicz_2002.mat';
[obs_Qinv, sigma_Qinv] = process_SeismicModels('Qinv', ...
    location, qfile, ifplot);


%% %%%%%%%%%%%%%%%%%% Get prior for State Variables %%%%%%%%%%%%%%%%%%% %%
% The prior probability distribution for the state variables can be      %
% assumed to be either uniform or normal across the input range given.   %
% The probability that the given state variable is actually correct.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preferably, load in a large, pre-calculated box
fname = 'data/plate_VBR/sweep.mat';
if ~exist(fname, 'file')
    sweep_params.T = 1100:50:1700; %[degrees C]
    sweep_params.phi = (0.0:0.005:0.03); % melt fraction
    sweep_params.gs = linspace(0.001,0.03,10)*1e6; % grain size [micrometres]
    % Set period range for the mask - used to define which calculated
    % velocities go into the returned average Vs for those conditions
    sweep_params.per_bw_max = 30; % max period of range of mask (s)
    sweep_params.per_bw_min = 10; % min period of range of mask (s)

    sweep = generate_parameter_sweep(sweep_params);
    clear sweep_params
    save(fname, 'sweep')
end

load(fname, 'sweep');
% Extract the relevant values for the input depth range.
% Need to choose the attenuation method used for anelastic calculations
%       see possible methods by running vbrListMethods()
q_method = 'andrade_psp';
[sweep.meanVs, sweep.z_inds] = extract_calculated_values_in_depth_range(...
    sweep, 'Vs', q_method, [location.z_min, location.z_max]);
sweep.meanQinv = extract_calculated_values_in_depth_range(sweep, ...
    'Qinv', q_method, [location.z_min, location.z_max]);

% For each of the variables in sweep, set the mean and std
% Default is to calculate these based on the ranges set in sweep_params
params = make_param_grid(sweep.state_names, sweep);
% Note - can manually set the expected value and standard deviation for
% each of your variables, e.g. params.T_mean = 1500; params.gs_std = 300;

% Calculate the prior for either a normal or uniform distribution
pdf_type = {'uniform'};
[prior_statevars, sigma_statevars] = priorModelProbs( ...
    params, sweep.state_names, pdf_type);


%% %%%%%%%%%%%%%%%%%%%%% Get likelihood for Vs, Q %%%%%%%%%%%%%%%%%%%%%% %%
% The likelihood p(D|A), e.g., P(Vs | T, phi, gs), is calculated using    %
% the residual (See manual, Menke book Ch 11):                            %
%       p(D|A) = 1 / sqrt(2 * pi * residual) * exp(-residual / 2)         %
% residual(k) here is a chi-squared residual. Given chi-square, the PDF   %
% of data with a normal distribution:                                     %
%       P = 1 / sqrt(2 * pi * sigma^2) * exp(-0.5 * chi-square)           %
% where sigma = std of data, chi-square=sum((x_obs - x_preds)^2 / sigma^2)%
% e.g. www-cdf.fnal.gov/physics/statistics/recommendations/modeling.html  %
% The probability of getting the observed Vs or Q given the assumed state %
% variable values.                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


likelihood_Vs = probability_distributions('likelihood from residuals', ...
    obs_Vs, sigma_Vs, sweep.meanVs);

likelihood_Qinv = probability_distributions('likelihood from residuals', ...
    obs_Qinv, sigma_Qinv, sweep.meanQinv);

%% %%%%%%%%%%%%%%%% Get posterior for State Variables %%%%%%%%%%%%%%%%%% %%
% The posterior probability distribution is calculated in a Bayesian way  %
%       p(S | D) ? p(D | S) * p(S)                                        %
% The probability of the state variables given the observed Q and Vs.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

posterior_S_given_Vs = probability_distributions('A|B', ...
    likelihood_Vs, prior_statevars, 1);

posterior_S_given_Qinv =  probability_distributions('A|B', ...
    likelihood_Qinv, prior_statevars, 1);


vs_str = sprintf(['Vs = %.3g ', 177, ' %.2g km/s'], obs_Vs, sigma_Vs);
plot_Bayes(posterior_S_given_Vs, sweep, vs_str, q_method)
q_str = sprintf(['Qinv = %.2g ', 177, ' %.2g '], obs_Qinv, sigma_Qinv);
plot_Bayes(posterior_S_given_Qinv, sweep, q_str, q_method)
plot_tradeoffs_posterior(posterior_S_given_Vs, sweep, vs_str, q_method)
plot_tradeoffs_posterior(posterior_S_given_Qinv, sweep, q_str, q_method)

%% %%%%%%% Get posterior for State Variables given both Vs and Q %%%%%%% %%
% The measurement uncertainties in Vs and Q are assumed to be             %
% uncorrelated.  As such, we can simplify                                 %
%    p(S | (Vs, Qinv)) ? (p(Vs | S) * p(Qinv | S) * p(S))                 %
% The probability of the state variables being correct given constraints  %
% from both Vs and Q.                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

posterior_S_given_Vs_and_Qinv = probability_distributions(...
    'C|A,B conditionally independent', likelihood_Vs, likelihood_Qinv, ...
    prior_statevars, 1);

vs_q_str = [vs_str, ', ', q_str];
plot_Bayes(posterior_S_given_Vs_and_Qinv, sweep, vs_q_str, q_method)

plot_tradeoffs_posterior(posterior_S_given_Vs_and_Qinv, sweep, vs_q_str, q_method)


% Identify the best combinations of T, g, and phi across T
[best_T_phi_g, posterior_T] = find_best_state_var_combo( ...
    posterior_S_given_Vs_and_Qinv, sweep);


end

% 
% %% %%%%%%%%%%%%%%%%%%%% Get prior for LAB from RFs %%%%%%%%%%%%%%%%%%%%% %%
% % LAB depth is constrained by seismic observations and their              %
% % uncertainties.                                                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% labfile = './data/LAB_models/HopperFischer2018.mat';
% [obs_LAB, sigma_LAB] = process_SeismicModels('LAB_Depth', ...
%     location, labfile);
% 
% 
% %% %%%%%%%%%%%%%%%% Get prior for Tpot, zPlate %%%%%%%%%%%%%%%%%%%%%%%%% %%
% % The prior probability distribution for these variables can be assumed   %
% % to be either uniform or normal across the input range given.            %
% % The probability that the given variable is actually correct.            %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%% Generate geotherms and VBR boxes  %%%%%%%%
% Files.SV_Box='./data/plate_VBR/BigBox.mat';
% Files.VBR_Box='./data/plate_VBR/thermalEvolution_1.mat';
% defaults = init_settings;
% 
% recalc_SV = 0;
% 
% % Thermal model sweep for thermodynamic state variables
% if ~exist(Files.SV_Box,'file') || recalc_SV==1
%     LABsweep.TpC = 1350:25:1500;
%     LABsweep.zPlatekm = 60:20:160;
%     generate_boxes_ThermalEvolution(Files.SV_Box, LABsweep.TpC,...
%         LABsweep.zPlatekm); % builds T(z,t)
% else
%     Box = load(Files.SV_Box, 'Box');
%     LABsweep.TpC = Box.Box(1, 1).info.var1range;
%     LABsweep.zPlatekm = Box.Box(1, 1).info.var2range;
%     LABsweep.state_names = {'TpC', 'zPlatekm'};
%     fprintf(['\nSV Box already exists, delete following file to ' ...
%         'recalculate\n    %s\n'], Files.SV_Box)
% end
% 
% % Calculating VBR on BigBox takes < 1 minute
% VBRSettings.freq = logspace(-2.8,-1,4);
% VBRSettings.recalc_VBR=1;
% sweep.TpC = sweep.T ...
%     - defaults.dTdz_ad * mean([location.z_min, location.z_max]) * 1e3;
% best_Tp_phi_g = best_T_phi_g;
% best_Tp_phi_g(:, 1) = sweep.TpC; % Replace T with Tp
% 
% if ~exist(Files.VBR_Box,'file') || VBRSettings.recalc_VBR==1
%     process_ThermalEvolution_vbr(Files,VBRSettings.freq, best_Tp_phi_g);
% else
%     fprintf(['\nVBR Box already exists, delete following file to ' ...
%         'recalculate\n    %s\n'], Files.VBR_Box)
% end
% %%%%%%%%%%%%%%%%%% end VBR box calculations   %%%%%%%%%%%%%%%%%%
% 
% % Set period filter for calculating the LAB depth
% LAB_settings.per_bw_max = 30;  % max period to use for fitting (s)
% LAB_settings.per_bw_min = 10;  % min period to use for fitting (s)
% LAB_settings.q_method = q_method; % match Q method above
% predicted_vals = calc_LAB_Vs(Files.VBR_Box, LAB_settings);
% 
% 
% % For each of the variables in sweep, set the mean and std
% % Default is to calculate these based on the ranges set in sweep_params
% params = make_param_grid(LABsweep.state_names, LABsweep);
% % Replace default PDF for TpC with the posterior calculated from Vs and Q
% params.TpC_pdf = repmat(...
%     interp1(sweep.TpC, posterior_T, LABsweep.TpC)', ...
%     1, length(LABsweep.zPlatekm));
% 
% pdf_types = {'input', 'uniform'};
% [prior_vars, sigma_vars] = priorModelProbs( ...
%     params, LABsweep.state_names, {'input', 'uniform'});
% 
% %% %%%%%%%%%%%%%%%%%%%%%% Get likelihood for LAB %%%%%%%%%%%%%%%%%%%%%%% %%
% % The likelihood p(D|A), e.g., P(Vs | T, phi, gs), is calculated using    %
% % the residual (See manual, Menke book Ch 11):                            %
% %       p(D|A) = 1 / sqrt(2 * pi * residual) * exp(-residual / 2)         %
% % residual(k) here is a chi-squared residual. Given chi-square, the PDF   %
% % of data with a normal distribution:                                     %
% %       P = 1 / sqrt(2 * pi * sigma^2) * exp(-0.5 * chi-square)           %
% % where sigma = std of data, chi-square=sum((x_obs - x_preds)^2 / sigma^2)%
% % e.g. www-cdf.fnal.gov/physics/statistics/recommendations/modeling.html  %
% % The probability of getting the observed LAB given the assumed state     %
% % variable values.                                                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % chi^2 = sum( (observed - predicted)^2) / sigma ^ 2)
% likelihood_LAB = probability_distributions('likelihood from residuals', ...
%     obs_LAB, sigma_LAB, predicted_vals.zLAB_Q);
% 
% %% %%%%%%%%%%%%%%%%%%% Get posterior for Variables %%%%%%%%%%%%%%%%%%%%% %%
% % The posterior probability distribution is calculated in a Bayesian way  %
% %       p(V | D) = p(D | V) * p(V) / p(D)                                 %
% % The probability of the state variables given the observed LAB and the   %
% % observed Vs and Q informing the prior for the potential temperature.    %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% posterior_vars_given_LAB = probability_distributions('A|B', ...
%     likelihood_LAB, prior_vars, prior_LAB);
% plot_Bayes(posterior_vars_given_LAB, LABsweep, 'Vs, Qinv, LAB')
