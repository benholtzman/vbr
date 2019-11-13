%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fit_plate
%
% Fit a LAB depth and estimate of asthenospheric temperature to a plate
% structure - defined by potential temperature (Tp) and the thermal
% lithospheric thickness, or depth of the intersection of the conductive
% geotherm with the adiabat (zPlate).
%
% Asthenospheric temperature cosntrains Tp only, so p(T) can be translated
% into p(Tp).  zLAB (estimates of the seismic LAB depth, e.g. from receiver
% functions) is controlled by both Tp and zPlate, though zPlate has a
% stronger influence.
%
% Given completely independent estimates of p(zLAB) and p(Tp), and a
% unique estimate of zPlate for every (zLAB, Tp) combination, we can just
% use Monte Carlo methods to approximate p(zPlate).
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


%% %%%%%%%%%%%%%%%%%%%%%% Load 1D thermal models %%%%%%%%%%%%%%%%%%%%%%% %%
% Equilibrium geotherms are calculated according to the parameters Tp and %
% zPlate.  This is relatively slow, so it is recommended to pre-calculate %
% a large box, save it, and re-use it.  The file name for this box is     %
% given by Files.SV_Box.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Generate geotherms and VBR boxes  %%%%%%%%
Files.SV_Box='./data/plate_VBR/BigBox.mat';
defaults = init_settings;

recalc_SV = 0;

% Thermal model sweep for thermodynamic state variables
if ~exist(Files.SV_Box,'file') || recalc_SV==1
    LABsweep.TpC = 1350:25:1500;
    LABsweep.zPlatekm = 60:20:160;
    generate_boxes_ThermalEvolution(Files.SV_Box, LABsweep.TpC,...
        LABsweep.zPlatekm); % builds T(z,t)
else
    Box = load(Files.SV_Box, 'Box');
    LABsweep.TpC = Box.Box(1, 1).info.var1range;
    LABsweep.zPlatekm = Box.Box(1, 1).info.var2range;
    LABsweep.state_names = {'TpC', 'zPlatekm'};
    fprintf(['\nSV Box already exists, delete following file to ' ...
        'recalculate\n    %s\n'], Files.SV_Box)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% Get prior for Tp %%%%%%%%%%%%%%%%%%%%%%%%%% %%
% p(Tp) can be set to whatever distribution is best constrained for your  %
% data.  Here, we assume that you have a posterior p(asthenospheric T)    %
% by fitting asthenospheric measurements of Vs (and/or Q).  This can be   %
% straightforwardly mapped to p(Tp) by converting T to Tp assuming an     %
% adiabatic gradient.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set range for Potential Temperature, Tp
sweep.TpC = sweep.T ...
    - defaults.dTdz_ad * mean([location.z_min, location.z_max]) * 1e3;

% Set the phi and g to use for every Tp in the VBR calculations
Tp_phi_g = [sweep.TpC', best_T_phi_g(:, 2:end)];


% Interpolate p(T) as a function of the asthenospheric sweep in T
% (sweep.TpC) into the sweep of TpC used for the geotherm modelling
% (LABsweep.TpC).
Tp_vals = linspace(min(LABsweep.TpC), max(LABsweep.TpC), 100);
p_Tp = interp1(sweep.TpC, posterior_T, Tp_vals);



%% %%%%%%%%%%%%%% Calculate zLAB for Tpot, zPlate %%%%%%%%%%%%%%%%%%%%%% %%
% The prior probability distribution for these variables can be assumed   %
% to be either uniform or normal across the input range given.            %
% The probability that the given variable is actually correct.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating VBR on 360 plate models takes ~ 30 s, but you can avoid
% recalculating if you set VBRSettings.recalc_VBR = 0 and give the file
% name of a previously calculated VBR box (Files.VBR_Box).
VBRSettings.recalc_VBR=1;
Files.VBR_Box='./data/plate_VBR/thermalEvolution_1.mat';
VBRSettings.freq = logspace(-2.8,-1,4);

if ~exist(Files.VBR_Box,'file') || VBRSettings.recalc_VBR==1
    process_ThermalEvolution_vbr(Files,VBRSettings.freq, Tp_phi_g);
else
    fprintf(['\nVBR Box already exists, delete following file to ' ...
        'recalculate\n    %s\n'], Files.VBR_Box)
end
%%%% end VBR box calculations   %%%

% Calculating the LAB depth also takes about 30 s
LAB_settings.per_bw_max = 30;  % max period to use for fitting (s)
LAB_settings.per_bw_min = 10;  % min period to use for fitting (s)
LAB_settings.q_method = q_method; % match Q method above
predicted_vals = calc_LAB_Vs(Files.VBR_Box, LAB_settings);

[zPlate_grid, Tp_grid] = meshgrid(LABsweep.zPlatekm, LABsweep.TpC);
zLAB_grid = predicted_vals.zLAB_Q;



%% %%%%%%%%%%%%%%%%%%%% Get prior for LAB from RFs %%%%%%%%%%%%%%%%%%%%% %%
% LAB depth is constrained by seismic observations and their              %
% uncertainties.                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labfile = './data/LAB_models/HopperFischer2018.mat';
[obs_LAB, sigma_LAB] = process_SeismicModels('LAB_Depth', ...
    location, labfile);

zLAB_vals = linspace(min(predicted_vals.zLAB_Q(:)), ...
    max(predicted_vals.zLAB_Q(:)), 100);
p_zLAB = probability_distributions('normal', zLAB_vals, obs_LAB, sigma_LAB);


%% %%%%%%%%%%%%%%%%%%%%%%%% Monte Carlo Time! %%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Randomly select values of Tp and zPlate from their respective pdfs      %
% and find the appropriate zPlate to construct p(zPlate).                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Best way to select a value from an arbitrary (and unnormalised)
% probability distribution is to construct the cdf and then randomly pick
% values between 0 and the largest value of the cdf (1 if pdf normalised)

cdf_zLAB = cumsum(p_zLAB * diff(zLAB_vals(1:2)));
cdf_Tp = cumsum(p_Tp * diff(Tp_vals(1:2)));

F = scatteredInterpolant(Tp_grid(:), zLAB_grid(:), zPlate_grid(:));

n_MC_trials = 1e5;
res_zPlate = zeros(1, n_MC_trials);
r_Tp = rand(1, n_MC_trials) .* max(cdf_Tp);
r_zLAB = rand(1, n_MC_trials) .* max(cdf_zLAB);
for n = 1:n_MC_trials
    [~, i_Tp] = min(abs(cdf_Tp - r_Tp(n)));
    [~, i_zLAB] = min(abs(cdf_zLAB - r_zLAB(n)));
    res_zPlate(n) = F(Tp_vals(i_Tp), LAB_vals(i_zLAB));
    
end
