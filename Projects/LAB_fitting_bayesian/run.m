


location.lat = 40; %lat; % degrees North
location.lon = 240; %lon; % degrees East
location.z_min = 100; % averaging min depth for asth.
location.z_max=150; % averaging max depth for asth.
location.smooth_rad = 5;


filenames.Vs = './data/vel_models/Shen_Ritzwoller_2016.mat';
filenames.Q = './data/Q_models/Gung_Romanowicz_2002.mat';
filenames.LAB = './data/LAB_models/HopperFischer2018.mat';

% Extract the relevant values for the input depth range.
% Need to choose the attenuation method used for anelastic calculations
%       see possible methods by running vbrListMethods()
q_method = 'andrade_psp';


posterior_A = fit_seismic_observations(filenames, location, q_method);

posterior_L = fit_plate(filenames, location, q_method, posterior_A);


