function generate_boxes(Work)
% This is time consuming as need to calculate the thermal evolution.
% Therefore, only do this step when sweeping through Tpot and zPlate
% (other variables do not impact thermal evolution).

clc
[calc_boxes_yn] = input('Do you want to calculate new boxes, y/[n]?  ','s');
if ~strcmp(calc_boxes_yn,'y'); return; end



% Parameter sweep
%  define parameter sweep here. var1name must match EXACTLY a field in
%  settings structure. var1 must be defined, var2 lines can be
%  commented/deleted if desired.
settings.Box.var1range = 1160:20:1600;
settings.Box.var1name = 'Tpot';
settings.Box.var1units =' C';

settings.Box.var2range = 60:10:250;
settings.Box.var2name = 'zPlate';
settings.Box.var2units =' km';


%  Overwrite any of the deafault settings (from init_settings) as desired
%    Mesh
settings.dz0=3; % grid cell size [km]
settings.Zinfo.asthenosphere_max_depth = 350; % adiabatic T from zmax to here [km]
settings.Z_moho_km = 30; % Moho depth [km]

%    Computational settings
%    time
settings.nt= 5000; % max number of time steps
settings.outk = settings.nt ; % frequency of output (output every outk steps)
% number of timesteps to save = outn = nt/outk
settings.t_max_Myrs=500; % max time to calculate [Myr]

%    for melt fraction calc
settings.sstol = 1e-5; % steady state target residual
settings.Flags.T_init='continental'; % 'continental' 'oceanic' or 'adiabatic'

master_DRIVE_SVFM_Plate(Work, settings)

end