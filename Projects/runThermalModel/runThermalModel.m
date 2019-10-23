%% ===================================================================== %%
%%                     runThermalModel.m
%% ===================================================================== %%
%  example of running a single forward model
%% ===================================================================== %%

clear

% put VBR in the path
path_to_top_level_vbr='../../';
addpath(path_to_top_level_vbr)
vbr_init

% set some  parameters
settings.Cs0_H20=300;
settings.Tpot=1400;

% Mesh
settings.dz0=2; % grid cell size [km]
settings.Zinfo.zmax = 200; % max z depth to model [km]
settings.Zinfo.asthenosphere_max_depth = 300; % adiabatic T from zmax to here [km]
settings.Z_moho_km = 10; % Moho depth [km]

% Computational settings
settings.nt= 500; % max number of time steps
settings.outk = 10 ; % frequency of output (output every outk steps)
settings.t_max_Myrs=50; % max time to calculate [Myr]
settings.sstol = 1e-5; % steady state target residual

% initial temperature condition
settings.Flags.T_init='oceanic'; % 'continental' 'oceanic' or 'adiabatic'

% Load Default Settings for the forward model runs (preserves those set above)
[settings]=init_settings(settings);

% build the mesh and initial conditions
settings.Zinfo.zmax=settings.zPlate;
settings.Zinfo.dz0 = settings.dz0;
settings.Zinfo = init_mesh(settings.Zinfo); % build the mesh!
[Info] = init_values(settings); % calculate initial values

% initalize  the boundary conditions
[Info.BCs]=init_BCs(struct(),'T','zmin','dirichlet',0);
[Info.BCs]=init_BCs(Info.BCs,'T','zmax','dirichlet',Info.init.T(end));

% Lithosphere temperature evolution (the actual calculation)
[Vars,Info]=Thermal_Evolution(Info,settings);

% Add on adiabatic asthenosphere below model domain
[Vars,Info]=postproc_append_astheno(Vars,Info,settings);


disp(' ');disp('--------------------------------------------------');disp(' ')

figure()
subplot(1,2,1)
plot(Vars(:,end).Tsol,Info.z_km,'--k')
hold on
for it = 1:2: numel(Info.tMyrs)
  cf=(it - 1) / (numel(Info.tMyrs)-1);
  plot(Vars.T(:,it),Info.z_km,'color',[0,0,cf])
end
box on
set(gca,'ydir','reverse')
xlabel('geotherms and solidii [C]')
ylabel('depth [km]')

subplot(1,2,2)
plot(Info.tMyrs,Info.zSOL/1000,'k','linewidth',2)
hold on
plot(Info.tMyrs,Info.zLAB/1000,'--k','linewidth',2)
for it = 1: numel(Info.tMyrs)
  cf=(it - 1) / (numel(Info.tMyrs)-1);
  plot(Info.tMyrs(it),Info.zSOL(it)/1000,'color',[0,0,cf])
  plot(Info.tMyrs(it),Info.zLAB(it)/1000,'color',[0,0,cf])
end
box on
xlabel('model time [Myr]')
ylabel('zSOL, zLAB [km]')
