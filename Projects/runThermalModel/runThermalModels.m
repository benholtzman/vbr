%% ===================================================================== %%
%%                     runThermalModels.m
%% ===================================================================== %%
%  example of running multiple forward models using the box framework
%% ===================================================================== %%
clear

% put VBR in the path
path_to_top_level_vbr='../../';
addpath(path_to_top_level_vbr)
addpath('./extraFun')
vbr_init

% Parameter sweep settings
%  define parameter sweep here. var1name must match EXACTLY a field in
%  settings structure. var1 must be defined, var2 lines can be
%  commented/deleted if desired.
settings.Box.var1range = [1200,1400];
settings.Box.var1name = 'Tpot';
settings.Box.var1units =' C';

settings.Box.var2range = [0,150,300];
settings.Box.var2name = 'Cs0_H2O';
settings.Box.var2units =' PPM';

% overwrite/set any of the parameters
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

%  Specify data reduction method for Box storage
settings.Box.DownSampleMeth='interp';
settings.Box.DownSampleFactor=1;
[Box,settings] = BuildBox(settings); %% Build The Box (do this once):

% Load Default Settings for the forward model runs
[settings]=init_settings(settings);

% loop over Parameter sweep settings settings, run forward model for each permuatation
Work.nBox=settings.Box.nvar1 * settings.Box.nvar2;
for iBox = 1:Work.nBox

  disp(' ');disp('--------------------------------------------------')
  disp(['Starting run ' num2str(iBox) ' of ' num2str(Work.nBox)])

  % update the current settings for this permutation
  settings.(settings.Box.var1name) = Box(iBox).info.var1val;
  disp([settings.Box.var1name '=' num2str(settings.(settings.Box.var1name))...
       settings.Box.var1units]);
  if isfield(settings.Box,'var2name')
    settings.(settings.Box.var2name) = Box(iBox).info.var2val;
    disp([settings.Box.var2name '=' num2str(settings.(settings.Box.var2name))...
       settings.Box.var2units]);
  end

  % build the initial conditions for this permutation
  settings.Zinfo.zmax=settings.zPlate;
  settings.Zinfo.dz0 = settings.dz0;
  settings.Zinfo = init_mesh(settings.Zinfo); % build the mesh!
  [Info] = init_values(settings); % calculate initial values

  % initalize  the boundary conditions
  [Info.BCs]=init_BCs(struct(),'T','zmin','dirichlet',0);
  [Info.BCs]=init_BCs(Info.BCs,'T','zmax','dirichlet',Info.init.T(end));

  % Lithosphere temperature evolution (the actual calculation)
  [Vars,Info]=Thermal_Evolution(Info,settings);

  % Add on adiabatic asthenosphere
  [Vars,Info]=postproc_append_astheno(Vars,Info,settings);

  % Put the run in the Box.
  [Box] = Put_in_Box(Box,Vars,Info,settings,iBox);
end

disp(' ');disp('--------------------------------------------------');disp(' ')

disp('Computations complete, plotting final temperature profile from each run')

figure()

subplot(2,2,1)
iBox=6;
hold on
plot(Box(iBox).Frames(end).Tsol,Box(iBox).run_info.Z_km,'--k')
for it = 1:2: numel(Box(iBox).run_info.timesteps_myrs)
  cf=(it - 1) / (numel(Box(iBox).run_info.timesteps_myrs)-1);
  plot(Box(iBox).Frames(it).T,Box(iBox).run_info.Z_km,'color',[0,0,cf])
end
box on
set(gca,'ydir','reverse')
xlabel('geotherms and solidii [C]')
ylabel('depth [km]')

subplot(2,2,2)
zIso=extractIsotherm(Box(iBox),1200);

plot(Box(iBox).run_info.timesteps_myrs,Box(iBox).run_info.zSOL/1000,'k','linewidth',2)
hold on
plot(Box(iBox).run_info.timesteps_myrs,Box(iBox).run_info.zLAB/1000,'--k','linewidth',2)
plot(Box(iBox).run_info.timesteps_myrs,zIso,'--r','linewidth',2)

for it = 1: numel(Box(iBox).run_info.timesteps_myrs)
  cf=(it - 1) / (numel(Box(iBox).run_info.timesteps_myrs)-1);
  plot(Box(iBox).run_info.timesteps_myrs(it),Box(iBox).run_info.zSOL(it)/1000,'color',[0,0,cf])
  plot(Box(iBox).run_info.timesteps_myrs(it),Box(iBox).run_info.zLAB(it)/1000,'color',[0,0,cf])
  plot(Box(iBox).run_info.timesteps_myrs(it),zIso(it),'color',[0,0,cf])
end
box on
xlabel('model time [Myr]')
ylabel('zSOL, zLAB [km]')


subplot(2,2,3)
hold on
for iv2 = 1: settings.Box.nvar2
  cf=(iv2 - 1) / (settings.Box.nvar2 -1);
  plot(Box(1,iv2).Frames(end).Tsol,Box(1,iv2).run_info.Z_km,'color',[0,cf,0],'linestyle','--')
end
for iv1 = 1:settings.Box.nvar1
  cf=(iv1 - 1) / (settings.Box.nvar1 -1);
  plot(Box(iv1,1).Frames(end).T,Box(iv1,1).run_info.Z_km,'color',[cf,0,0])
end
box on
set(gca,'ydir','reverse')
xlabel('final geotherms and solidii [C]')
ylabel('depth [km]')

subplot(2,2,4)
hold on
for iBox=1:Work.nBox
  cf=(iBox - 1) / (Work.nBox  -1);
  plot(Box(iBox).run_info.timesteps_myrs,Box(iBox).run_info.zSOL/1000,'color',[cf,0,0],'marker','.')
end
box on
ylim([0,settings.Zinfo.asthenosphere_max_depth+10])
xlabel('model time [Myr]')
ylabel('solidus-geotherm intersection [km]')
