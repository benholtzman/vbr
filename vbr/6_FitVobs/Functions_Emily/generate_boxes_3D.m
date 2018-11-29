function generate_boxes_3D(Work)
clc
[calc_boxes_yn] = input('Do you want to calculate new boxes, y/[n]?  ','s');
if ~strcmp(calc_boxes_yn,'y'); return; end



% Parameter sweep
%  define parameter sweep here. var1name must match EXACTLY a field in
%  settings structure. var1 must be defined, var2 lines can be
%  commented/deleted if desired.
settings.Box.Tpot_range = 1150:20:1550;
settings.Box.Tpot_units =' C';

settings.Box.zPlate_range = 60:10:250;
settings.Box.zPlate_units =' km';

settings.Box.phi_range = (0.0:0.002:0.03).*1e2; % melt fraction
settings.Box.phi_units = ' %';

settings.Box.gs_range = (10.^(0:0.5:4)).*1e6; % grain size
settings.Box.gs_units = ' m';


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

drive_plate(Work, settings)

end

function drive_plate(Work, settings_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script builds the initial settings to call TwoPhase. Commonly
% changed parameters are set here, less frequently changed parameters are
% in initialize_settings. Anything set here can overwrite anything set in
% settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp(' ');disp('--------------------------------------------------')
   addpath([Work.cwd '/00_init_conds']);
   addpath([Work.cwd '/02_box_functions']);
   addpath([Work.cwd './01_functions']);
   cur_dir = pwd; cd(Work.cwd); 

%%% --------------------------------------------------------------------- %%
%%%                       Problem Setup
%%% --------------------------------------------------------------------- %%

% settings common to all runs
%  Save settings
   Work.savedir0=[Work.hmdir '/Boxes/']; % the closet to store the box in
   Work.saveindi = 'no'; % save all output? yes or no. [bh: needs more explanation]

%  Load Default Settings
   [settings]=init_settings;

%  Overwrite any of the default settings as desired
    overwrite_settings = fieldnames(settings_in);
    for is = 1:length(overwrite_settings)
        if isstruct(settings_in.(overwrite_settings{is}))
           ovwr_set2 = fieldnames(settings_in.(overwrite_settings{is}));
           for iss = 1:length(ovwr_set2)
               settings.(overwrite_settings{is}).(ovwr_set2{iss}) ...
                   = settings_in.(overwrite_settings{is}).(ovwr_set2{iss});
           end
        else
            settings.(overwrite_settings{is}) = settings_in.(overwrite_settings{is});
        end
    end


% settings for parameter variations
%  Box settings

%
%   Specify data reduction method for Box storage
    settings.Box.DownSampleMeth='interp';
    settings.Box.DownSampleFactor=1;

%%% --------------------------------------------------------------------- %%
%%%                        Data Storage
%%% --------------------------------------------------------------------- %%
    Work.savedir = [Work.savedir0 '/' Work.Box_base_name];
    if exist(Work.savedir,'dir')~=7
        disp('Save directory does not exist, creating it...')
        mkdir(Work.savedir);
    end
    Work.cwd = pwd; cd(Work.savedir); Work.savedirfull = pwd; cd(Work.cwd);
    disp(' ')
    disp(['Saving data in directory ' Work.savedirfull])

    Work.mfile_name = [mfilename('fullpath') '.m'];
    
%%% --------------------------------------------------------------------- %%
%%% --------------  thermodynamic state forward model ------------------- %%
%%% --------------------------------------------------------------------- %%

   [Box,settings] = BuildBox(settings); %% Build The Box (do this once):
   Work.nBox = settings.n_zPlate * settings.n_Tpot * settings.n_phi * ...
       settings.n_gs;

   for iBox = 1:Work.nBox

      disp(' ');disp('--------------------------------------------------')
      disp(['Starting run ' num2str(iBox) ' of ' num2str(Work.nBox)])

%%% ---------------------------- %
%%%    set the Box parameters    %
%%% ---------------------------- %

      settings.zPlate = Box(iBox).info.zPlate_val;
      settings.Tpot   = Box(iBox).info.Tpot_val;
      settings.phi0   = 0.0; % calculate thermal profile w/o melt initially
      settings.phi2   = Box(iBox).info.phi_val;
      settings.grain0 = Box(iBox).info.gs;
      
      frpintf(['\n\tPlate Thickness: %g %s\n\tPotential temperature: '...
          '%g %s\n\tMelt Fraction: %g %s\n\tGrain size: %g micro%s'],...
          settings.zPlate, settings.Box.zPlate_units,...
          settings.Tpot, settings.Box.Tpot_units, ...
          Box(iBox).info.phi_val, settings.Box.phi_units, ...
          settings.grain0.*1e-6, settings.Box.gs_units)

%%% ------------------------------ %%
%%%    Initialize Thermal Solve    %%
%%% ------------------------------ %%

%      build mesh and initial conditions
       settings.Zinfo.zmax=settings.zPlate;
       settings.Zinfo.dz0 = settings.dz0;
       settings.Zinfo = init_mesh(settings.Zinfo); % build the mesh!
       [Info] = init_values(settings); % calculate initial values

%      set the boundary conditions
       [Info.BCs]=init_BCs(struct(),'T','zmin','dirichlet',0);
       [Info.BCs]=init_BCs(Info.BCs,'T','zmax','dirichlet',Info.init.T(end));

%%% ----------------------------- %%
%%%     Forward model solution    %%
%%% ----------------------------- %%

%     Lithosphere temperature evolution
       [Vars,Info]=Thermal_Evolution(Info,settings);

%     Add on adiabatic asthenosphere
       [Vars,Info]=postproc_append_astheno(Vars,Info,settings);

%%% ------------------- %%
%%%     Data Storage    %%
%%% ------------------- %%

%      Put the run in the Box.
       [Box] = Put_in_Box(Box,Vars,Info,settings,iBox);

%      Put the box in the closet (could do this just once at the end)
       Work.savename = [Work.savedir '/Box_' Work.Box_base_name];
       save(Work.savename,'Box')

%      save individual run
       if strcmp(Work.saveindi,'yes')
         if exist([Work.savedir '/individual_runs'],'dir')~=7
             mkdir([Work.savedir '/individual_runs']);
         end
         Work.indifile=[Work.Box_base_name '_' num2str(iBox)];
         Work.savename = [Work.savedir '/individual_runs/' Work.indifile];
         save(Work.savename,'Vars','settings','Info')
         disp(['Run name: ' Work.savedirfull '/individual_runs/' Work.indifile])
       end

%      update the user
       disp(['Completed Run ' num2str(iBox) ' of ' num2str(Work.nBox)])
       disp(['Box name: ' Work.savedirfull '/Box_' Work.Box_base_name])
   end

   Work.savename = [Work.savedir '/Box_' Work.Box_base_name];
   save(Work.savename,'Box')

   clear iBox Vars
disp(' ');disp('--------------------------------------------------');disp(' ')

cd(cur_dir);
end


function [Box,settings] = BuildBox(settings)


% vars = {'zPlate','Tpot','phi','gs'};

n_zPlate = numel(settings.Box.zPlate_range);
n_Tpot   = numel(settings.Box.Tpot_range);
n_phi    = numel(settings.Box.phi_range);
n_gs     = numel(settings.Box.gs_range);

Box = struct();
Box.settings = settings.Box;

for i_zPlate = 1:n_zPlate
    for i_Tpot = 1:n_Tpot
        for i_phi = 1:n_phi
            for i_gs = 1:n_gs
                
                Box(i_zPlate,i_Tpot, i_phi, i_gs).info.zPlate_val = ...
                    settings.Box.zPlate_range(i_zPlate);
                Box(i_zPlate,i_Tpot, i_phi, i_gs).info.Tpot_val = ...
                    settings.Box.Tpot_range(i_Tpot);
                Box(i_zPlate,i_Tpot, i_phi, i_gs).info.phi_val = ...
                    settings.Box.phi_range(i_phi);
                Box(i_zPlate,i_Tpot, i_phi, i_gs).info.gs_val = ...
                    settings.Box.gs_range(i_gs);
                
            end
        end
    end
end

settings.n_zPlate = n_zPlate;
settings.n_Tpot   = n_Tpot;
settings.n_phi    = n_phi;
settings.n_gs     = n_gs;



end
