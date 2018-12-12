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

drive_svfm_plate(Work, settings)

end

function drive_svfm_plate(Work, settings_in)
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

%  Overwrite any of the deafault settings as desired
    ovwr_set = fieldnames(settings_in);
    for is = 1:length(ovwr_set)
        if isstruct(settings_in.(ovwr_set{is}))
           ovwr_set2 = fieldnames(settings_in.(ovwr_set{is}));
           for iss = 1:length(ovwr_set2)
               settings.(ovwr_set{is}).(ovwr_set2{iss}) ...
                   = settings_in.(ovwr_set{is}).(ovwr_set2{iss});
           end
        else
            settings.(ovwr_set{is}) = settings_in.(ovwr_set{is});
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
    
    %Work.cp_src=system(['04_scripts/sh_source_copy.sh ' ...
    %                                    Work.savedir ' ' Work.mfile_name]);

%%% --------------------------------------------------------------------- %%
%%% --------------  thermodynamic state forward model ------------------- %%
%%% --------------------------------------------------------------------- %%

   [Box,settings] = BuildBox(settings); %% Build The Box (do this once):
   Work.nBox=settings.Box.nvar1 * settings.Box.nvar2;

   for iBox = 1:Work.nBox

      disp(' ');disp('--------------------------------------------------')
      disp(['Starting run ' num2str(iBox) ' of ' num2str(Work.nBox)])

%%% ---------------------------- %
%%%    set the Box parameters    %
%%% ---------------------------- %

      settings.(settings.Box.var1name) = Box(iBox).info.var1val;
      disp([settings.Box.var1name '=' num2str(settings.(settings.Box.var1name))...
           settings.Box.var1units]);
      if isfield(settings.Box,'var2name')
        settings.(settings.Box.var2name) = Box(iBox).info.var2val;
        disp([settings.Box.var2name '=' num2str(settings.(settings.Box.var2name))...
           settings.Box.var2units]);
      end

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
