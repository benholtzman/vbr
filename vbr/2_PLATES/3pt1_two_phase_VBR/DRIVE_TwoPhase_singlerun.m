%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script builds the initial settings to call TwoPhase. Commonly
% changed parameters are set here, less frequently changed parameters are
% in initialize_settings.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp(' ');disp('--------------------------------------------------')
   clear; close all   
   addpath ./01_init_conds_and_material ./02_box_functions ./02_functions
  
%% --------------------------------------------------------------------- %%
%%                       Problem Setup      
%% --------------------------------------------------------------------- %%
   
%% settings common to all runs   
%  Save settings
   Work.cwd=pwd; cd ~; Work.hmdir=pwd; cd(Work.cwd);
   Work.savedir0=[Work.hmdir '/Desktop']; % the closet to store the box in 
   Work.savebase = '2016-02-16-TEST'; % boxes end up named ['Box_' savebase]   
   Work.saveindi = 'yes'; % save all output? yes or no. 
   
%  Load Default Settings  
   [settings]=init_settings;
   
%  Overwrite any of the deafault settings as desired
%    Mesh
     settings.dz0=.1; 
     settings.Zinfo.zmin =0; % min depth for model domain [km]
     settings.Zinfo.zmax = 200; % max z depth [km]
     
%    Computational settings 
%    time
     settings.nt= 100; % max time steps 
     settings.outk =5; % output frequency
     settings.t_max_Myrs=0.5; % max time to calculate [Myr]     
     settings.z_thin=settings.Z_moho_km; % quit if lith is thinned to here [km]
     
%    for melt fraction calc          
     settings.phimax = 0.15; % maximum phi     
     settings.dt_max = 10;  % max step for advection [Myrs] if advective velo is 0
       
%  Problem Flags 
%    .problem   'One_Phase_T' = single phase energy conservation only
%               'Two_Phase_Sd' = dike-flux melt infiltration
%    .PropType  specifies dependencies of Kc, rho and Cp.           
%               'con'     Constant rho, Cp and Kc
%               'P_dep'   pressure dependent rho, constant Cp and Kc
%               'T_dep'   temperature dependent rho, Cp and Kc
%               'PT_dep'  temperature and pressure dependent rho, Cp and Kc
%    .LABdef    sets the method for calculating the "LAB." The melt 
%               fraction evolution always uses the phi-LAB, this flags is 
%               for the upwelling velocity.              
%               'visc'    eta/eta_astheno > 10
%               'phi'     first point above solidus
%    .VbzFlag   sets the method for caculating or setting Vbg
%               'constant'   constant value (i.e., V +S = Vo)
%               'variable'   tapers to 0 at LAB, recalculates as LAB moves
%               'var_z_con_t' tapers to 0 at LAB, doesn't recalculate
%    .XtalFactor 'constant' or 'variable'
%    .problemkill  'stop_on_no_melt' or 'keepgoing_on_no_melt'
%    .phikill    'stop_on_max_phi' or 'keepgoing'
%    .H2O       'constant' or 'variable' (KEEP AT CONSTANT)
%    .parg       'yes' or 'no' (KEEP AT NO)
     settings.Flags.PropType='PT_dep';         
     settings.Flags.problem='Two_Phase_Sd';
     settings.Flags.LABdef='phi';
     settings.Flags.VbzFlag = 'var_z_con_t';    
     settings.Flags.progress_plot = 'no';
     settings.Flags.XtalFactor = 'variable';
     settings.Flags.problemkill='stop_on_no_melt';
     settings.Flags.phikill='stop_on_max_phi'; 
     settings.Flags.H2O = 'constant'; % KEEP AT CONSTANT
     settings.Flags.parg = 'no'; % KEEP AT NO
     
%  Settings for Initial conditions 
     settings.Tpot = 1325; % 
     settings.Tpot_excess = 1450; % 
     settings.T_init='continental'; % 'continental' 'oceanic' or 'adiabatic'
     settings.T_init_Zlab=80; % [km] 
     settings.age0 = 100; % [Myrs] only used if T_init is oceanic. 
     settings.Vbg = 15; % [cm/yr]  
     settings.DBL_m = 10e3; % [m] mechancial boundary layer for upwelling V
     settings.phi_init = 0.005;
     settings.grain0 = 0.01; % grain size [m]
     settings.mufo = 1; % fluid viscosity [Pa s]
     settings.Gu_olv = 66; % ref shear mod [GPa]
     settings.Cs0 = 0 * 1e-4; % H2O concentration [wt % -- needs to be wt % for solidus]
     settings.Cs0_CO2 = 0 *1e-4; % CO2 concentration [wt %]
     settings.y_sl = 1; % surface tension energy term (not used?)
     settings.Sd_coefficient = .001; % controls LAB "leakiness", 0 should make it impermeable (watch out!)
     settings.DikingBL_km = 5; % [km] diking process zone boundary layer     
     settings.dTdz_ad = .45*1e-3; % [C/m] adiabatic gradient
     settings.XtalFactor = 1; % fraction of infiltrated melt that freezes
     settings.XtalFactor_z0 = settings.Z_moho_km; % controls offset [km]
     settings.XtalFactor_z1 = 15; % controls gradient [km]
     
%% settings for parameter variations     
%  Box settings
%  define parameter sweep here. var1name must match EXACTLY a field in 
%  settings structure. var1 must be defined, var2 lines can be 
%  commented/deleted if desired. Will overwrite previous var1/var2
%  definitions above. 
% 
    settings.Box.var1range = [80];
    settings.Box.var1name = 'T_init_Zlab';
    settings.Box.var1units =' km';    
     
    settings.Box.var2range = [1525];
    settings.Box.var2name = 'Tpot_excess';
    settings.Box.var2units =' C';
%      
%   Specify data reduction method for Box storage     
    settings.Box.DownSampleMeth='interp';
    settings.Box.DownSampleFactor=1;      
   
%% --------------------------------------------------------------------- %%   
%%                        Data Storage   
%% --------------------------------------------------------------------- %% 
    Work.savedir = [Work.savedir0 '/' Work.savebase];
    if exist(Work.savedir,'dir')~=7
        disp('Save directory does not exist, creating it...')
        mkdir(Work.savedir);
    end
    Work.cwd = pwd; cd(Work.savedir); Work.savedirfull = pwd; cd(Work.cwd);
    disp(' ')
    disp(['Saving data in directory ' Work.savedirfull])
    
    Work.mfile_name = [mfilename('fullpath') '.m'];
    Work.cp_src=system(['04_scripts/sh_source_copy.sh ' ...
                                        Work.savedir ' ' Work.mfile_name]);
       
%% --------------------------------------------------------------------- %%
%% --------------  thermodynamic state forward model ------------------- %%
%% --------------------------------------------------------------------- %%
        
   [Box,settings] = BuildBox(settings); %% Build The Box (do this once):  
   Work.nBox=settings.Box.nvar1 * settings.Box.nvar2; 
   
   for iBox = 1:Work.nBox
       
      disp(' ');disp('--------------------------------------------------')        
      disp(['Starting run ' num2str(iBox) ' of ' num2str(Work.nBox)])  
      
%% ---------------------------- %
%%    set the Box parameters    %
%% ---------------------------- %

      settings.(settings.Box.var1name) = Box(iBox).info.var1val;
      disp([settings.Box.var1name '=' num2str(settings.(settings.Box.var1name))...
           settings.Box.var1units]);
      if isfield(settings.Box,'var2name');
        settings.(settings.Box.var2name) = Box(iBox).info.var2val;
        disp([settings.Box.var2name '=' num2str(settings.(settings.Box.var2name))...
           settings.Box.var2units]);
      end    
      
%% ------------------------------ %%      
%%    Initialize Thermal Solve    %%
%% ------------------------------ %%

%      build mesh and initial conditions
       settings.Zinfo.zmax=settings.Zinfo.zmax; 
       settings.Zinfo.dz0 = settings.dz0; 
       settings.Zinfo = init_mesh(settings.Zinfo); % build the mesh!
       [Info] = init_values(settings); % calculate initial values
       
%      set the boundary conditions             
       Info.BCs = init_BCs(Info,settings);
       
%% ----------------------------- %%
%%     Forward model solution    %%
%% ----------------------------- %%

%%     Lithosphere temperature evolution
       [Vars,Info]=TwoPhase_SSC(Info,settings);
       
%% ------------------- %%
%%     Data Storage    %%
%% ------------------- %%

%      Put the run in the Box.
       [Box] = Put_in_Box(Box,Vars,Info,settings,iBox);
              
%      Put the box in the closet (could do this just once at the end)
       Work.savename = [Work.savedir '/Box_' Work.savebase];
       save(Work.savename,'Box')   
                    
%      save individual run       
       if strcmp(Work.saveindi,'yes')
         if exist([Work.savedir '/individual_runs'],'dir')~=7
             mkdir([Work.savedir '/individual_runs']);
         end
         Work.indifile=[Work.savebase '_' num2str(iBox)];
         Work.savename = [Work.savedir '/individual_runs/' Work.indifile];
         save(Work.savename,'Vars','settings','Info')
         disp(['Run name: ' Work.savedirfull '/individual_runs/' Work.indifile])
       end
       
%      update the user       
       disp(['Completed Run ' num2str(iBox) ' of ' num2str(Work.nBox)])
       disp(['Box name: ' Work.savedirfull '/Box_' Work.savebase])       
   end
   
   Work.savename = [Work.savedir '/Box_' Work.savebase];
   save(Work.savename,'Box')
   
   clear iBox Vars 
disp(' ');disp('--------------------------------------------------');disp(' ')
