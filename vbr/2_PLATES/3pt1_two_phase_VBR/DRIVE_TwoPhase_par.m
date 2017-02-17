%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script builds the initial settings to call TwoPhase. Commonly
% changed parameters are set here, less frequently changed parameters are
% in initialize_settings.  
%
% a slightly parallelized version...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp(' ');disp('--------------------------------------------------')
   clear; close all
   addpath ./01_init_conds_and_material   
   addpath ./02_box_functions
   addpath ./02_functions
   
%% --------------------------------------------------------------------- %%
%%                       Problem Setup      
%% --------------------------------------------------------------------- %%

   disp('Setting up the calculation')

%  Parallel settings   
   Work.procs_to_use =7;  % processors to use
   Work.par_run_multiple='no'; % 'no' will close parallel pool, 'yes' will keep it open
   
%  Save settings
   Work.cwd=pwd; cd ~; Work.hmdir=pwd; cd(Work.cwd);
   Work.savedir0=[Work.hmdir '/Dropbox/Research/0_Boxes']; % the closet to store the box in 
   Work.savebase = '2016-07-22-re_Vbg_zLAB_varXf'; % boxes end up named ['Box_' savebase]            
   Work.saveindi = 'yes'; % save all output? yes or no. 
   
%  Load Default Settings  
   [settings]=init_settings;
   
%  Overwrite any of the deafault settings as desired
%    material properties     
     settings.Z_moho_km = 10; % Moho depth [km] 
     
%    Mesh
     settings.dz0 = 0.5; % target node spacing in [km]
     settings.Zinfo.zmin =0; % min depth for model domain [km]
     settings.Zinfo.zmax = 200; % max z depth [km]
     
%    Computational settings 
%    time
     settings.nt= 50000; % max time steps 
     settings.outk =250; % output frequency          
%      settings.nt= 4; % max time steps 
%      settings.outk =2; % output frequency  
     settings.t_max_Myrs=1.05; % max time to calculate [Myr]     
     settings.z_thin=settings.Z_moho_km; % quit if lith is thinned to here [km]
     
%    for melt fraction calc          
     settings.phimax = 0.15; % maximum phi     
     settings.dt_max = 10;  % max step for advection [Myrs] if advective velo is 0
       
%  Flags
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
%    .H2O       'constant' or 'variable' 
%    .XtalFactor 'constant' or 'variable'
%    .problemkill  'stop_on_no_melt' or 'keepgoing_on_no_melt'
%    .phikill    'stop_on_max_phi' or 'keepgoing'
     settings.Flags.PropType='PT_dep';         
     settings.Flags.problem='Two_Phase_Sd';
     settings.Flags.LABdef='phi';
     settings.Flags.VbzFlag = 'variable';
     settings.Flags.H2O = 'constant';
     settings.Flags.progress_plot = 'no';
     settings.Flags.parg = 'no';
     settings.Flags.XtalFactor = 'variable';
     settings.Flags.problemkill='stop_on_no_melt';
     settings.Flags.phikill='keepgoing';
     
%  Settings for Initial conditions 
     settings.Tpot = 1325; % 
     settings.Tpot_excess = 1475; % 
     settings.T_init='continental'; % 'continental' 'oceanic' or 'adiabatic'
     settings.T_init_Zlab=80; % [km] 
     settings.age0 = 100; % [Myrs] only used if T_init is oceanic. 
     settings.Vbg = 10; % [cm/yr]  
     settings.DBL_m = 10e3; % [m] mechancial boundary layer for upwelling V
     settings.phi_init = 0.005;
     settings.grain0 = 0.01; % grain size [m]
     settings.mufo = 1; % fluid viscosity [Pa s]
     settings.Gu_olv = 66; % ref shear mod [GPa]
     settings.Cs0 = 0 * 1e-4; % H2O concentration [wt % -- needs to be wt % for solidus]
     settings.Cs0_CO2 = 0 *1e-4; % CO2 concentration [wt %]
     settings.y_sl = 1; 
     settings.Sd_coefficient = .001;      
     settings.DikingBL_km = 5; % [km] diking process zone boundary layer     
     settings.dTdz_ad = .45*1e-3; % [C/m] adiabatic gradient
     settings.XtalFactor = 1; % fraction of infiltrated melt that freezes
     settings.XtalFactor_z0 = settings.Z_moho_km; % controls offset [km]
     settings.XtalFactor_z1 = 15; % controls gradient [km]
     
%  Box settings
%  define parameter sweep here. var1name must match EXACTLY a field in 
%  settings structure. var1 must be defined (for now - change this?),
%  var2 lines can be commented/deleted if desired.

     settings.Box.var2range = 50:5:110;
     settings.Box.var2name = 'T_init_Zlab';
     settings.Box.var2units =' km';
   
%     settings.Box.var2range = 1350:10:1550; 
%     settings.Box.var2name = 'Tpot_excess';
%     settings.Box.var2units =' ^oC';    

      settings.Box.var1range = 2.5:2.5:25;
      settings.Box.var1name = 'Vbg';
      settings.Box.var1units =' cm/yr';
   

%    Specify data reduction method for Box storage     
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
       Work.currentD = pwd; cd(Work.savedir); Work.savedirfull = pwd; cd(Work.currentD); 
       disp(' ')
       disp(['Saving data in directory ' Work.savedirfull])
       
       Work.mfile_name = [mfilename('fullpath') '.m'];
       Work.sys=system(['04_scripts/sh_source_copy.sh ' Work.savedir ' ' Work.mfile_name]);
     
%% --------------------------------------------------------------------- %%
%% -----------------  box and paralell setup --------------------------- %%
%% --------------------------------------------------------------------- %%
      
   [Box,settings] = BuildBox(settings); %% Build The Box (do this once):  
   Work.nvar1=settings.Box.nvar1; Work.nvar2=settings.Box.nvar2;  
   Work.nBox=Work.nvar1*Work.nvar2;
%  vectorize var1,var2 for parallel efficiency
   InitValuesVar1=zeros(Work.nBox,1);
   InitValuesVar2=zeros(Work.nBox,1);
   
   Work.ivec = 1; 
   for ivar1=1:Work.nvar1
       for ivar2=1:Work.nvar2
           InitValuesVar1(Work.ivec)=settings.Box.var1range(ivar1);
           if isfield(settings.Box,'var2name')
               InitValuesVar2(Work.ivec)=settings.Box.var2range(ivar2);
           end
           Work.ivec = Work.ivec+1;  
           BoxVec(Work.ivec)=struct();
       end
   end   
   
%  initialize parallel pool   
   if isempty(gcp('nocreate'))
       Work.c = parcluster('local'); % build the 'local' cluster object
       Work.nw = Work.c.NumWorkers;        % get the number of workers
       if Work.procs_to_use > Work.nw; 
           Work.procs_to_use = Work.nw; 
           disp(['specificed processors exceeds max (' num2str(Work.nw) '), using max']);
       end
       Work.pool = parpool(Work.procs_to_use);
   end
   
   settings0=settings; % save it so each proc can modify it
   rng('shuffle') 
   Work.parallel_monitor_dir=['parallel_monitor_' num2str(round(rand(1)*1000))];
   mkdir(Work.parallel_monitor_dir)
%% --------------------------------------------------------------------- %%
%% --------------------  thermodynamic state --------------------------- %%
%% --------------------------------------------------------------------- %%   
   parfor ivec=1:Work.nBox
        
      settings=settings0;
      disp(' ');disp('--------------------------------------------------')        
      disp(['Starting run ' num2str(ivec) ' of ' num2str(Work.nBox)])  

%     set the parameters        
      settings.(settings.Box.var1name) = InitValuesVar1(ivec);
      disp([settings.Box.var1name '=' num2str(settings.(settings.Box.var1name))...
           settings.Box.var1units]);
      if isfield(settings.Box,'var2name');
        settings.(settings.Box.var2name) = InitValuesVar2(ivec);
        disp([settings.Box.var2name '=' num2str(settings.(settings.Box.var2name))...
           settings.Box.var2units]);
      end    
      
%    calculate initial conditions 
     settings.Zinfo.dz0 =  settings.dz0;
     settings.Zinfo = init_mesh(settings.Zinfo); % build the mesh!
     [Info] = init_values(settings); % for T, phi, Vbg
%    set the boundary conditions      
     Info.BCs = init_BCs(Info,settings);
      
%    solve it      
     [Vars,Info_out]=TwoPhase_SSC(Info,settings);
     
%    monitor progress     
     disp(['Completed Run ' num2str(ivec) ' of ' num2str(Work.nBox)])
     ThisFile=fopen([Work.parallel_monitor_dir '/Run' num2str(ivec) '.txt'],'w');
     fwrite(ThisFile,['Completed Run ' num2str(ivec)]);
     fclose(ThisFile);
     
     D = dir(Work.parallel_monitor_dir);
     Num = length(D(not([D.isdir])));
     
     disp('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\')
     disp(['   TOTAL PROGRESS ---> ' num2str(Num) ' of ' num2str(Work.nBox)])
     disp('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\')
%    put it in the 1-D box      
     BoxVec(ivec).Info=Info_out;
     BoxVec(ivec).settings=settings;
     BoxVec(ivec).Vars=Vars;
    
   end
   
%  empty and remove the parallel monitoring folder   
   delete([Work.parallel_monitor_dir '/*'])
   rmdir(Work.parallel_monitor_dir)
   
%  move the 1D box back to 2D  
   ivec = 1;
   for ivar1=1:Work.nvar1
       for ivar2=1:Work.nvar2
           Box0(ivar1,ivar2)=BoxVec(ivec);
           ivec = ivec+1;
       end
   end
   
%  close up parallel
   Work.poolobj = gcp('nocreate');
   if strcmp(Work.par_run_multiple,'yes')==0
       delete(Work.poolobj);
   end
   
%  store the temp box in a real box   
   [Box] = par_Put_in_Box(Box0,Box);
   
%  Put the box in the closet         
   Work.savename = [Work.savedir '/Box_' Work.savebase];
   save(Work.savename,'Box')  
   disp(['Parallel runs complete, Box name: ' Work.savedirfull '/Box_' Work.savebase])
       
   clear ivec D Num ivar1 ivar2 settings0 InitValuesVar1 InitValuesVar2 ...
         BoxVec Box0 settings
disp(' ');disp('--------------------------------------------------');disp(' ')
