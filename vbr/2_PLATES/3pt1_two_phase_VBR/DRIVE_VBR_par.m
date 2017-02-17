%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel VBR driver
% batch-parallel VBR processing for each Box index. 
% 
% Only saves a subset of VBR output -- see VBRsave structure in parfor
% loop. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
tstart = cputime; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parallel and background settings %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  procs_to_use = 3; % number of processors to use (if no processor pool  
                    % exists already)
  RunningMore = 'no'; % are you going to run more in this session? If yes, 
                       % the matlab parallel pool is not closed at the end 
                       % to reduce run time on subsequent runs. 
                     
  RunningInBackground = 'yes'; % if yes, forces matlab to exit at completion,  
                              % meant for running matlab in background via 
                              % command line. Use with bash script,
                              % RunMatBG, set to 'no' otherwise.

%%%%%%%%%%%%%%       
% file setup %
%%%%%%%%%%%%%%

  Box_base_name='2016-04-21-g0_vs_muf_con'; % the box to load
    
  cwd=pwd;cd ~;hmdir=pwd;cd(cwd)
  Box_dir =[hmdir '/Dropbox/Research/0_Boxes/']; % box directory
  Box_dir = [Box_dir Box_base_name '/'];
  Box_name_IN = ['Box_' Box_base_name];  
  Box_name_IN = [Box_dir Box_name_IN]; % full path to load
  Box_name_OUT = [Box_name_IN , '_VBR']; % full path to write
  load(Box_name_IN); % actually load the box
  
%%%%%%%%%%%%%
% VBR setup %
%%%%%%%%%%%%%
  freq_vec = logspace(-2.2,-1.3,4); % frequencies to calculate at [Hz]
  MELT = 1; % 0 or 1 to zero or keep melt fraction, phi
  
  addpath('./02_box_functions', './03_plotting')
  VBR_version = 'VBR_v0p94' ;
  addpath(['../../4_VBR/' VBR_version ],...
          ['../../4_VBR/' VBR_version '/functions'],...
          ['../../4_VBR/' VBR_version '/params'])
      
  Params_QMV_scaling; % initialize VBR
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  initialize parallel pool %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   if isempty(gcp('nocreate'))
       c = parcluster('local'); % build the 'local' cluster object
       nw = c.NumWorkers;        % get the number of workers
       if procs_to_use > nw; 
           procs_to_use = nw; 
           disp(['specificed processors exceeds max (' num2str(nw) '), using max']);
       end
       pool = parpool(procs_to_use);
   end   
   if exist('parallel_monitor','dir')==0; mkdir('parallel_monitor'); end   
   
%%%%%%%%%%%%%%%%%%%%%
% set up the VBRBox %
%%%%%%%%%%%%%%%%%%%%%
  VBR0 = VBR; 
  VBRBox = struct('VBRsetup',[],'Movie',[]);
  [n1,n2]=size(Box);
  for i1=1:n1
      for i2=1:n2
          VBRBox(i1,i2).VBRsetup=VBR0;
          VBRBox(i1,i2).VBRsetup.f=freq_vec;
      end
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run VBR calculator for each Box index %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nBox = numel(Box);
  parfor iBox = 1:nBox
        
%       Extract current box	values	
		Movie = Box(iBox).Movie; 
		nt = numel(Movie.info.timesteps) ; 
			
%       Initialize current node's VBR and Movie structure        
        VBR = VBR0; 
		VBR.ISV.f = freq_vec ;        
        MovieNode = struct(); 
        
%       Loop over timesteps, calling VBR calculator        
        for it = 1:nt            
%           pull out the current frame            
            VBR.ISV.P_GPa = (Movie.Frames(it).P)./1e9 ;
            VBR.ISV.T_K = Movie.Frames(it).T +273;
            VBR.ISV.rho = Movie.Frames(it).rho ;
            VBR.SSV.phi = MELT * Movie.Frames(it).phi ; 
            VBR.SSV.dg_um = Movie.Frames(it).dg_um ;
            VBR.ISV.sig_MPa = Movie.Frames(it).sig_MPa ;
            VBR.Gu_0 = Movie.Frames(it).Gu_0_GPa*1e9 ; 
            VBR.CSV.Ch2o = Movie.Frames(it).Cs; 
            
%           It's VBR time!
            [VBR] = VBR_spine(VBR) ;        
            
%           Store the VBR structure in the current frame -- only store the
%           variables that will be saved. 
            VBRsave=struct(); 
            VBRsave.Vave_eBurgers=VBR.eBurgers.Vave; 
            VBRsave.Vave_AndradePsP=VBR.AndradePsP.Vave; 

%           Store it in the current node's Movie structure  
            MovieNode.Frames(it).VBR = VBRsave ;               
        end
        
        VBRBox(iBox).Movie=MovieNode; % put it in the Box

%      check parallel progress (counts number of files)
       ThisFile=fopen(['parallel_monitor/Run' num2str(iBox) '.txt'],'w');
       fwrite(ThisFile,['Completed Run ' num2str(iBox)]);
       fclose(ThisFile);
       D = dir('parallel_monitor');
       Num = length(D(not([D.isdir])));
       disp(['   TOTAL PROGRESS ---> ' num2str(Num) ' of ' num2str(nBox)])     
  end
  
% shutdown parallel pool  
  if strcmp(RunningMore,'no') || strcmp(RunningInBackground,'yes');      
      poolobj = gcp('nocreate');
      delete(poolobj);
  end

% save it   
  save(Box_name_OUT,'VBRBox')

% empty and remove the parallel monitoring folder   
  delete parallel_monitor/*
  rmdir('parallel_monitor')

% endgame  
  tend = cputime; 
  display(' ')
  display('    _____________________________')
  display('   |                             |')
  display('   |   Computations complete     |')
  display('   |_____________________________|')
  display(' ')
  display(['     Total CPU Time: ' num2str((tend-tstart)/60) ' min']);

% Quit matlab if running in background from command line
  if strcmp(RunningInBackground,'yes')==1    
      disp(' ')
      for it = 10:-1:1
          display(['     Matlab shutting down in ' num2str(it)]);
          pause(1)
      end
      display('     Matlab shutting down NOW!!');
      quit force
  end