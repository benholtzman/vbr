clear
addpath('../../6_FitVobs/Functions_Plotting/')
addpath('../../6_FitVobs/')
addpath ./02_box_functions 
close all;

% load Box data
  Box_dir = '../../../../0_BOXES/';
  Box_name0 ='Box_halfspace_VBR';
  Box_name = [Box_dir Box_name0 '.mat']
  load(Box_name)
  
% WHICH TIME STEP to look at? 
  age_myr = 80; % age to plot [myr]
  

% get VarInfo
  VarInfo.Var1_name=Box(1,1).info.var1name;
  VarInfo.Var1_units=Box(1,1).info.var1units; 
  VarInfo.Var1_range = Box(1,1).info.var1range; 
%   if isfield(Box(1,1).info,'var2range')
    VarInfo.Var2_name=Box(1,1).info.var2name;
    VarInfo.Var2_units=Box(1,1).info.var2units;
    VarInfo.Var2_range = Box(1,1).info.var2range;
%   end
  VarInfo.Var1_n=numel(Box(:,1));
  VarInfo.Var2_n=numel(Box(1,:)); 
  
%   for i_Var1 = 1:VarInfo.Var1_n
%       VarInfo.Var1_range(i_Var1) = Box(i_Var1,1).info.T_pot_C;      
%   end
%   if VarInfo.Var2_n > 1
%   for j_Var2 = 1:VarInfo.Var2_n
%       VarInfo.Var2_range(j_Var2) = Box(1,j_Var2).info.ZLAB/1e3;      
%   end
%   end 
  for i = 1:VarInfo.Var1_n
      for j = 1:VarInfo.Var2_n
      [val, tsnap(i,j)]=min(abs(Box(i,j).Movie.info.timesteps_myrs-age_myr));
      end
  end
% set depth range for all plots
  ylimits = [0 350]; % NOMELT Vs only resolved to 325 km. Below that, the 
                     % observations are smoothed to a reference model.
  xlimits = [min(VarInfo.Var2_range) max(VarInfo.Var2_range)];                     
  maxZ_km = max(Box(1,1).Movie.info.Z_km) ;                    
% load, modify NOMELT Vs profile
  velfile = '../../6_FitVobs/velocity_models/Nomlet_Vs';
  nomeltobs = load(velfile)
  obsdepth = nomeltobs.Z';
  obsVs    = nomeltobs.Vs';
  obsdepth = obsdepth(3:end); % remove water layer
  obsVs = obsVs(3:end);
  obsVs = obsVs(obsdepth < max(maxZ_km)); % remove points deeper than model domain
  obsdepth = obsdepth(obsdepth < max(maxZ_km));
% choose which VBR to plot
  Vel_ExptSet = 3
   % 1 = Andrade.Va 
   % 2 = Andrade.Va_comp
   % 3 = AndradePsP.Va
   % 4 = AndradePsP.Va_comp
   % 5 = eBurgers.V

% set depth weighting (weighting for all depths initialized to 1, set changes here)  
  depthrange = [0 50 ; 50 100; 100 600 ];
  weighted = [ 0 0 1  ];
  
% calculate and plot the best fitting Vs


[i_Var_b,j_Var_b,msfit,wtz]=Calc_bestFit_ch(Box,Vel_ExptSet,tsnap,depthrange,...
                                   weighted,obsVs,obsdepth);  


[F1]=PLOT_bestFit_ch(Box,VarInfo,tsnap,obsVs,obsdepth,...
                                  ylimits,xlimits,i_Var_b,j_Var_b,msfit,wtz);
                               
% plot the full range of tests
F2=PLOT_wholeBox_ch(Box,VarInfo,i_Var_b,j_Var_b,obsVs,obsdepth,tsnap,ylimits,'fullbox');

% plot the best fitting profiles  
F3=PLOT_bestProfiles_ch(Box,VarInfo,i_Var_b,j_Var_b,tsnap,obsVs,obsdepth,ylimits);

BN = [Box_dir Box_name0];
saveas(F1,[BN '_Fit1.eps'],'epsc')
saveas(F2,[BN '_Fit2.eps'],'epsc')
saveas(F3,[BN '_Fit3.eps'],'epsc')









