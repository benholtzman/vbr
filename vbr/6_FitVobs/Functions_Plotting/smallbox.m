%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smallbox.m
% function smallbox(Box_name_IN)
% CJH, 02/2015. 
% Downsamples box frames to single time step. I really just made this so 
% that I can transfer files between orestes and my laptop efficiently for 
% making plots... I like to run it using the bash script smallboxbash. 
%
% It first saves the Box as Smallbox but then overwrites the Frames(:) to 
% a single step Frames(i), where i is chosen based on the closest frame to
% the age_myr variable. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function smallbox(Box_name_IN)
tstart = cputime; 
RunningInBackground = 'yes'; % if yes, forces matlab to exit at completion,  
                             % meant for running matlab in background via 
                             % command line. Use with bash script,
                             % RunMatBG, et to 'no' otherwise.
% chose the frame to save
  age_myr = 80; % the frame time to plot [Myr]

% Path and file info
  Box_dir = '../../../0_BOXES/';
  Box_name_OUT = ['Small' Box_name_IN];
  Box_name_IN = [Box_dir Box_name_IN] 
  Box_name_OUT = [Box_dir Box_name_OUT]

% Load the big box
  load([Box_name_IN]) ; 
  SmallBox = Box; % copy the full box into the small box, then shrink the box...
  Box
  sz = size(Box)
  ni_V1 = sz(1) ; 
  nj_V2 = sz(2) ; 

% LOOP OVER BOX size... 
for i_V1 = 1:ni_V1 % 1:ni_V1
 for j_V2 = 1:nj_V2 % 1:nj_V2
    display('-------------------------------------------------------- ')
    display(['Box number ' num2str(j_V2+(i_V1-1)*nj_V2) ' of ' num2str(ni_V1*nj_V2)])
        
    % pick out the frame
     t_Myr = SmallBox(i_V1,j_V2).Movie.info.timesteps_myrs;
     [val, frame_vec]=min(abs(t_Myr-age_myr));
     i = frame_vec; 

    % overwrite the frame vars with a single time step
    SmallBox(i_V1,j_V2).Movie.info.timesteps=SmallBox(i_V1,j_V2).Movie.info.timesteps(i); 
    SmallBox(i_V1,j_V2).Movie.info.timesteps_myrs=SmallBox(i_V1,j_V2).Movie.info.timesteps_myrs(i); 
    SmallBox(i_V1,j_V2).Movie.info.i_LAB_vec=SmallBox(i_V1,j_V2).Movie.info.i_LAB_vec(i); 
    SmallBox(i_V1,j_V2).Movie.info.q_surf_t=SmallBox(i_V1,j_V2).Movie.info.q_surf_t(i); 
    SmallBox(i_V1,j_V2).Movie.info.q_lab_t=SmallBox(i_V1,j_V2).Movie.info.q_lab_t(i); 
    SmallBox(i_V1,j_V2).Movie.Frames=SmallBox(i_V1,j_V2).Movie.Frames(i); 

 end
end
Box = SmallBox; 
save([Box_name_OUT],'Box')

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
 display(['     Matlab shutting down... Farewell!']);
 quit force
end
end
