function [depth,obsVs] = load_velprofile(velfile,depth_column,vel_column,Z_km_calc)
    
%    vel=load([velfile,'.txt']);
disp velfile

vel=load(velfile, '-ascii');
display(['num vals in Z_km_calc: ' num2str(length(Z_km_calc))])
display(['depth of Z_km_calc: ' num2str(Z_km_calc(end))])

Z_km_obs = vel(:,depth_column) ; 
display(['num vals in Z_km_obs: ' num2str(length(Z_km_obs))])
display(['depth of Z_km_obs: ' num2str(Z_km_obs(end))])

Vs_obs = vel(:,vel_column) ; 
    
    %set(gca,'YDir','reverse');hold on;
    %plot(obsVs,depth,'r:','LineWidth', 2);hold on;

% SHAPE THE OBSERVED and CALCULATED to FIT each other... 

% obsdepth is the location of the observations, Z_km is the model z.  
[max_z,which_is_longer] = max([Z_km_obs(end) Z_km_calc(end)]) ; 

if which_is_longer==1
	disp('Obs is longer')
	Z_max = Z_km_calc(end) ;
	% get rid of any repeating z values at steps-- they kill interp1
	epsilon = 0.01 ; 
	for iz = 2:length(Z_km_obs) 
		if Z_km_obs(iz)==Z_km_obs(iz-1)
			Z_km_obs(iz) = Z_km_obs(iz-1)+epsilon ; 
		end
	end
	Z_km_obs ; 
	%  find the first element that is GREATER than Z_max of the shallower depth vector.
	i_max = find(Z_km_obs>Z_max,1);

	Z_km_obs_shorter = Z_km_obs(1:i_max) ;
	length(Z_km_obs_shorter) ;
	% interpolate to the length and depth of the shallower depth vector 
	% and n_pts of the calculation 
	% vq = interp1(x,v,xq) 
	Vs_obs_shorter = Vs_obs(1:i_max);
	length(Vs_obs_shorter) ;
	Vs_obs_new = interp1(Z_km_obs_shorter,Vs_obs_shorter,Z_km_calc) ;
		
	for iz = 1:length(Vs_obs_new) 
		if isnan(Vs_obs_new(iz))==1
			Vs_obs_new(iz)=1 ; 
		end
	end


elseif which_is_longer==2
	disp('Calc is longer')
		Z_max = Z_km_obs(end)
	% get rid of any repeating z values at steps-- they kill interp1
	epsilon = 0.01 ; 
	for iz = 2:length(Z_km_obs) 
		if Z_km_obs(iz)==Z_km_obs(iz-1)
			Z_km_obs(iz) = Z_km_obs(iz-1)+epsilon ; 
		end
	end
	Z_km_obs ; 
	%  find the first element that is GREATER than Z_max of the shallower depth vector.
	i_max = find(Z_km_calc>Z_max,1);

	Z_km_calc_shorter = Z_km_calc(1:i_max) ;
	length(Z_km_calc_shorter)
	% interpolate to the length and depth of the shallower depth vector 
	% and n_pts of the calculation
	Vs_calc_shorter = Vs_calc(1:i_max);
	length(Vs_calc_shorter)
	Vs_obs_new = interp1(Z_km_obs,Vs_obs,Z_km_calc_shorter) ;


	for iz = 1:length(Vs_obs_new) 
		if isnan(Vs_obs_new(iz))==1
			Vs_obs_new(iz)=1 ; 
		end
	end
end

display(['num vals in Vs_obs_new: ' num2str(length(Vs_obs_new))])

depth = Z_km_calc ;

obsVs = Vs_obs_new ;

%[Z_km_obs_good,i_good] = Z_km_obs(Z_km_obs<max(Z_km_calc))
%vel  = interp1(Z_km,meanV,))   


end
