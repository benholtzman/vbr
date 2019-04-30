function Observations = process_SeismicModels(Files,Coords)

  % load the models
    load(Files.Vs_Model_file);
    Vs_Model=checkErrors(Vs_Model,'Vs',0.01);
    if isfield(Files,'LAB_Model_file')
      load(Files.LAB_Model_file);
      LAB_Model=checkErrors(LAB_Model,'LAB_Depth',1);
    else
      LAB_Model=struct();
    end

  % check model coordinate overlap and then limit measurements to coordinate ranges
    [LAB_Model,Vs_Model,rCode]=checkOverlap(LAB_Model,Vs_Model,Coords);
    Vs_Model=limitByCoords(Vs_Model,'Vs',Coords);
    if ~isfield(Files,'LAB_Model_file')
      LAB_Model=limitByCoords(LAB_Model,'LAB_Depth',Coords);
    end

  % find median Vs, Vs Error
    Vs_Model=findMedian(Vs_Model,'Vs');
    Vs_Model=findMedian(Vs_Model,'Error');

  % find Moho, LAB from Vs profile
    Min_Moho=20;  Max_Moho=60;
    Moho=findMoho(Vs_Model,Min_Moho,Max_Moho) ;

  % find median LAB depth
    LAB_median = median(LAB_Model.LAB_Depth(~isnan(LAB_Model.LAB_Depth)));
  %   % LAB_from_Vs=find_LAB_from_Vs(Vs_Model); % (allVs, Moho, Depth)

  % find asthenosphere velocity
  %   % Coords  = pick_depth_range(medianVs, Vs_Model.Depth, Moho, LAB, Work);
    depth_mask=(Vs_Model.Depth >= Coords.z_min & Vs_Model.Depth <= Coords.z_max);
    asth_v = mean(Vs_Model.Vs_median(depth_mask));
    asth_v_error = mean(Vs_Model.Error_median(depth_mask));

  % % format final structure
  Observations.asth_v = asth_v;
  Observations.asth_v_error = asth_v_error;
  Observations.medianVs = Vs_Model.Vs_median;
  Observations.medianVs_error = Vs_Model.Error_median;
  Observations.depth = Vs_Model.Depth;
  Observations.depthrange = [Coords.z_min Coords.z_max];
  Observations.Moho = Moho;
  Observations.LAB = LAB_median;
  Observations.Vs_Model=Vs_Model;
  Observations.LAB_Model=LAB_Model;

end

function [LAB_Model,Vs_Model,returnCode] = checkOverlap(LAB_Model,Vs_Model,Coords)
  % make sure there's overlap in Vs and LAB models
    returnCode=1;

  % check for negative Longitudes, convert to positive
    Vs_Model.Longitude(Vs_Model.Longitude<0)=Vs_Model.Longitude(Vs_Model.Longitude<0)+360;
    if isfield(LAB_Model,'LAB_Depth')
      LAB_Model.Longitude(LAB_Model.Longitude<0)=LAB_Model.Longitude(LAB_Model.Longitude<0)+360;
    end

  % global Vs limits
    vs_lims = [min(Vs_Model.Latitude) max(Vs_Model.Latitude);...
        min(Vs_Model.Longitude) max(Vs_Model.Longitude)];
    global_lims=vs_lims;

  % check overlap between Vs and LAB models
    if isfield(LAB_Model,'LAB_Depth')
      lab_lims = [min(LAB_Model.Latitude) max(LAB_Model.Latitude);...
          min(LAB_Model.Longitude) max(LAB_Model.Longitude)];
      if (vs_lims(1,1) > lab_lims(1,2)) || (vs_lims(2,1) > lab_lims(2,2)) ...
              || (vs_lims(1,2) < lab_lims(1,1)) || (vs_lims(2,2) < lab_lims(2,1))
        disp('Your velocity and LAB models don''t overlap!! Try again!');
        returnCode=0;
      else
        for i=1:2
          for j=1:2
            global_lims(i,j)=min(global_lims(i,j),lab_lims(i,j));
          end
        end
      end
    end

  % check that chosen coordinates are within measurements
    if ~(Coords.lat > global_lims(1,1) && Coords.lat < global_lims(1,2) ...
          && Coords.lon > global_lims(2,1) && Coords.lon < global_lims(2,2))
          disp(['Your coordinates are not within your chosen model bounds'])
          returnCode=0
    end

end

function Model = checkErrors(Model,field_name,constant_error)
  % checks model structure for Error field. If it's not there and
  % constant_error>0, will use a constant value equal to constant_error.
  if ~isfield(Model,'Error')
    disp(['Error field does not exist for ',field_name])
    if constant_error > 0
      disp(['using at constant value of ',num2str(constant_error)])
      Model.Error=constant_error * ones(size(Model.(field_name)));
    end
  end
end

function Model = limitByCoords(Model,field_name,Coords)
  % limits to measurements within (lat +/- smooth_rad, lon +/- smooth_rad)

  % create lat/lon masks
    d_deg=Coords.smooth_rad;
    lat=Coords.lat;
    lon=Coords.lon;
    lat_mask=(Model.Latitude >= lat-d_deg & Model.Latitude <= lat+d_deg);
    lon_mask=(Model.Longitude >= lon-d_deg & Model.Longitude <= lon+d_deg);

  % apply masks
    Model.Longitude=Model.Longitude(lon_mask);
    Model.Latitude=Model.Latitude(lat_mask);
    Model.(field_name) = Model.(field_name)(lat_mask,lon_mask,:);
    if isfield(Model,'Error')
      Model.Error = Model.Error(lat_mask,lon_mask,:);
    end

end

function Model = findMedian(Model,field_name)
  % calculate median value at each depth
    allVals=Model.(field_name);
    medianVal = nan(size(Model.Depth));
    for id = 1:length(Model.Depth)
        medianVal(id) = median(reshape(allVals(:,:,id),1,numel(allVals(:,:,id))));
    end
    Model.([field_name,'_median'])=medianVal;
end

function Moho = findMoho(Vs_Model,min_depth,max_depth)
  % Find the Moho: depth is location of max positive gradient in Vs

  % find Moho for every lat/lon
  moho_dep_range =  find(Vs_Model.Depth>min_depth & Vs_Model.Depth<max_depth);
  allVs = Vs_Model.Vs;
  allVs = reshape(allVs,size(allVs,1)*size(allVs,2),size(allVs,3));
  Moho = nan(size(Vs_Model.Vs(:,:,1)));
  for iv = 1:size(allVs,1)
     [~,i_moho] = max(diff(allVs(iv,moho_dep_range)));
     if ~isempty(i_moho)
         Moho(iv) = Vs_Model.Depth(moho_dep_range(1) + i_moho(1));
     end
  end

  % find median of that value
  Moho = median(Moho);
end


% function  LAB = find_LAB_from_Vs(allVs, Moho, Depth)
%
%   % Identify the LAB as the maximum negative velocity gradient below the Moho
%   mantle_inds = find(Depth > Moho & Depth < 300);
%
%   LABs = nan(size(allVs,1),1);
%
%   for i_v = 1:size(allVs,1)
%       % See if Vs profile has nice shape (i.e. high velocity lid with LVZ below)
%       [~,i_high_vel] = findpeaks(allVs(i_v,mantle_inds), ...
%           'minpeakprominence',0.05);
%       [~,i_low_vel] = findpeaks(-allVs(i_v,mantle_inds), ...
%           'minpeakprominence',0.05);
%
%
%       if isempty(i_high_vel); i_high_vel = 1; end
%       i_low_vel(i_low_vel<i_high_vel(1)) = []; % remove any lows beneath first peak
%       if isempty(i_low_vel);  i_low_vel = length(mantle_inds); end
%
%
%       % Find the steepest negative velocity gradient in this depth range
%       [~, i_max_nvg] = max(-diff(allVs(i_v,...
%           mantle_inds(i_high_vel(1):i_low_vel(1)))));
%
%       LABs(i_v) = Depth(mantle_inds(i_high_vel(1)+i_max_nvg -1));
%   end
%
%   LAB = nanmedian(LABs);
%
% end
