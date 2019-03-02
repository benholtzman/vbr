function Observations = process_SiesmicModels(Files)

  % load the models
  LAB_Model=load(Files.LAB_Model_file);
  Vs_Model=load(Files.Vs_Model_file);

  Coords=struct();
  Coords.lat=30;
  Coords.lon=-100;
  Coords.smooth_rad = 0.5;  % radius over which to average RFs and vel model (degrees)
  Coords.z_min=100;
  Coords.z_max=150;

  % check model coordinate overlap and then limit measurements to coordinate ranges
    [LAB_Model,Vs_Model]=checkOverlap(LAB_Model,Vs_Model,Coords);
    Vs_Model=limitByCoords(Vs_Model,'Vs',Coords);
    LAB_Model=limitByCoords(LAB_Model,'LAB_Depth',Coords);

  % find median Vs, Vs Error, LAB Depth and LAB Depth error
    Vs_Model=findMedian(Vs_Model,'Vs');
    Vs_Model=findMedian(Vs_Model,'Error');
    LAB_Model=findMedian(LAB_Model,'LAB_Depth');

    LAB_Model=findMedian(LAB_Model,'Error');


  % find Moho, LAB from Vs profile
    Moho=findMoho(Vs_Model) ;
    LAB_from_Vs=find_LAB_from_Vs(Vs_Model); % (allVs, Moho, Depth)

  % find asthenosphere velocity
    % Coords  = pick_depth_range(medianVs, Vs_Model.Depth, Moho, LAB, Work);
    depth_mask=(Vs_Model.Depth >= Coords.z_min & Vs_Model.Depth <= Coords.z_max);
    asth_v = mean(VS_model.Vs_median(depth_mask));
    asth_v_error = mean(VS_model.Error_median(depth_mask));

  % format final structure
  Observations.asth_v = asth_v;
  Observations.asth_v_error = asth_v_error;
  Observations.medianVs = medianVs;
  Observations.medianVs_error = medianVs_error;
  Observations.depth = Vs_Model.Depth;
  Observations.depthrange = depthrange;
  Observations.Moho = Moho;
  Observations.LAB = LAB;
  Observations.Vs_Model=Vs_Model;
  Observations.LAB_MOdel=LAB_Model;

end

function returnCode = checkOverlap(LAB_Model,Vs_Model)
  % maker sure there's overlap in Vs and LAB models
  returnCode=1;
  vs_lims = [min(Vs_Model.Latitude) max(Vs_Model.Latitude);...
      min(Vs_Model.Longitude) max(Vs_Model.Longitude)];
  lab_lims = [min(LAB_Model.Latitude) max(LAB_Model.Latitude);...
      min(LAB_Model.Longitude) max(LAB_Model.Longitude)];
  if (vs_lims(1,1) > lab_lims(1,2)) || (vs_lims(2,1) > lab_lims(2,2)) ...
          || (vs_lims(1,2) < lab_lims(1,1)) || (vs_lims(2,2) < lab_lims(2,1))
    disp('Your velocity and LAB models don''t overlap!! Try again!');
    returnCode=0;
  end

  if ~(Coords.lat > min([vs_lims(1,1), lab_lims(1,1)]) && ...
       Coords.lat < max([vs_lims(1,2), lab_lims(1,2)]) && ...
       Coords.lon > min([vs_lims(2,1), lab_lims(2,1)]) && ...
       Coords.lon < max([vs_lims(2,2), lab_lims(2,2)]))
          disp(['Your coordinates are not within your chosen model bounds'])
          returnCode=0
  end

  if returnCode==1

  end
end

function Model = limitByCoords(Model,field_name,Coords)

    smooth_rad=Coords.smooth_rad;
    lat=Coords.lat;
    lon=Coords.lon;

  % limits to measurements within (lat +/- smooth_rad, lon +/- smooth_rad)
    Model.(fieldVal) = Model.(fieldVal)(...
      Model.Latitude >= lat-smooth_rad & Model.Latitude <= lat+smooth_rad,...
      Model.Longitude >= lon-smooth_rad & Model.Longitude <= lon+smooth_rad,...
      :);

    Model.Error = Model.Error(...
       Model.Latitude >= lat-smooth_rad & Model.Latitude <= lat+smooth_rad,...
       Model.Longitude >= lon-smooth_rad & Model.Longitude <= lon+smooth_rad,...
       :);

    Model.Longitude=Model.Longitude(Model.Longitude >= lon-smooth_rad & ...
                                    Model.Longitude <= lon+smooth_rad)
    Model.Latitude=Model.Latitude(Model.Latitude >= lat-smooth_rad & ...
                                    Model.Latitude <= lat+smooth_rad)
end


function Model = findMedian(Model,field_name)
  % calculate median value at each depth
  allVals=Model.(field_name)
  medianVal = nan(size(Model.Depth));
  for id = 1:length(Model.Depth)
      medianVal(id) = median(reshape(allVals(:,:,id),1,numel(allVals(:,:,id))));
  end
  Model.([field_name,'_median'])=medianVal;
end

function Moho = findMoho(Vs_Model)

  % Find the Moho
  allVs = reshape(allVs,size(allVs,1)*size(allVs,2),size(allVs,3));
  Moho = nan(size(allVs,1),1);
  moho_dep_range =  find(Vs_Model.Depth>20 & Vs_Model.Depth<60);
  for iv = 1:size(allVs,1)
     [~,i_moho] = max(diff(allVs(iv,moho_dep_range)));
     if ~isempty(i_moho)
         Moho(iv) = Vs_Model.Depth(moho_dep_range(1) + i_moho(1));
     end
  end
  Moho = median(Moho);
end


function  LAB = find_LAB_from_Vs(allVs, Moho, Depth)

% Identify the LAB as the maximum negative velocity gradient below the Moho
mantle_inds = find(Depth > Moho & Depth < 300);

LABs = nan(size(allVs,1),1);

for i_v = 1:size(allVs,1)
    % See if Vs profile has nice shape (i.e. high velocity lid with LVZ below)
    [~,i_high_vel] = findpeaks(allVs(i_v,mantle_inds), ...
        'minpeakprominence',0.05);
    [~,i_low_vel] = findpeaks(-allVs(i_v,mantle_inds), ...
        'minpeakprominence',0.05);


    if isempty(i_high_vel); i_high_vel = 1; end
    i_low_vel(i_low_vel<i_high_vel(1)) = []; % remove any lows beneath first peak
    if isempty(i_low_vel);  i_low_vel = length(mantle_inds); end


    % Find the steepest negative velocity gradient in this depth range
    [~, i_max_nvg] = max(-diff(allVs(i_v,...
        mantle_inds(i_high_vel(1):i_low_vel(1)))));

    LABs(i_v) = Depth(mantle_inds(i_high_vel(1)+i_max_nvg -1));
end

LAB = nanmedian(LABs);

end
