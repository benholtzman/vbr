function prior = process_SeismicModels(obs_name, location, model_file)

% Load in seismic observations saved in model_file and return prior.
% 
% For a seismic observable 'obs_name' (Vs, Q, LAB), extract the observed
% value at a location gived in 'location' from the .mat file 'model_file'.
% Then calculate the prior probability distribution for that observable.
%
% Note that the prior distribution is assumed to be normal, with the 
% observation defining the mean of the distribution and the uncertainty
% on the observation defining the standard deviation of the distribution.

[obs_value, obs_error] = load_seismic_data(model_file, location, obs_name);
fprintf('\n%s:  %g, %g\n', obs_name, obs_value, obs_error)

prior = normal_probability(obs_value, obs_error);


end


function [obs_value, obs_error] = ...
    load_seismic_data(model_file, location, fieldname)

  % Load the models
    Model = load(model_file);
    varname = fieldnames(Model); % rename variable from its saved name
    Model = Model.(varname{1});  % e.g. 'Vs_Model' to just 'Model'
    
  % Check if there are any uncertainties saved in the model - if there
  % are none, put in a constant error
    Model=check_errors(Model, fieldname);
    

  % Check model coordinate overlap and then limit measurements to 
  % coordinate ranges
    Model = check_overlap(Model, location);
    Model = limit_by_coords(Model, fieldname, location);

  % Find median and error
    obs_value = find_median(Model, fieldname);
    median_error = find_median(Model, 'Error');
    lateral_error = find_lateral_error(Model, fieldname);
    %fprintf('\n%s:  %g, %g\n', fieldname, median_error, lateral_error)
    obs_error = max(median_error, lateral_error);   
    

end

function Model = check_errors(Model, field_name)
  % checks model structure for Error field. If it's not there and
  % constant_error>0, will use a constant value equal to constant_error.
  if isfield(Model, 'Error')
      return
  end
  
  switch field_name
      case 'Vs'
          constant_error = 0.05;
      case 'LAB'
          constant_error = 5;
      case 'Qinv'
          constant_error = 0.005;
  end

  fprintf(['\nError field does not exist for %s - \n\t\t', ...
      'using a constant value of %g\n'], field_name, constant_error)
  
  Model.Error = constant_error * ones(size(Model.(field_name)));

end

function [Model, returnCode] = check_overlap(Model, location)
  % make sure there's overlap in Vs and LAB models
    returnCode=1;

  % check for negative Longitudes, convert to positive
    Model.Longitude(Model.Longitude < 0) = ...
        Model.Longitude(Model.Longitude < 0) + 360;

  % global model limits
    model_lims = [min(Model.Latitude), max(Model.Latitude);...
        min(Model.Longitude), max(Model.Longitude)];

  % check that chosen coordinates are within measurements
    if ~(...
            location.lat > model_lims(1,1) ...
            && location.lat < model_lims(1,2) ...
            && location.lon > model_lims(2,1) ...
            && location.lon < model_lims(2,2) ...
         )
          disp('Your coordinates are not within your chosen model bounds')
          returnCode=0;
    end

end

function Model = limit_by_coords(Model, field_name, location)
  % limits to measurements within (lat +/- smooth_rad, lon +/- smooth_rad)

  % create lat/lon masks
    if isfield(location, 'smooth_rad')
        d_deg = location.smooth_rad;
    else
        d_deg = 0.5; % default radius over which to average seismic model 
    end
    
    lat = location.lat;
    lon = location.lon;
    lat_mask = (Model.Latitude >= lat - d_deg ...
        & Model.Latitude <= lat + d_deg);
    lon_mask = (Model.Longitude >= lon - d_deg ...
        & Model.Longitude <= lon + d_deg);
    if isfield(Model, 'Depth')
        depth_mask = (Model.Depth >= location.z_min ...
            & Model.Depth <= location.z_max);
    else
        depth_mask = 1;
    end

  % apply masks
    Model.Longitude=Model.Longitude(lon_mask);
    Model.Latitude=Model.Latitude(lat_mask);
    Model.(field_name) = Model.(field_name)(lat_mask,lon_mask,depth_mask);
    Model.Error = Model.Error(lat_mask,lon_mask,depth_mask);
    
end

function medianVal = find_median(Model, field_name)
  % calculate median value at each depth
    allVals=Model.(field_name);
    medianVal = nanmedian( ...
            reshape(allVals, 1, numel(allVals)) ...
            );
    
end

function lateral_error = find_lateral_error(Model, field_name)

    allVals = Model.(field_name);
    lateral_error = nanstd( ...
            reshape(allVals, 1, numel(allVals)) ...
            );

end

function normal_pdf = normal_probability(x, sigma)
   
    % Assume x, the observed value, coincides with mu, the mean value
    mu = x;

    normal_pdf = normpdf(x, mu, sigma);

end