function seismic_obs = get_seismic_data(Work)
% Pick your seismic models and locations, and extract relevent data
% Note, can hardwire in which models & location to use:
%       (saved in Functions_Emily/hardwired_values/)
%  - use_this_vel_model.mat:  string of file name e.g. USA_ShenRitzwoller2016
%  - use_this_lab_model.mat:  string of file name e.g. USA_HopperFischer2018
%  - use_this_location.mat:   vector [latitude, longitude]
%  - use_this_depthrange.mat: vector [min depth, max depth]

vel_name = choose_vs_model(Work);
lab_name = choose_LAB_model(Work,vel_name);

seismic_obs = getSeismicData(Work, vel_name, lab_name);


end

function use_this_vel_model = choose_vs_model(Work)

% All velocity models should be saved in Functions_Emily/vel_models 
% Following format required:
%    - filename:  [MODEL REGION]_[REFERENCE]    e.g. USA_ShenRitzwoller2016
%    - variable name: Vs_Model
%    - variable should be a structure with the following fields
%        - Name:       optional, e.g. string of reference info
%        - Latitude:   n_lat x 1 vector of latitude points
%        - Longitude:  n_lon x 1 vector of longitude points
%        - Depth:      n_dep x 1 vector of depth points
%        - Vs:         n_lat x n_lon x n_dep matrix of Vs values
%        - Error:      n_lat x n_lon x n_dep matrix of uncertainty on Vs


clc

% Look for hardwired velocity model choice
if exist([Work.valsdir '/use_this_vel_model.mat'],'file')
    load([Work.valsdir '/use_this_vel_model.mat']);
    if ~exist([Work.veldir '/' use_this_vel_model '.mat'],'file')
        disp('Velocity model not found!  Please choose saved model!');
        pause(0.5); clc; clear use_this_vel_model
    end
end

if ~exist('use_this_vel_model','var')
    % List the possible velocity models
    vels = dir(Work.veldir); clear model_names; model_names{length(vels),2} = '';
    for i_v = length(vels):-1:1
        if vels(i_v).name(1) == '.'
            vels(i_v) = []; model_names(i_v,:) = []; 
        else  model_names(i_v,:) = strsplit(vels(i_v).name,'_');
        end
    end
    
    % Pick the region of interest
    regions = unique(model_names(:,1));
    if length(regions)>1
        fprintf('\tWhich region are you interested in?\n')
        for ir = 1:length(regions)
            fprintf('\t\t%.0f.  %s\n',ir,regions{ir});
        end
        ir = input('');
        model_names = model_names(strcmp(model_names(:,1),regions{ir}));
    end
    
    % Pick the velocity model of interest
    if length(model_names)>1
        fprintf('\tWhich model are you interested in?\n')
        for im = 1:length(model_names)
            fprintf('\t\t%.0f.  %s\n',im,model_names{im,2}(1:end-4));
        end
        im = input('');
    else im = 1;
    end
    use_this_vel_model = [model_names{im,1} '_' model_names{im,2}(1:end-4)];

    fprintf('\n\nYou are using %s as your shear velocity model.\n',...
        use_this_vel_model);
    if_save = input('Save this choice for future use y/[n]?','s');
    if strcmp(if_save,'y')
        save([Work.valsdir 'use_this_vel_model.mat'], 'use_this_vel_model');
    end
end

end

function use_this_lab_model = choose_LAB_model(Work, vel_model_name)
clc
% Pick the LAB depth model

% All LAB data should be saved in Functions_Emily/LAB_models 
% Following format required:
%    - filename:  [MODEL REGION]_[REFERENCE]    e.g. USA_HopperFischer2018
%    - variable name: LAB_Model
%    - variable should be a structure with the following fields
%        - Name:       optional, e.g. string of reference info
%        - Latitude:   n_lat x 1 vector of latitude points
%        - Longitude:  n_lon x 1 vector of longitude points
%        - LAB_Depth:  n_lat x n_lon matrix of depths
%        - Error:      n_lat x n_lon matrix of uncertainty on LAB_Depth

clc

% Look for hardwired LAB model choice
if exist([Work.valsdir '/use_this_lab_model.mat'],'file')
    load([Work.valsdir '/use_this_lab_model.mat']);
    if isempty(use_this_lab_model)
        disp('Using your shear velocity model...'); pause(0.5); clc;
    elseif ~exist([Work.labdir '/' use_this_lab_model '.mat'],'file')
        disp('LAB model not found!  Please choose saved model!');
        pause(0.5); clc; clear use_this_lab_model
    end
end

if ~exist('use_this_lab_model','var')
    use_vel_model = ...
        input('Use shear velocity model to define LAB depth y/[n]?','s');
    if strcmp(use_vel_model,'y')
        if_save = input('Save this choice for future use y/[n]?','s');
        use_this_lab_model = [];
        if strcmp(if_save,'y')
            save([Work.labdir 'use_this_lab_model.mat'], 'use_this_lab_model');
        end
        return
    end


    labs = dir(Work.labdir); clear lab_names; lab_names{length(labs),2} = '';
    for i_l = length(labs):-1:1
        if labs(i_l).name(1) == '.'
            labs(i_l) = []; lab_names(i_l,:) = []; 
        else lab_names(i_l,:) = strsplit(labs(i_l).name,'_');
        end
    end

    % Pick the region of interest
    regions = unique(lab_names(:,1));
    vel_region = strsplit(vel_model_name,'_'); vel_region = vel_region{1};
    ir = find(strcmp(regions,vel_region));
    if isempty(ir)
        fprintf(['/n/tThere are no LAB models for your velocity model ', ...
            'region!\n\t\t']);
        use_vel_model = input('Just use your velocity model y/[n]?','s');
        if strcmp(use_vel_model,'y'); return; end
        
        if length(regions)>1
            fprintf('\tWhich region are you interested in?\n')
            for ir = 1:length(regions)
                fprintf('\t\t%.0f.  %s\n',ir,lab_names{ir,1});
            end
            ir = input('');
        end
    end
    lab_names = lab_names(strcmp(lab_names(:,1),regions{ir}),:);

    % Pick the velocity model of interest
    if size(lab_names,1)>1
        fprintf('\tWhich model are you interested in?\n')
        for im = 1:length(lab_names)
            fprintf('\t\t%.0f.  %s\n',im,lab_names{im,2}(1:end-4));
        end
        im = input('');
    else; im = 1;
    end
    use_this_lab_model = [lab_names{im,1} '_' lab_names{im,2}(1:end-4)];
    
    
    fprintf('\n\nYou are using %s for LAB depth.\n', use_this_lab_model);
    if_save = input('Save this choice for future use y/[n]?','s');
    if strcmp(if_save,'y')
        save([Work.valsdir 'use_this_lab_model.mat'], 'use_this_lab_model');
    end
end
   
end

function [seismic_obs] = getSeismicData(Work, vel_name, lab_name)

% Extract seismic information for a location

smooth_rad = 0.5;  % radius over which to average RFs and vel model (degrees)

load([Work.veldir vel_name '.mat']);
if isempty(lab_name)
    LAB_Model.Name = Vs_Model.Name;
    LAB_Model.Latitude = Vs_Model.Latitude;
    LAB_Model.Longitude = Vs_Model.Longitude;
else; load([Work.labdir lab_name '.mat']);
end

vs_lims = [min(Vs_Model.Latitude) max(Vs_Model.Latitude);...
    min(Vs_Model.Longitude) max(Vs_Model.Longitude)];
lab_lims = [min(LAB_Model.Latitude) max(LAB_Model.Latitude);...
    min(LAB_Model.Longitude) max(LAB_Model.Longitude)];

if (vs_lims(1,1) > lab_lims(1,2)) || (vs_lims(2,1) > lab_lims(2,2)) ...
        || (vs_lims(1,2) < lab_lims(1,1)) || (vs_lims(2,2) < lab_lims(2,1))
    clc; disp('Your velocity and LAB models don''t overlap!! Try again!');
    seismic_obs = nan; return;
end
 

% Pick the location
% Look for hardwired velocity model choice
if exist([Work.valsdir 'use_this_location.mat'],'file')
    load([Work.valsdir 'use_this_location.mat']);
     loc = use_this_location;
    
    % Is this location within the limits of our models?
    if ~(loc(1) > min([vs_lims(1,1), lab_lims(1,1)]) && ...
            loc(1) < max([vs_lims(1,2), lab_lims(1,2)]) && ...
            loc(2) > min([vs_lims(2,1), lab_lims(2,1)]) && ...
            loc(2) < max([vs_lims(2,2), lab_lims(2,2)]))
        fprintf(['Your saved location is not within your chosen model'...
            ' bounds!\n\tPlease pick a new location on the map!'])
        clear loc
    end  
end

if ~exist('loc','var'); loc = Pick_location(Vs_Model, LAB_Model, Work); end


% Extract the Vs Model
allVs = Vs_Model.Vs(Vs_Model.Latitude >= loc(1)-smooth_rad & ...
    Vs_Model.Latitude <= loc(1)+smooth_rad,...
    Vs_Model.Longitude >= loc(2)-smooth_rad...
    & Vs_Model.Longitude <= loc(2)+smooth_rad,:);
medianVs = nan(size(Vs_Model.Depth));
for id = 1:length(Vs_Model.Depth)
    medianVs(id) = median(reshape(allVs(:,:,id),1,numel(allVs(:,:,id))));
end

% Find the uncertainty
if ~isfield(Vs_Model,'Error')
    fprintf(['\n\nNo errors saved for velocity model!!\n\n\t' ...
        'Use a constant value for the error?\n\t\t(Enter a numeric value ' ...
        'if so, or any letter to quit)\n']);
    constant_error = str2double(input('','s'));
    if isnan(constant_error)
        error('Input velocity model should have error field!');
    else
        medianVs_error = constant_error*ones(size(medianVs));
    end
else
    allVs_error = Vs_Model.Error(Vs_Model.Latitude >= loc(1)-smooth_rad & ...
        Vs_Model.Latitude <= loc(1)+smooth_rad,...
        Vs_Model.Longitude >= loc(2)-smooth_rad...
        & Vs_Model.Longitude <= loc(2)+smooth_rad,:);
    medianVs_error = nan(size(Vs_Model.Depth));
    for id = 1:length(Vs_Model.Depth)
        medianVs_error(id) = median(reshape(...
            allVs_error(:,:,id),1,numel(allVs_error(:,:,id))));
    end
end


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

% Find the LAB
if isfield(LAB_Model, 'LAB_Depth')
    allLABs = LAB_Model.LAB_Depth(LAB_Model.Latitude >= loc(1)-smooth_rad & ...
        LAB_Model.Latitude <= loc(1)+smooth_rad,...
        LAB_Model.Longitude >= loc(2)-smooth_rad...
        & LAB_Model.Longitude <= loc(2)+smooth_rad);
    LAB = nanmedian(allLABs(:));
    
    if isnan(LAB); LAB = find_LAB_from_Vs(allVs, Moho, Vs_Model.Depth); end
else;  LAB = find_LAB_from_Vs(allVs, Moho, Vs_Model.Depth);
end


% Look for hardwired depth range
if exist([Work.valsdir 'use_this_depthrange.mat'],'file')
    load([Work.valsdir 'use_this_depthrange.mat']);
    depthrange = use_this_depthrange;
else depthrange = pick_depth_range(medianVs, Vs_Model.Depth, Moho, LAB, Work);
end

asth_v = mean(medianVs(Vs_Model.Depth >= depthrange(1) & ...
    Vs_Model.Depth <= depthrange(2)));
asth_v_error = mean(medianVs_error(Vs_Model.Depth >= depthrange(1) & ...
    Vs_Model.Depth <= depthrange(2)));



fprintf(['\n\nAverage Vs (%.0f - %.0f km)' ...
    ': %.2f %s %.2f km/s\n'],depthrange(1), depthrange(2), asth_v,...
    char(177), asth_v_error);


seismic_obs.asth_v = asth_v;
seismic_obs.asth_v_error = asth_v_error;
seismic_obs.medianVs = medianVs;
seismic_obs.medianVs_error = medianVs_error;
seismic_obs.depthrange = depthrange;
seismic_obs.Moho = Moho;
seismic_obs.LAB = LAB;

    
end

function loc = Pick_location(Vs_Model, LAB_Model, Work)
close all
% Find dataset limits
vs_lims = [min(Vs_Model.Latitude) max(Vs_Model.Latitude);...
    min(Vs_Model.Longitude) max(Vs_Model.Longitude)];
lab_lims = [min(LAB_Model.Latitude) max(LAB_Model.Latitude);...
    min(LAB_Model.Longitude) max(LAB_Model.Longitude)];
lims = [max([vs_lims(:,1), lab_lims(:,1)],[],2)-1,...
    min([vs_lims(:,2), lab_lims(:,2)],[],2)+1];

coast = load('coast');

figure('color','w','position',[50 50 900 600]);
for iff = 1:2
    if iff == 1
        a1 = subplot(1,2,1); hold on;
        [~,id] = min(abs(Vs_Model.Depth - 100));
        contourf(Vs_Model.Longitude, Vs_Model.Latitude, Vs_Model.Vs(:,:,id),...
            'linestyle','none')
        cstr = 'Shear Velocity at 100 km depth (km/s)';
    else
        a2 = subplot(1,2,2); hold on;
        [glon,glat] = meshgrid(LAB_Model.Longitude, LAB_Model.Latitude);
        inds = find(~isnan(LAB_Model.LAB_Depth));
        scatter(glon(inds),glat(inds),50, LAB_Model.LAB_Depth(inds), ...
            's','filled');
        cstr = 'LAB Depth (km)';
    end
    
    plot(coast.long,coast.lat,'k-')
    daspect([111.16, 111.16*distance(mean(lims(1,:)),0,mean(lims(1,:)),1), 1]);
    xlim(lims(2,:)); ylim(lims(1,:));
    box on; set(gca,'layer','top'); 
    c=colorbar('location','southoutside'); colormap(copper);
    xlabel(c,cstr);
end

clc
yn = input('Zoom in y/[n]? ','s');
if strcmp(yn,'y')
    fprintf(['\n\tPick the northwest and southeast limits of region \n' ...
        'to zoom in on'])
    [lon,lat] = ginput(2);
    set(a1, 'xlim', sort(lon), 'ylim', sort(lat));
    set(a2, 'xlim', sort(lon), 'ylim', sort(lat));
end

fprintf('\n\nChoose a good location on the map!\n\n');
[lon,lat] = ginput(1);

loc = [lat, lon];

if_save = input('Save this choice for future use y/[n]?','s');
if strcmp(if_save,'y')
    use_this_location = loc;
    save([Work.valsdir 'use_this_location.mat'], 'use_this_location');
end

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

figure;
plot(median(allVs),Depth); axis ij; hold on
xl = get(gca,'xlim');
plot(xl, LAB*[1 1],'b--');
pause; close;



end

function asth_v_deps = pick_depth_range(medianVs, Depth, Moho, LAB, Work)

clc; close all
fprintf(['\n\nChoose the depth range over which you want to model\n' ...
    'asthenospheric velocities...\n']);

figure; fp = get(gcf, 'position'); 
fp(2) = 50; fp(3:4) = fp(3:4).*[1/1.4,1.4];
set(gcf, 'position', fp);

% Plot the velocity profile
axes('position',[0.15 0.15 0.8 0.8]); 
plot(medianVs,Depth,'r-','linewidth',2); axis ij; hold on
xlabel('Shear Velocity (km/s)'); ylabel('Depth (km)');
xl = get(gca,'xlim');
if xl(1) < 3.8; xl(1) = 3.8; end; xlim(xl);
yl = get(gca,'ylim');
if yl(2) > 300; yl(2) = 350; end; ylim(yl);

% Plot on the Moho
plot(xl,Moho*[1 1],'--','color',0.6*[1 1 1]);
text(mean([xl,xl(1)]), Moho - 5, 'Moho'); 

% Plot on the LAB
plot(xl,LAB*[1 1],'b--');
text(mean([xl,xl(1)]), LAB - 5, 'LAB'); 

title('Pick your depth range!');
[~, asth_v_deps] = ginput(2); 
asth_v_deps = sort(asth_v_deps);


% Check correct depths chosen
while 1 == 1
    p1 = patch(xl([1 2 2 1 1]), asth_v_deps([1 1 2 2 1]),'r');
    alpha(0.3);
    yn = input('Keep this depth range y/[n]?  ', 's');
    if strcmp(yn,'y'); break; end
    delete(p1);
    [~, asth_v_deps] = ginput(2);
    asth_v_deps = sort(asth_v_deps);
end

clc
fprintf(['\n\nYou have chosen %.0f-%.0f km depth.\n\t%.0f to %.0f km '...
    'below the Moho\n\t%.0f to %.0f km below the LAB\n\n'], asth_v_deps(1), ...
    asth_v_deps(2), asth_v_deps(1) - Moho, asth_v_deps(2) - Moho, ...
    asth_v_deps(1) - LAB, asth_v_deps(2) - LAB);

if_save = input('Save this depth range for future use y/[n]?','s');
if strcmp(if_save,'y')
    use_this_depthrange = asth_v_deps;
    save([Work.valsdir 'use_this_depthrange.mat'], 'use_this_depthrange');
end
    

end