function [probs] = multivariate_Bayesian_inversion(...
    Work, seismic_obs, zPlate, ifnormal)


% Parameter sweep - define ranges for phi and grain size!
sweep.phi_range = (0.0:0.005:0.03).*1e2; % melt fraction
sweep.phi_units = ' %';
sweep.gs_range = [0.5 1 5:5:100].*1e3; % grain size
sweep.gs_units = ' microns';

sweepBox = generate_parameter_sweep(Work, zPlate, sweep, seismic_obs);
[vs_vals, residual] = ...
    extract_Vs(sweepBox, seismic_obs);


% Bayes' theorem:  P(A|B) = ( P(B|A)P(A) ) / P(B)
%     In this case, A is [phi, Tp, grain size];  B is observed Vs
%
% P(T, phi, gs | Vs) = __P(_Vs_|_T,_phi,_gs)_*_P(_T,_phi,_gs_)___
%                                        P(Vs)
% (the RHS there is a fraction, btw...)

% Let's calculate some probabilities!
% Hardwire all PDFs in this function
% ifnormal:  0: uniform probability; 1: normal distribution
probs = calculate_probabilities(sweepBox, seismic_obs, residual, ifnormal);

% And do some plotting!
plot_probs(probs, seismic_obs, vs_vals, sweepBox);

end

function sweepBox = generate_parameter_sweep(Work, zPlate, sweep, seismic_obs)


Work.boxname = ['Box_' Work.Box_base_name '_VBR_py.mat'];
Work.boxpath = [Work.hmdir 'Boxes/' Work.Box_base_name '/'];
load([Work.boxpath Work.boxname]);
Work.sweep_boxname = ['Box_' Work.Box_base_name '_sweep_' ...
    num2str(zPlate.zPlate) 'km.mat'];
clc
fprintf('\n\nPlate thickness: %g\n', ...
    Box(1,1).info.var2range(zPlate.zPlate_ind))

if exist([Work.boxpath Work.sweep_boxname],'file')
    fprintf('\nYou already have a saved parameter sweep\n\t%s\n', ...
        Work.sweep_boxname(1:end-4));
    if_recalculate = input('Recalculate y/[n]? ', 's');

    if ~strcmp(if_recalculate,'y')
        load([Work.boxpath Work.sweep_boxname])
        return
    end
end

% For a given zPlate, we want to run through all combinations of Tpot, melt
% fraction, and grain size.
sweep.Tpot_range = Box(1,1).info.var1range;
sweep.Tpot_units = Box(1,1).info.var1units;
sweep.zPlate = Box(1,1).info.var2range(zPlate.zPlate_ind);
sweep.zPlate = Box(1,1).info.var2units;

n_Tpot   = numel(sweep.Tpot_range);
n_phi    = numel(sweep.phi_range);
n_gs     = numel(sweep.gs_range);

tstart = cputime;

for i_Tpot = n_Tpot:-1:1  % run backwards so structure is preallocated

    Tbox = Box(i_Tpot, zPlate.zPlate_ind);
    VBR_in = Tbox.Frames(end).VBR;

    % Extract relevant parameters from the initial run
    sweep.zLAB = Tbox.run_info.zLAB(end);
    sweep.Z_km = Tbox.run_info.Z_km;
    sweep.T_z  = Tbox.Frames(end).T(:,1);
    sweep.Tpot = sweep.Tpot_range(i_Tpot);

    compare_ind = round(mean(find(sweep.Z_km >= seismic_obs.depthrange(1) ...
        & sweep.Z_km <= seismic_obs.depthrange(2))));


    for i_phi = n_phi:-1:1

        sweep.phi = sweep.phi_range(i_phi);

        for i_gs = n_gs:-1:1

            fprintf(['-----------------------------------------\n' ...
                'Run %g of %g\n'], n_gs*n_phi*n_Tpot + 1 - ...
                (i_gs + (i_phi-1)*n_gs + (i_Tpot-1)*n_phi*n_gs),...
                n_gs*n_phi*n_Tpot);

            sweep.gs = sweep.gs_range(i_gs);

            % Put values in the VBR and run it for asthenospheric depth
            VBR.in = VBR_in.in;
            fn = fieldnames(VBR.in.SV);
            for iff = 1:length(fn)
                if size(VBR.in.SV.(fn{iff}),1) == length(sweep.Z_km)
                    VBR.in.SV.(fn{iff}) = VBR.in.SV.(fn{iff})(compare_ind);
                end
            end
            VBR.in.SV.phi = sweep.phi.*1e-2;
            VBR.in.SV.dg_um = sweep.gs;

            VBR = VBR_spine(VBR);
            anelastic_methods = fieldnames(VBR.out.anelastic);
            if isfield(VBR.out.anelastic,'YT_maxwell')
                VBR.out.anelastic.YT_maxwell.Vave = ...
                    VBR.out.anelastic.YT_maxwell.V;
            end
            for i_an = 1:length(anelastic_methods)
                vbr.(anelastic_methods{i_an}).Vave = ...
                    VBR.out.anelastic.(anelastic_methods{i_an}).Vave;
                vbr.(anelastic_methods{i_an}).Qinv = ...
                   VBR.out.anelastic.(anelastic_methods{i_an}).Qinv;
            end

            sweepBox(i_Tpot, i_phi, i_gs).info = orderfields(sweep);
            sweepBox(i_Tpot, i_phi, i_gs).VBR = vbr;

        end
    end
    %save([Work.boxpath Work.sweep_boxname], 'sweepBox');
end
tend = cputime;

save([Work.boxpath Work.sweep_boxname], 'sweepBox');

fprintf(['-----------------------------------------\n' ...
    'Computations complete!\n\tTotal VBR time: %.1f min\n' ...
    'Box saved to %s.'],(tend-tstart)/60, Work.sweep_boxname);

end

function [vs_vals, residual] = extract_Vs(sweepBox, seismic_obs)

vs_vals = zeros(size(sweepBox));

for k = 1:numel(vs_vals)

    if strcmp(seismic_obs.q_method, 'YT_maxwell')
        vs_vals(k) = mean(sweepBox(k).VBR.(seismic_obs.q_method).Vave.*1e-3);
    else
        vs_vals(k) = sweepBox(k).VBR.(seismic_obs.q_method).Vave.*1e-3;
    end

end

residual = (vs_vals - seismic_obs.asth_v).^2 / seismic_obs.asth_v;

end

function probs = calculate_probabilities(sweepBox, seismic_obs, ...
    residual, ifnormal)

% Remember the PDF for a normal distribution
%    PDF  = 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2))


P_mod = zeros(size(sweepBox));

% P(Vs)
% This is the probability of the observed Vs being correct
% We'll assume a normal distribution, mean = observed Vs, std = obs error
sigma = seismic_obs.asth_v_error;
P_Vs = 1/sqrt(2*pi*sigma^2); % assume mean = observed Vs so (x-mu) = 0
sigmaVs = median(seismic_obs.medianVs_error);

for k = 1:numel(sweepBox)

    values = sweepBox(k).info;

    % P(T, phi, gs)  == P(T) * P(phi) * P(gs)
    % Assume that all of these are independent of one another.......
    % And have pretty broad normal distributions - can either hard wire in
    % values or guesstimate from the size of the box calculated


    if ifnormal

        % P(T)
        sigmaT = abs(diff(values.Tpot_range([1 end])));
        mu     = mean(values.Tpot_range);
        x      = values.Tpot;
        P_T    = 1/sqrt(2*pi*sigmaT^2) * exp((-(x-mu)^2)/(2*sigmaT^2));

        % P(phi)
        sigmaP = abs(diff(values.phi_range([1 end])));
        mu     = mean(values.phi_range);
        x      = max(values.phi);
        P_phi  = 1/sqrt(2*pi*sigmaP^2) * exp((-(x-mu)^2)/(2*sigmaP^2));

        % P(grain size)
        sigmaG = abs(diff(values.gs_range([1 end])));
        mu     = mean(values.gs_range);
        x      = values.gs;
        P_gs   = 1/sqrt(2*pi*sigmaG^2) * exp((-(x-mu)^2)/(2*sigmaG^2));

    else

        P_T   = 1;
        P_phi = 1;
        P_gs  = 1;

    end

    % P(Vs | T, phi, gs)
    % We have a residual between observed and calculated Vs
    % So... probability = 1/sqrt(2*pi*residual) * exp(-residual/2)
    % See manual, Menke book Ch 11
    P_Vs_given_mod = (2*pi*residual(k)).^-0.5 .* exp(-0.5*residual(k));


    % All together now!
    P_mod(k) = (P_Vs_given_mod * P_T * P_phi * P_gs)./(P_Vs);


end

norm_for_residual = (2*pi*0.0001).^-0.5 .* exp(-0.5*0.0001);


if ifnormal
    probs.P_mod = P_mod .* (2*pi*sigmaT*sigmaP*sigmaG/sigmaVs) ...
        / norm_for_residual;
else
    probs.P_mod = P_mod ./ sqrt(2*pi*sigmaVs^2) / norm_for_residual;
end
probs.Tpot  = sweepBox(k).info.Tpot_range;
probs.phi   = sweepBox(k).info.phi_range;
probs.gs    = sweepBox(k).info.gs_range;
probs.ifnormal = ifnormal;

end


function plot_probs(probs, seismic_obs, vs_vals, sweepBox)

% probs.P_mod = probs.P_mod(:,:,4:end);
% probs.P_mod = probs.P_mod./max(probs.P_mod(:));
% probs.gs = probs.gs(4:end);


[~, i_max] = max(probs.P_mod(:));
[i_Tpot, i_phi, i_gs] = ind2sub(size(probs.P_mod),i_max);


figure('color','w','position',[100 100 1100 600]);

axes('visible','off','position',[0 0 1 1]);
text(0.47, 0.95, strrep(seismic_obs.q_method,'_',' '),...
    'fontweight','bold','fontsize',24);

tstr = 'Potential Temperature (\circC)';
phistr = 'Melt Fraction (%)';
gstr = 'Grain size (mm)';
probs.gs = probs.gs.*1e-3;
ms = 10;

% Plot probabilities
cols = [0 114 178; ... % best fit: blue
    86 180 233; ... % smaller value: light blue
    0 158 115]./255; % larger value: blue green
%contour_levels = [50 70 90 95 99];
if probs.ifnormal
    contour_levels = [99 99.5 99.9];
else
    contour_levels = [99 99.5 99.9];
end


% Probability in Tpot/phi space at constant gs
axes('position',[0.1 0.35 0.225 0.55],'layer','top'); hold on; box on;
xlabel(tstr); ylabel(phistr);
% Find other grain sizes to plot - test extremes first
% Small grain size
i_small = 1;
while (max(max(probs.P_mod(:,:,i_small))) < contour_levels(1)/100 ...
        && i_small < length(probs.gs)) ||...
        (i_small == i_gs && i_small < length(probs.gs))
    i_small = i_small + 1;
end
contour_plot(probs.P_mod(:,:,i_small)', probs.Tpot, probs.phi, ...
    contour_levels, cols(2,:))
ps = probs.P_mod(:,:,i_small); [~, s_i_max] = max(ps(:));
[s_i_Tpot, s_i_phi] = ind2sub(size(ps),s_i_max);
% Large grain size
i_large = length(probs.gs);
while (max(max(probs.P_mod(:,:,i_large))) < contour_levels(1)/100 ...
        && i_large > 1) || (i_large == i_gs && i_large > 1)
    i_large = i_large - 1;
end
contour_plot(probs.P_mod(:,:,i_large)', probs.Tpot, probs.phi, ...
    contour_levels, cols(3,:))
pl = probs.P_mod(:,:,i_large); [~, l_i_max] = max(pl(:));
[l_i_Tpot, l_i_phi] = ind2sub(size(pl),l_i_max);

% Best fitting grain size
contour_plot(probs.P_mod(:,:,i_gs)', probs.Tpot, probs.phi, ...
    contour_levels, cols(1,:))

plot(probs.Tpot(s_i_Tpot),probs.phi(s_i_phi),'kp','markerfacecolor',cols(2,:));
plot(probs.Tpot(l_i_Tpot),probs.phi(l_i_phi),'kp','markerfacecolor',cols(3,:));
plot(probs.Tpot(i_Tpot),probs.phi(i_phi),'kp','markerfacecolor',cols(1,:));

axes('position',[0.1 0.15 0.225 0.1]); hold on;
plot(probs.gs, zeros(size(probs.gs)), 'k-');
plot(probs.gs(i_small),0,'ko','markerfacecolor',cols(2,:),'markersize',ms);
plot(probs.gs(i_large),0,'ko','markerfacecolor',cols(3,:),'markersize',ms);
plot(probs.gs(i_gs),0,'ko','markerfacecolor',cols(1,:),'markersize',ms);
xlabel(gstr); ylim([0 eps]); daspect([1 1 1])



% Probability in Tpot/gs space at constant phi
axes('position',[0.4 0.35 0.225 0.55],'layer','top'); hold on; box on;
xlabel(tstr); ylabel(gstr);
% Find other grain sizes to plot - test extremes first
% Small grain size
i_small = 1;
while (max(max(squeeze(probs.P_mod(:,i_small,:)))) < contour_levels(1)/100 ...
        && i_small < length(probs.phi)) || ...
        (i_small == i_phi && i_small < length(probs.phi))
    i_small = i_small + 1;
end
contour_plot(squeeze(probs.P_mod(:,i_small,:))', probs.Tpot, probs.gs, ...
    contour_levels, cols(2,:))
ps = squeeze(probs.P_mod(:,i_small,:)); [~, s_i_max] = max(ps(:));
[s_i_Tpot, s_i_gs] = ind2sub(size(ps),s_i_max);

% Large grain size
i_large = length(probs.phi);
while (max(max(squeeze(probs.P_mod(:,i_large,:)))) < contour_levels(1)/100 ...
        && i_large > 1) || ...
        (i_large == i_phi && i_large > 1)
    i_large = i_large - 1;
end
contour_plot(squeeze(probs.P_mod(:,i_large,:))', probs.Tpot, probs.gs, ...
    contour_levels, cols(3,:))
pl = squeeze(probs.P_mod(:,i_large,:)); [~, l_i_max] = max(pl(:));
[l_i_Tpot, l_i_gs] = ind2sub(size(pl),l_i_max);

% Best fitting grain size
contour_plot(squeeze(probs.P_mod(:,i_phi,:))', probs.Tpot, probs.gs, ...
    contour_levels, cols(1,:))

plot(probs.Tpot(s_i_Tpot),probs.gs(s_i_gs),'kp','markerfacecolor',cols(2,:));
plot(probs.Tpot(l_i_Tpot),probs.gs(l_i_gs),'kp','markerfacecolor',cols(3,:));
plot(probs.Tpot(i_Tpot),probs.gs(i_gs),'kp','markerfacecolor',cols(1,:));



axes('position',[0.4 0.15 0.225 0.1]); hold on
plot(probs.phi, zeros(size(probs.phi)), 'k-');
plot(probs.phi(i_small),0,'ko','markerfacecolor',cols(2,:),'markersize',ms);
plot(probs.phi(i_large),0,'ko','markerfacecolor',cols(3,:),'markersize',ms);
plot(probs.phi(i_phi),0,'ko','markerfacecolor',cols(1,:),'markersize',ms);
xlabel(phistr); ylim([0 eps]); daspect([1 1 1])





% Probability in phi/gs space at constant Tpot
axes('position',[0.7 0.35 0.225 0.55],'layer','top'); hold on; box on;
xlabel(phistr); ylabel(gstr);
% Find other grain sizes to plot - test extremes first
% Small grain size
i_small = 1;
while (max(max(squeeze(probs.P_mod(i_small,:,:)))) < contour_levels(1)/100 ...
        && i_small < length(probs.gs)) || i_small == i_Tpot
    i_small = i_small + 1;
end
contour_plot(squeeze(probs.P_mod(i_small,:,:))', probs.phi, probs.gs, ...
    contour_levels, cols(2,:))
% Large grain size
i_large = length(probs.Tpot);
while (max(max(squeeze(probs.P_mod(i_large,:,:)))) < contour_levels(1)/100 ...
        && i_large > 1) || i_large == i_Tpot
    i_large = i_large - 1;
end
contour_plot(squeeze(probs.P_mod(i_large,:,:))', probs.phi, probs.gs, ...
    contour_levels, cols(3,:))
% Best fitting grain size
contour_plot(squeeze(probs.P_mod(i_Tpot,:,:))', probs.phi, probs.gs, ...
    contour_levels, cols(1,:))

axes('position',[0.7 0.15 0.225 0.1]); hold on
plot(probs.Tpot, zeros(size(probs.Tpot)), 'k-');
plot(probs.Tpot(i_small),0,'ko','markerfacecolor',cols(2,:),'markersize',ms);
plot(probs.Tpot(i_large),0,'ko','markerfacecolor',cols(3,:),'markersize',ms);
plot(probs.Tpot(i_Tpot),0,'ko','markerfacecolor',cols(1,:),'markersize',ms);
xlabel(tstr); ylim([0 eps]); daspect([1 1 1])





% % And plot the possible velocity models
% figure('color','w','position',[100 50 800 600]);
% axes('visible','off','position',[0 0 1 1]);
% text(0.47, 0.95, strrep(seismic_obs.q_method,'_',' '),...
%     'fontweight','bold','fontsize',24);
%
% axes('position',[0.1 0.1 0.4 0.8]); hold on;
% obs_vel_c = [81 209 70]./255;
% patch([seismic_obs.medianVs + seismic_obs.medianVs_error; ...
%     flipud(seismic_obs.medianVs - seismic_obs.medianVs_error)], ...
%     [seismic_obs.depth; flipud(seismic_obs.depth)], ...
%     obs_vel_c,'edgecolor','none','facealpha',0.25)
% plot(seismic_obs.medianVs, seismic_obs.depth,'color',...
%     obs_vel_c,'linewidth',2); axis ij;
% yl = get(gca,'ylim'); xl = [3.8 5.0]; xlim(xl);
% if yl(2) > 300; yl(2) = 350; end; ylim(yl);
% ylabel('Depth (km)'); xlabel('Vs (km/s)');
% set(gca,'xaxislocation','top'); box on
%
% % Plot on the Moho
% plot(xl,seismic_obs.Moho*[1 1],'--','color',0.6*[1 1 1]);
% %text(mean([xl,xl(1)]), seismic_obs.Moho - 5, 'Moho');
%
% % Plot on the LAB
% plot(xl,seismic_obs.LAB*[1 1],'b--');
% %text(mean([xl,xl(1)]), seismic_obs.LAB - 5, 'LAB');
%
%
% % Plot on all models within top 5% probability
% %inds = find(probs.P_mod > contour_levels(end));
% [~,inds]=sort(probs.P_mod(:),'descend');
% top_vals = zeros(50,3);
% for k = 1:size(top_vals,1)
%     Vave = mean(sweepBox(inds(k)).VBR.(seismic_obs.q_method).Vave,2).*1e-3;
%     plot(Vave, sweepBox(inds(k)).info.Z_km,'-','color',0.6*[1 1 1]);
%     top_vals(k,:) = [sweepBox(inds(k)).info.Tpot ...
%         max(sweepBox(inds(k)).info.phi) sweepBox(inds(k)).info.gs];
%     plot(vs_vals(inds(k))*[1 1], seismic_obs.depthrange, ':',...
%         'color',[0.6 0 0],'linewidth',1);
% end
%
% plot(mean(sweepBox(i_max).VBR.(seismic_obs.q_method).Vave,2).*1e-3, ...
%     sweepBox(i_max).info.Z_km,'k-','linewidth',2);
%
%
% patch(xl([1 2 2 1 1]), seismic_obs.depthrange([1 1 2 2 1]),...
%     'r','facealpha',0.3);
% plot(vs_vals(i_max)*[1 1], seismic_obs.depthrange, '-',...
%     'color',[0.6 0 0],'linewidth',2);
% plot(seismic_obs.asth_v*[1 1], seismic_obs.depthrange,'r--',...
%     'linewidth',2);
%
%
% axes('position',[0.6 0.2 0.35 0.6]); hold on; box on
% scatter(top_vals(:,1), top_vals(:,2), 50, log10(top_vals(:,3)),...
%     'filled','s','markerfacealpha',0.5)
% xlabel(tstr); ylabel(phistr); c=colorbar('location','southoutside');
% xlabel(c,gstr);
% xlim([min(probs.Tpot) max(probs.Tpot)]);
% ylim([min(probs.phi) max(probs.phi)]);
% caxis([min(log10(probs.gs)) max(log10(probs.gs))])

end


function contour_plot(P_mod, x, y, contour_levels, col)

nc = length(contour_levels);

% linestyles = {':','-.','--','-','-'};
% if nc < 5; linestyles = linestyles(end-nc+1:end);
% else; eval(['linestyles = {' repmat(''':'',',1,nc-5+1) ...
%         '''-.'',''--'',''-'',''-''};']);
% end



for ic = 1:nc
    [~, h] = contour(x, y, P_mod,contour_levels(ic)*[.01 .01]);
    %set(h,'linestyle',linestyles{ic},'color',col);
    set(h,'linewidth',ic,'color',col);
    if ic == nc; set(h,'fill','on','facecolor',col); end

end

end



function sweepBox = old_generate_parameter_sweep(Work, zPlate, sweep)


Work.boxname = ['Box_' Work.Box_base_name '_VBR_py.mat'];
Work.boxpath = [Work.hmdir 'Boxes/' Work.Box_base_name '/'];
load([Work.boxpath Work.boxname]);
Work.sweep_boxname = ['Box_' Work.Box_base_name '_sweep_' ...
    num2str(zPlate.zPlate) 'km.mat'];
clc
fprintf('\n\nPlate thickness: %g\n', ...
    Box(1,1).info.var2range(zPlate.zPlate_ind))

if exist([Work.boxpath Work.sweep_boxname],'file')
    fprintf('\nYou already have a saved parameter sweep\n\t%s\n', ...
        Work.sweep_boxname(1:end-4));
    if_recalculate = input('Recalculate y/[n]? ', 's');

    if ~strcmp(if_recalculate,'y')
        load([Work.boxpath Work.sweep_boxname])
        return
    end
end


% add VBR paths
Work.VBR_version = 'VBR_v0p95';
addpath([Work.fundir '4_VBR/'  Work.VBR_version ],...
    [Work.fundir '4_VBR/'  Work.VBR_version '/functions'],...
    [Work.fundir '4_VBR/'  Work.VBR_version '/params'])
addpath([Work.fundir '2_PLATES/4pt0_1d_plates/01_functions/'])


% For a given zPlate, we want to run through all combinations of Tpot, melt
% fraction, and grain size.
sweep.Tpot_range = Box(1,1).info.var1range;
sweep.Tpot_units = Box(1,1).info.var1units;
sweep.zPlate = Box(1,1).info.var2range(zPlate.zPlate_ind);
sweep.zPlate = Box(1,1).info.var2units;

n_Tpot   = numel(sweep.Tpot_range);
n_phi    = numel(sweep.phi_range);
n_gs     = numel(sweep.gs_range);

tstart = cputime;

for i_Tpot = n_Tpot:-1:1  % run backwards so structure is preallocated

    Tbox = Box(i_Tpot, zPlate.zPlate_ind);
    VBR_in = Tbox.Frames(end).VBR;

    % Extract relevant parameters from the initial run
    sweep.zLAB = Tbox.run_info.zLAB(end);
    sweep.Z_km = Tbox.run_info.Z_km;
    sweep.T_z  = Tbox.Frames(end).T(:,1);
    sweep.Tpot = sweep.Tpot_range(i_Tpot);

    for i_phi = n_phi:-1:1

        % Find melt fraction profile (melt in asthenosphere)
        Z_phiBumpCenter_km = sweep.zLAB./1e3; % LAB depth in km
        ind_phiBumpCenter = find(sweep.Z_km > Z_phiBumpCenter_km,1)-1;
        N_z_steps = floor(length(sweep.Z_km)/10);
        [step_vec]= make_meltStep(length(sweep.Z_km), N_z_steps, ...
            ind_phiBumpCenter);
        sweep.phi = step_vec.*sweep.phi_range(i_phi);

        for i_gs = n_gs:-1:1

            fprintf(['-----------------------------------------\n' ...
                'Run %g of %g\n'], n_gs*n_phi*n_Tpot + 1 - ...
                (i_gs + (i_phi-1)*n_gs + (i_Tpot-1)*n_phi*n_gs),...
                n_gs*n_phi*n_Tpot);

            sweep.gs = sweep.gs_range(i_gs);

            % Put values in the VBR and run it
            VBR.in = VBR_in.in;
            VBR.in.SV.phi = sweep.phi.*1e-2;
            VBR.in.SV.dg_um = sweep.gs;

            VBR = VBR_spine(VBR);
            anelastic_methods = fieldnames(VBR.out.anelastic);
            if isfield(VBR.out.anelastic,'YT_maxwell')
                VBR.out.anelastic.YT_maxwell.Vave = ...
                    VBR.out.anelastic.YT_maxwell.V;
            end
            for i_an = 1:length(anelastic_methods)
                vbr.(anelastic_methods{i_an}).Vave = ...
                    VBR.out.anelastic.(anelastic_methods{i_an}).Vave;
                %vbr.(anelastic_methods{i_an}).Qinv = ...
                 %   VBR.out.anelastic.(anelastic_methods{i_an}).Qinv;
            end

            sweepBox(i_Tpot, i_phi, i_gs).info = orderfields(sweep);
            sweepBox(i_Tpot, i_phi, i_gs).VBR = vbr;

        end
    end
    save([Work.boxpath Work.sweep_boxname], 'sweepBox');
end
tend = cputime;

save([Work.boxpath Work.sweep_boxname], 'sweepBox');

fprintf(['-----------------------------------------\n' ...
    'Computations complete!\n\tTotal VBR time: %.1f min\n' ...
    'Box saved to %s.'],(tend-tstart)/60, Work.sweep_boxname);

end
