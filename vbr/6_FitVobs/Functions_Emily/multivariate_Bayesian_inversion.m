function [probs] = multivariate_Bayesian_inversion(...
    Work, seismic_obs, zPlate)


% Parameter sweep - define ranges for phi and grain size!
sweep.phi_range = (0.0:0.002:0.03).*1e2; % melt fraction
sweep.phi_units = ' %';
sweep.gs_range = 10.^(0:0.5:4); % grain size
sweep.gs_units = ' microns';

sweepBox = generate_parameter_sweep(Work, zPlate, sweep);
[vs_vals, normalised_residual] = ...
    extract_Vs(sweepBox, seismic_obs);


% Bayes' theorem:  P(A|B) = ( P(B|A)P(A) ) / P(B)
%     In this case, A is [phi, Tp, grain size];  B is observed Vs
%
% P(T, phi, gs | Vs) = __P(_Vs_|_T,_phi,_gs)_*_P(_T,_phi,_gs_)___
%                                        P(Vs)
% (the RHS there is a fraction, btw...)

% Let's calculate some probabilities!
% Hardwire all PDFs in this function
probs = calculate_probabilities(sweepBox, seismic_obs, normalised_residual);

% And do some plotting!
plot_probs(probs, seismic_obs, vs_vals, sweepBox);

end

function sweepBox = generate_parameter_sweep(Work, zPlate, sweep)


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

function [vs_vals, normalised_residual] = extract_Vs(sweepBox, seismic_obs)

vs_vals = zeros(size(sweepBox));

for k = 1:numel(vs_vals)
    
    z = sweepBox(k).info.Z_km;
    vs  = sweepBox(k).VBR.(seismic_obs.q_method).Vave;
    
    vs_vals(k) = mean(vs(z >= seismic_obs.depthrange(1) & ...
        z <= seismic_obs.depthrange(2)))./1e3;
    
end

residual = vs_vals - seismic_obs.asth_v;
normalised_residual = residual./max(abs(residual(:)));

end

function probs = calculate_probabilities(...
    sweepBox, seismic_obs, normalised_residual)

% Remember the PDF for a normal distribution
%    PDF  = 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2))


P_mod = zeros(size(sweepBox));

% P(Vs)
% This is the probability of the observed Vs being correct
% We'll assume a normal distribution, mean = observed Vs, std = obs error
sigma = seismic_obs.asth_v_error;
P_Vs = 1/sqrt(2*pi*sigma^2); % assume mean = observed Vs so (x-mu) = 0

for k = 1:numel(sweepBox)
    
    values = sweepBox(k).info;
    
    % P(T, phi, gs)  == P(T) * P(phi) * P(gs)
    % Assume that all of these are independent of one another.......
    % And have pretty broad normal distributions - can either hard wire in
    % values or guesstimate from the size of the box calculated
    
    % P(T)
    sigma = 0.5*abs(diff(values.Tpot_range([1 end])));
    mu    = mean(values.Tpot_range);
    x     = values.Tpot;
    P_T   = 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2));
    
    % P(phi)
    sigma = 0.5*abs(diff(values.phi_range([1 end])));
    mu    = mean(values.phi_range);
    x     = max(values.phi);
    P_phi = 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2));
    
    % P(grain size)
    sigma = 1e3;%0.5*abs(diff(values.gs_range([1 end])));
    mu    = 1e3; %mean(values.gs_range);
    x     = values.gs;
    P_gs = 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2));
    
    
    % P(Vs | T, phi, gs)
    % We have a residual between observed and calculated Vs
    % So... probability = 1 - residual^2 ???
    P_Vs_given_mod = 1 - normalised_residual(k).^2;
    
    
    % All together now!
    P_mod(k) = (P_Vs_given_mod * P_T * P_phi * P_gs)./(P_Vs);
    

end

probs.P_mod = P_mod./max(P_mod(:));
probs.Tpot  = sweepBox(k).info.Tpot_range;
probs.phi   = sweepBox(k).info.phi_range;
probs.gs    = sweepBox(k).info.gs_range;

end


function plot_probs(probs, seismic_obs, vs_vals, sweepBox)

% probs.P_mod = probs.P_mod(:,:,4:end);
% probs.P_mod = probs.P_mod./max(probs.P_mod(:));
% probs.gs = probs.gs(4:end);


[~, i_max] = max(probs.P_mod(:));
[i_Tpot, i_phi, i_gs] = ind2sub(size(probs.P_mod),i_max);


figure('color','w','position',[100 100 1100 600]);
tstr = 'Potential Temperature (\circC)';
phistr = 'Melt Fraction (%)';
gstr = 'Grain size (log_1_0(\mum))';

% Plot probabilities
cols = [142 255 255; ... % best fit: blue
    250 135 251; ... % smaller value: pink
    250 255 0]./255; % larger value: yellow
contour_levels = [50 70 90 95 99];



% Probability in Tpot/phi space at constant gs
axes('position',[0.1 0.35 0.225 0.55]); hold on; box on;
xlabel(tstr); ylabel(phistr);
% Find other grain sizes to plot - test extremes first
% Small grain size
i_small = 1; 
while (max(max(probs.P_mod(:,:,i_small))) < contour_levels(1)/100 ...
        && i_small < length(probs.gs)) || i_small == i_gs
    i_small = i_small + 1;
end
contour_plot(probs.P_mod(:,:,i_small)', probs.Tpot, probs.phi, ...
    contour_levels, cols(2,:))
% Large grain size
i_large = length(probs.gs);
while (max(max(probs.P_mod(:,:,i_large))) < contour_levels(1)/100 ...
        && i_large > 1) || i_large == i_gs
    i_large = i_large - 1;
end
contour_plot(probs.P_mod(:,:,i_large)', probs.Tpot, probs.phi, ...
    contour_levels, cols(3,:))
% Best fitting grain size
contour_plot(probs.P_mod(:,:,i_gs)', probs.Tpot, probs.phi, ...
    contour_levels, cols(1,:))

axes('position',[0.1 0.15 0.225 0.1]); hold on
plot(log10(probs.gs), zeros(size(probs.gs)), 'k-');
plot(log10(probs.gs(i_small)),0,'ko','markerfacecolor',cols(2,:));
plot(log10(probs.gs(i_large)),0,'ko','markerfacecolor',cols(3,:));
plot(log10(probs.gs(i_gs)),0,'ko','markerfacecolor',cols(1,:));
xlabel(gstr); ylim([0 eps]); daspect([1 1 1])



% Probability in Tpot/gs space at constant phi
axes('position',[0.4 0.35 0.225 0.55]); hold on; box on;
xlabel(tstr); ylabel(gstr);
% Find other grain sizes to plot - test extremes first
% Small grain size
i_small = 1; 
while (max(max(squeeze(probs.P_mod(:,i_small,:)))) < contour_levels(1)/100 ...
        && i_small < length(probs.gs)) || i_small == i_phi
    i_small = i_small + 1;
end
contour_plot(squeeze(probs.P_mod(:,i_small,:))', probs.Tpot, log10(probs.gs), ...
    contour_levels, cols(2,:))
% Large grain size
i_large = length(probs.phi);
while (max(max(squeeze(probs.P_mod(:,i_large,:)))) < contour_levels(1)/100 ...
        && i_large > 1) || i_large == i_phi
    i_large = i_large - 1;
end
contour_plot(squeeze(probs.P_mod(:,i_large,:))', probs.Tpot, log10(probs.gs), ...
    contour_levels, cols(3,:))
% Best fitting grain size
contour_plot(squeeze(probs.P_mod(:,i_phi,:))', probs.Tpot, log10(probs.gs), ...
    contour_levels, cols(1,:))

axes('position',[0.4 0.15 0.225 0.1]); hold on
plot(probs.phi, zeros(size(probs.phi)), 'k-');
plot(probs.phi(i_small),0,'ko','markerfacecolor',cols(2,:));
plot(probs.phi(i_large),0,'ko','markerfacecolor',cols(3,:));
plot(probs.phi(i_phi),0,'ko','markerfacecolor',cols(1,:));
xlabel(phistr); ylim([0 eps]); daspect([1 1 1])





% Probability in phi/gs space at constant Tpot
axes('position',[0.7 0.35 0.225 0.55]); hold on; box on;
xlabel(phistr); ylabel(gstr);
% Find other grain sizes to plot - test extremes first
% Small grain size
i_small = 1; 
while (max(max(squeeze(probs.P_mod(i_small,:,:)))) < contour_levels(1)/100 ...
        && i_small < length(probs.gs)) || i_small == i_Tpot
    i_small = i_small + 1;
end
contour_plot(squeeze(probs.P_mod(i_small,:,:))', probs.phi, log10(probs.gs), ...
    contour_levels, cols(2,:))
% Large grain size
i_large = length(probs.Tpot);
while (max(max(squeeze(probs.P_mod(i_large,:,:)))) < contour_levels(1)/100 ...
        && i_large > 1) || i_large == i_Tpot
    i_large = i_large - 1;
end
contour_plot(squeeze(probs.P_mod(i_large,:,:))', probs.phi, log10(probs.gs), ...
    contour_levels, cols(3,:))
% Best fitting grain size
contour_plot(squeeze(probs.P_mod(i_Tpot,:,:))', probs.phi, log10(probs.gs), ...
    contour_levels, cols(1,:))

axes('position',[0.7 0.15 0.225 0.1]); hold on
plot(probs.Tpot, zeros(size(probs.Tpot)), 'k-');
plot(probs.Tpot(i_small),0,'ko','markerfacecolor',cols(2,:));
plot(probs.Tpot(i_large),0,'ko','markerfacecolor',cols(3,:));
plot(probs.Tpot(i_Tpot),0,'ko','markerfacecolor',cols(1,:));
xlabel(tstr); ylim([0 eps]); daspect([1 1 1])





% And plot the possible velocity models
figure('color','w','position',[100 100 400 600]); hold on; 
obs_vel_c = [81 209 70]./255;
patch([seismic_obs.medianVs + seismic_obs.medianVs_error; ...
    flipud(seismic_obs.medianVs - seismic_obs.medianVs_error)], ...
    [seismic_obs.depth; flipud(seismic_obs.depth)], ...
    obs_vel_c,'edgecolor','none','facealpha',0.25)
plot(seismic_obs.medianVs, seismic_obs.depth,'color',...
    obs_vel_c,'linewidth',2); axis ij; 
yl = get(gca,'ylim'); xl = [3.5 5.0]; xlim(xl);
if yl(2) > 300; yl(2) = 350; end; ylim(yl);
ylabel('Depth (km)'); xlabel('Vs (km/s)');
set(gca,'xaxislocation','top'); box on

% Plot on the Moho
plot(xl,seismic_obs.Moho*[1 1],'--','color',0.6*[1 1 1]);
%text(mean([xl,xl(1)]), seismic_obs.Moho - 5, 'Moho'); 

% Plot on the LAB
plot(xl,seismic_obs.LAB*[1 1],'b--');
%text(mean([xl,xl(1)]), seismic_obs.LAB - 5, 'LAB'); 


% Plot on all models within top 5% probability
inds = find(probs.P_mod > 0.95);
for k = 1:length(inds)
    plot(sweepBox(k).VBR.(seismic_obs.q_method).Vave.*1e-3, ...
        sweepBox(k).info.Z_km,'-','color',0.6*[1 1 1]);
end

plot(sweepBox(i_max).VBR.(seismic_obs.q_method).Vave.*1e-3, ...
    sweepBox(i_max).info.Z_km,'k-','linewidth',2);


patch(xl([1 2 2 1 1]), seismic_obs.depthrange([1 1 2 2 1]),...
    'r','facealpha',0.3);
plot(seismic_obs.asth_v*[1 1], seismic_obs.depthrange,'r:',...
    'linewidth',2);
plot(vs_vals(i_max)*[1 1], seismic_obs.depthrange, '--',...
    'color',[0.6 0 0],'linewidth',2);


end

function contour_plot(P_mod, x, y, contour_levels, col)

alphas = [0.1 0.2 0.4 0.6 0.8];

for ic = 1:length(contour_levels)
    c = contourc(x, y, P_mod,contour_levels(ic)*[.01 .01]);
    if isempty(c); break; end
    c(:,c(1,:)==c(1,1))=[]; pc = c; 
    if c(1,end) == max(x); pc = [c [c(1,end); c(2,1)]]; end
    if c(2,end) == max(y); pc = [c [c(1,1); c(2,end)]]; end
    patch(pc(1,:),pc(2,:),col,'facealpha',alphas(ic), ...
        'edgecolor','none')
end
xlim(x([1 end])); ylim(y([1 end]));

end
