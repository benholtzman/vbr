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

function [vs_vals, normalised_residual] = extract_Vs(sweepBox, seismic_obs)

vs_vals = zeros(size(sweepBox));

for k = 1:numel(vs_vals)
    
    z = sweepBox(k).info.Z_km;
    vs  = sweepBox(k).VBR.(seismic_obs.q_method).Vave;
    
    vs_vals(k) = mean(vs(z >= seismic_obs.depthrange(1) & ...
        z <= seismic_obs.depthrange(2)))./1e3;
    
end

residual = abs(vs_vals - seismic_obs.asth_v);
normalised_residual = residual./max(residual(:));

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

normalised_residual = normalised_residual.^2;
%normalised_residual = normalised_residual - min(normalised_residual(:))+0.1;
%normalised_residual = normalised_residual./max(normalised_residual(:));

for k = 1:numel(sweepBox)
    
    values = sweepBox(k).info;
    
    % P(T, phi, gs)  == P(T) * P(phi) * P(gs)
    % Assume that all of these are independent of one another.......
    % And have pretty broad normal distributions - can either hard wire in
    % values or guesstimate from the size of the box calculated
    
    % P(T)
    sigma = 1.5*abs(diff(values.Tpot_range([1 end])));
    mu    = mean(values.Tpot_range);
    x     = values.Tpot;
    P_T   = 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2));
    
    % P(phi)
    sigma = 1.5*abs(diff(values.phi_range([1 end])));
    mu    = mean(values.phi_range);
    x     = max(values.phi);
    P_phi = 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2));
    
    % P(grain size)
    sigma = 1e4;%0.5*abs(diff(values.gs_range([1 end])));
    mu    = 1e3; %mean(values.gs_range);
    x     = values.gs;
    P_gs  = 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2));
    
    
    % P(Vs | T, phi, gs)
    % We have a residual between observed and calculated Vs
    % So... probability = 1 - residual^2 ???
    P_Vs_given_mod = 1 - normalised_residual(k);
    
    
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
%contour_levels = [50 70 90 95 99];
contour_levels = [93 95 97 99];


% Probability in Tpot/phi space at constant gs
axes('position',[0.1 0.35 0.225 0.55]); hold on; box on;
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
% Large grain size
i_large = length(probs.gs);
while (max(max(probs.P_mod(:,:,i_large))) < contour_levels(1)/100 ...
        && i_large > 1) || (i_large == i_gs && i_large > 1)
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
figure('color','w','position',[100 50 800 600]);
subplot(1,2,1); hold on; 
obs_vel_c = [81 209 70]./255;
patch([seismic_obs.medianVs + seismic_obs.medianVs_error; ...
    flipud(seismic_obs.medianVs - seismic_obs.medianVs_error)], ...
    [seismic_obs.depth; flipud(seismic_obs.depth)], ...
    obs_vel_c,'edgecolor','none','facealpha',0.25)
plot(seismic_obs.medianVs, seismic_obs.depth,'color',...
    obs_vel_c,'linewidth',2); axis ij; 
yl = get(gca,'ylim'); xl = [3.8 5.0]; xlim(xl);
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
%inds = find(probs.P_mod > contour_levels(end));
[~,inds]=sort(probs.P_mod(:),'descend');
top_vals = zeros(50,3);
for k = 1:size(top_vals,1)
    plot(sweepBox(inds(k)).VBR.(seismic_obs.q_method).Vave.*1e-3, ...
        sweepBox(inds(k)).info.Z_km,'-','color',0.6*[1 1 1]);
    top_vals(k,:) = [sweepBox(inds(k)).info.Tpot ...
        max(sweepBox(inds(k)).info.phi) sweepBox(inds(k)).info.gs];
    plot(vs_vals(inds(k))*[1 1], seismic_obs.depthrange, ':',...
        'color',[0.6 0 0],'linewidth',1);
end

plot(sweepBox(i_max).VBR.(seismic_obs.q_method).Vave.*1e-3, ...
    sweepBox(i_max).info.Z_km,'k-','linewidth',2);


patch(xl([1 2 2 1 1]), seismic_obs.depthrange([1 1 2 2 1]),...
    'r','facealpha',0.3);
plot(vs_vals(i_max)*[1 1], seismic_obs.depthrange, '-',...
    'color',[0.6 0 0],'linewidth',2);
plot(seismic_obs.asth_v*[1 1], seismic_obs.depthrange,'r--',...
    'linewidth',2);


subplot(1,2,2); hold on; box on
scatter(top_vals(:,1), top_vals(:,2), 50, log10(top_vals(:,3)),...
    'filled','s','markerfacealpha',0.5)
xlabel(tstr); ylabel(phistr); c=colorbar('location','southoutside');
xlabel(c,gstr);
xlim([min(probs.Tpot) max(probs.Tpot)]); 
ylim([min(probs.phi) max(probs.phi)]);
caxis([min(log10(probs.gs)) max(log10(probs.gs))])

end

function contour_plot(P_mod, x, y, contour_levels, col)

alphas = [0.1 0.2 0.4 0.6 0.8];

for ic = 1:length(contour_levels)
    c = contourc(x, y, P_mod,contour_levels(ic)*[.01 .01]);
    if isempty(c); continue; end
    inds = find(c(1,:) == c(1,1));
    if ~isempty(inds); c(:,inds(2:end)) = []; c(2,1) = size(c,2)-1; end
    c = contourdata(c); 
    cx = c.xdata'; cy = c.ydata';
    if ~c.isopen
        px = cx; py = cy;
    else
        
        if cx(1)>cx(end); cx = fliplr(cx); cy = fliplr(cy); end
        
        c1 = contourc(x, y, P_mod,(contour_levels(ic)-0.01)*[.01 .01]);
        inds = find(c1(1,:) == c1(1,1));
        if ~isempty(inds); c1(:,inds(2:end)) = []; c1(2,1) = size(c1,2)-1; end
        c1 = contourdata(c1);
        
        if cx([1 end]) == x([1 end])
            if mean(cy) > mean(c1.ydata)
                py = [cy max(y) max(y)];
                px = [cx cx([end 1])];
            else
                py = [cy min(y) min(y)];
                px = [cx cx([end 1])];
            end
        elseif cx(1) == cx(end) || cy(1) == cy(end)
            py = cy; px = cx;
            
        else
            if mean(cx) > mean(c1.xdata)
                if cx(1) == min(x)
                    if cy(end) == min(y)
                        py = [cy min(y) max(y) max(y)];
                        px = [cx max(x) max(x) min(x)];
                    elseif cy(end) == max(y)
                        py = [cy max(y) min(y) min(y)];
                        px = [cx max(x) max(x) min(x)];
                    end
                elseif cx(end) == max(x)
                    if cy(1) == min(y) || cy(1) == max(y)
                        py = [cy cy(1)];
                        px = [cx cx(end)];
                    end
                elseif sort(cy([1 end])) == y([1 end])
                    py = [cy cy([end 1])];
                    px = [cx max(x) max(x)];
                end
            else
                if cx(end) == max(x)
                    if cy(1) == min(y)
                        py = [cy max(y) max(y) min(y)];
                        px = [cx max(x) min(x) min(x)];
                    elseif cy(1) == max(y)
                        py = [cy min(y) min(y) max(y)];
                        px = [cx max(x) min(x) min(x)];
                    end
                elseif cx(1) == min(x)
                    if cy(end) == min(y) || cy(end) == max(y)
                        py = [cy cy(end)];
                        px = [cx cx(end)];
                    end
                elseif sort(cy([1 end])) == y([1 end])
                    py = [cy cy([end 1])];
                    px = [cx min(x) min(x)];
                end
            end
        end 
    end
    
    fill(px,py,col,'facealpha',alphas(ic), ...
        'edgecolor','none')
end
xlim(x([1 end])); ylim(y([1 end]));

end


function s = contourdata(c)
% version 1.0.0.0 (1.88 KB) by Duane Hanselman
% Matlab File Exchange 2018.12.03

%CONTOURDATA Extract Contour Data from Contour Matrix C.
% CONTOUR, CONTOURF, CONTOUR3, and CONTOURC all produce a contour matrix
% C that is traditionally used by CLABEL for creating contour labels.
%
% S = CONTOURDATA(C) extracts the (x,y) data pairs describing each contour
% line and other data from the contour matrix C. The vector array structure
% S returned has the following fields:
%
% S(k).level contains the contour level height of the k-th line.
% S(k).numel contains the number of points describing the k-th line.
% S(k).isopen is True if the k-th contour is open and False if it is closed.
% S(k).xdata contains the x-axis data for the k-th line as a column vector.
% S(k).ydata contains the y-axis data for the k-th line as a column vector.
%
% For example: PLOT(S(k).xdata,S(k).ydata)) plots just the k-th contour.
%
% See also CONTOUR, CONTOURF, CONTOUR3, CONTOURC.
% From the help text of CONTOURC:
%   The contour matrix C is a two row matrix of contour lines. Each
%   contiguous drawing segment contains the value of the contour,
%   the number of (x,y) drawing pairs, and the pairs themselves.
%   The segments are appended end-to-end as
%
%       C = [level1 x1 x2 x3 ... level2 x2 x2 x3 ...;
%            pairs1 y1 y2 y3 ... pairs2 y2 y2 y3 ...]
% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2007-05-22
if nargin<1 || ~isfloat(c) || size(c,1)~=2 || size(c,2)<4
    error('CONTOURDATA:rhs',...
        'Input Must be the 2-by-N Contour Matrix C.')
end
tol=1e-12;
k=1;     % contour line number
col=1;   % index of column containing contour level and number of points
while col<size(c,2) % while less than total columns in c
    s(k).level = c(1,col); %#ok
    s(k).numel = c(2,col); %#ok
    idx=col+1:col+c(2,col);
    s(k).xdata = c(1,idx).'; %#ok
    s(k).ydata = c(2,idx).'; %#ok
    s(k).isopen = abs(diff(c(1,idx([1 end]))))>tol || ...
        abs(diff(c(2,idx([1 end]))))>tol; %#ok
    k=k+1;
    col=col+c(2,col)+1;
end
end