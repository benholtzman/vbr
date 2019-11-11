function plot_tradeoffs_posterior(posterior, sweep, obs_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_tradeoffs_posterior(posterior, sweep, obs_name)
%
% Plot the posterior probability distribution across the parameter ranges
% in sweep.
%
% Parameters:
% -----------
%      posterior           (size(sweep.Box)) matrix of posterior
%                          probability for each parameter combination
%
%      sweep               structure with the following fields
%               state_names     cell of the names of the varied parameters
%               [param name]    vector of the range of values that were
%                               calculated
%               Box             output of VBR calculation
%               (also other fields recording values relevant to the
%               calculation)
%       
%       q_method            string of the method to use for calculating
%                           the anelastic effects on seismic properties
%
%       obs_name            string of the observation name, e.g. Vs, Qinv
%                           Only used to label the figure.
%
% Output:
% -------
%       Plot of the posterior probability distribution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f = figure('position', [400, 200, 1300, 700], 'color', 'w');
posterior = posterior ./ sum(posterior(:));

if isfield(sweep, 'gs')
    sweep.gs = sweep.gs * 1e-3; % put in mm rather than microns
end

% Plot up panels to show the trade-offs

plot_box(posterior, sweep, 2, 3, 1);
plot_box(posterior, sweep, 1, 3, 2);
plot_box(posterior, sweep, 1, 2, 3);





end

function plot_box(posterior, sweep, i1, i2, i3)


xpos = 0.1 + 0.3 * (i3 - 1);
ax = axes('position', [xpos, 0.55, 0.225, 0.4]);

sh = size(posterior);
sh(i3) = 1;


% Calculate limits for the color scales so all subplots will be on 
% the same scale
marg_sc = max([...
    max(sum(sum(posterior, 1), 2)), ...
    max(sum(sum(posterior, 2), 3)), ...
    max(sum(sum(posterior, 1), 3))]);
joint_sc = max([...
    max(median(sum(posterior, 1))), ...
    max(median(sum(posterior, 2))), ...
    max(median(sum(posterior, 3))), ...
    ]) * marg_sc;

% For each pair of parameters, plot the joint probability
% i.e. p(T, phi) = sum_over_g(p(T, phi | g) p(g))
p_marginal = sum(sum(posterior, i1), i2);
p_marginal_box = repmat(p_marginal, sh(1), sh(2), sh(3));
p_joint = sum(posterior .* p_marginal_box, i3);
imagesc(sweep.(sweep.state_names{i2}), sweep.(sweep.state_names{i1}), ...
    reshape(p_joint, sh(i1), sh(i2)));
xlabel(sweep.state_names{i2})
ylabel(sweep.state_names{i1});
set(ax, 'ydir', 'normal')
caxis([0, joint_sc])

ax2 = axes('position', [xpos, 0.425, 0.225, 0.05]);
plot(reshape(sweep.(sweep.state_names{i3}), 1, []), ...
    reshape(p_marginal, 1, []))
set(ax2, 'color', 'none', 'ycolor', 'none', 'box', 'off');
xlabel(sweep.state_names{i3});
ylim([0, marg_sc])

ax3 = axes('position', [xpos, 0.3, 0.225, 0.05]);
plot_P = zeros(size(sweep.P_GPa));
plot_P(sweep.z_inds) = 1;
plot(sweep.P_GPa, plot_P)
set(ax3, 'color', 'none', 'ycolor', 'none', 'box', 'off');
xlabel('Pressure (GPa)');


end