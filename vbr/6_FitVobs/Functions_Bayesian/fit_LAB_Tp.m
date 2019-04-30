function zPlate_fit = fit_LAB_Tp(Work, seismic_obs, set_Tp)
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %
%  Fit the measured LAB and Vs to find reasonable models       %
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %
% Note that this does do a joint fit between temperature and seismic LAB
% depth, but really we are only looking for Plate Thickness (zPlate).
% If there is melt present (i.e. seismic velocities are low), we will not
% get good constraints on the potential temperature from this inversion
% alone.  However, zPlate estimates are much more strongly dependent on
% seismic LAB depth than on temperature estimates.

Q_LAB = 800;
per_bw_max = 30;  % max period to use for fitting (s)
per_bw_min = 10;  % min period to use for fitting (s)

Work.boxname = ['Box_' Work.Box_base_name '_VBR_py.mat'];
Work.boxpath = [Work.hmdir 'Boxes/' Work.Box_base_name '/'];


%    %    %    %    %    %    %
%    READ PRIOR MODELS IN     %
%    %    %    %    %    %    %
times.t0 = datetime;
load([Work.boxpath Work.boxname]);

% Pull some vectors to make processing easier
zPlate_vec = Box(1,1).info.var2range;
Tpot_vec = Box(1,1).info.var1range;
f_band = Box(1,1).Frames(end).VBR.in.SV.f;

clc;
fprintf(['\n\nTpot: %.0f - %.0f C\nzPlate: %.0f - %.0f km\n' ...
    'Freq. band: %.3f - %.1f Hz (%.1f - %.1f s)\n'],...
    Tpot_vec(1), Tpot_vec(end), zPlate_vec(1), zPlate_vec(end),...
    f_band(1), f_band(end), 1/f_band(1), 1/f_band(end));

checkQ = checkForSelectedQMethod(Box, seismic_obs.q_method);
if ~checkQ; return; end

% Find the values and indices for the frequency range of interest
i_fmin = getnearest(f_band,1/per_bw_max); %fmin = f_band(i_fmin);
i_fmax = getnearest(f_band,1/per_bw_min); %fmax = f_band(i_fmax);
times.t1 = datetime;

%    %    %    %    %    %    %    %
%  FIND BEST MODELS FOR LAB DEPTH  %
%    %    %    %    %    %    %    %
dont_use_input_Qlab = 1;
[fits.Res_LAB, fits.ind, fits.LAB] = ...
    find_LAB_Q_Res(Box,Q_LAB, seismic_obs.LAB, ...
    i_fmin, i_fmax, seismic_obs.q_method, dont_use_input_Qlab);
fits.Rmin = min(fits.Res_LAB(:));

fprintf('\nPlate best fit - \n  min residual: %f\n',...
    fits.Rmin);

% Find average Vs best fitting
[fits.Res_Vs_adavg_mat, fits.Vs_adavg_mat] = ...
    find_Vs_adavg_Res(Box,seismic_obs.asth_v, i_fmin, i_fmax, ...
    seismic_obs.q_method);

% Find the joint best fitting model!
fits.Res_Vs_N = fits.Res_Vs_adavg_mat./...
    max(fits.Res_Vs_adavg_mat(:));
fits.Res_zLAB_N = fits.Res_LAB./...
    max(fits.Res_LAB(:));

% See manual, Menke book Ch 11
fits.P_Vs = (2*pi*fits.Res_Vs_N).^-0.5 .* ...
    exp(-0.5*fits.Res_Vs_N);
fits.P_zPlate = (2*pi*fits.Res_zLAB_N).^-0.5 .* ...
    exp(-0.5*fits.Res_zLAB_N);
fits.P_JOINT = (fits.P_Vs.*fits.P_zPlate).^2;

[~, fits.i_best] = max(fits.P_JOINT(:));
%[r, c] = ind2sub(size(Box), fits.i_best);

% Actually return the best zPlate at input potential temperature
[~, zPlate_fit.Tpot_ind] = min(abs(Tpot_vec - set_Tp));
zPlate_fit.Tpot = Tpot_vec(zPlate_fit.Tpot_ind);
[~, zPlate_fit.zPlate_ind] = max(fits.P_zPlate(zPlate_fit.Tpot_ind,:));
zPlate_fit.zPlate = zPlate_vec(zPlate_fit.zPlate_ind);
fits.ij_best = [zPlate_fit.Tpot_ind, zPlate_fit.zPlate_ind];

times.t2 = datetime;
fprintf('\n\n%.1f seconds to do the finding\n',seconds(times.t2-times.t1));


% Add in a bit here using surface wave kernels...?!


%    %    %    %    %    %    %    %
%          PLOT EVERYTHING         %
%    %    %    %    %    %    %    %
Layout = buildLayout(1);
figure('color','w','position',[50 30 1200 600]);

% Q(Z) plot
plot_Q_profiles(Layout.Qz, Box, fits, seismic_obs, i_fmin, i_fmax);
% Vs(Z) plot
plot_Vs_profiles(Layout.Vsz, Box, fits, seismic_obs, i_fmin, i_fmax);
% Residuals (can plot up to two 'temps' at a time)
layout = Layout.rh;  clr = [1 0 0];

residPlots(layout, Tpot_vec, zPlate_vec, fits, clr)
ij_best = fits.ij_best;
fprintf('\n\n    Tpot: %.0f; zPlate: %.0f\n',...
    Tpot_vec(ij_best(1)),zPlate_vec(ij_best(2)));


zPlate_fit.LAB_obs = seismic_obs.LAB;
zPlate_fit.Vs_obs = seismic_obs.asth_v;

end

% FITTING FUNCTIONS
function q_method_exists = checkForSelectedQMethod(box,q_method)

frame = box(1,1).Frames(end);
methods = fieldnames(frame.VBR.out.anelastic);

if any(strcmp(methods,q_method))
    q_method_exists = true;
else
    fprintf('\nAvailable Q methods:\n');
    for k=1:length(methods); fprintf('%s\n',methods{k}); end
    fprintf('\n\nQ method %s does not exist!\n',q_method)
    q_method_exists = false;
end

end

function Pp_f_mnstd = meanstd_fband(Pprop_fband)
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %
%  Function to get mean and standard deviation at each depth   %
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %

Pp_f_mnstd = [mean(Pprop_fband,2) ...
    std(Pprop_fband,0,2)];

end

function idx = getnearest(array,val)
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %
%  Function to get index of nearest value in array to val      %
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %

[~,idx] = min(abs(array-val));

end

function [Res_lab_Q_mat, ind_zLAB_Q_mat, Z_LAB_Q_mat] =  find_LAB_Q_Res(Box, ...
    Q_LAB, zLAB_obs_km, i_fmin, i_fmax, q_method, method)
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %
%  Function to find best fitting models (calculate residuals)  %
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %

[n_var1, n_var2] = size(Box);
Res_lab_Q_mat = zeros(size(Box));
ind_zLAB_Q_mat = zeros(size(Box));
Z_LAB_Q_mat = zeros(size(Box));

for i_var1 = 1:n_var1
    for i_var2 = 1:n_var2

        Z_km = Box(i_var1,i_var2).run_info.Z_km;
        frame = Box(i_var1,i_var2).Frames(end);

        switch q_method
            case 'AndradePsP'
                Qs_f_mat = frame.VBR.out.anelastic.(q_method).Qa;
            otherwise
                Qs_f_mat = frame.VBR.out.anelastic.(q_method).Q;
        end
        Qs_f_band = Qs_f_mat(:,i_fmin:i_fmax);
        Qs_mnstd = meanstd_fband(Qs_f_band);

        % Interpolate Qz and Z_km to higher resolution
        nn_pts = 1000;
        Z_km_interp = linspace(Z_km(1),Z_km(end),nn_pts);
        Qz_interp = interp1(Z_km,Qs_mnstd(:,1),Z_km_interp);

        % Identify the LAB depth in each model
        %   Method 0: absolute value of Q (input variable Q_LAB)
        %   Method 1: a factor above the adiabatic average

        if method == 1
            fQlab = 20;

            % Find the average Q in the adiabatic part
            z_plate = Box(i_var1,i_var2).info.var2val;
            ind_zplate = getnearest(Z_km_interp,z_plate);
            Q_adbt_vec = Qz_interp(ind_zplate:end);
            Q_adbt_mean = mean(Q_adbt_vec);

            % Calculate Q at the LAB
            Q_LAB = fQlab*Q_adbt_mean;
        end

        % Find the position of the nearest Q value to Q_LAB
        ind_Qlab = find(Qz_interp<=Q_LAB,1);
        Z_LAB_Qs_km = Z_km_interp(ind_Qlab);

        % Calculate the residual between the predicted and measured
        Res = (Z_LAB_Qs_km - zLAB_obs_km).^2 / zLAB_obs_km;
        Res_lab_Q_mat(i_var1,i_var2)  = Res;
        ind_zLAB_Q_mat(i_var1,i_var2) = ind_Qlab;
        Z_LAB_Q_mat(i_var1,i_var2)    = Z_LAB_Qs_km;

    end
end

end

function [Res_Vs_adavg_mat, Vs_adavg_mat] = find_Vs_adavg_Res(Box, ...
                                        Vs_adavg_obs, i_fmin, i_fmax, q_method)
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %
%  Function to find model that best fits average Vs in the     %
%  adiabatic part                                              %
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %

[n_var1, n_var2] = size(Box);
Res_Vs_adavg_mat = zeros(size(Box));
Vs_adavg_mat = zeros(size(Box));

for i_var1 = 1:n_var1
    for i_var2 = 1:n_var2

        Z_km = Box(i_var1,i_var2).run_info.Z_km;
        frame = Box(i_var1,i_var2).Frames(end);

        switch q_method
            case 'AndradePsP'
                Vs_f_mat = frame.VBR.out.anelastic.(q_method).Va;
            otherwise
                Vs_f_mat = frame.VBR.out.anelastic.(q_method).V;
        end
        Vs_f_mat = Vs_f_mat.*1e-3;
        Vs_f_band = Vs_f_mat(:,i_fmin:i_fmax);
        Vs_mnstd = meanstd_fband(Vs_f_band);

        % Calculate average Vs in the adiabatic region
        z_plate = Box(i_var1,i_var2).info.var2val;
        ind_zplate = getnearest(Z_km, z_plate);
        Z_ad_int = 20;
        ind_bot_asth = min(ind_zplate + Z_ad_int-1,length(Z_km));
        Vs_adbt_vec = Vs_mnstd(ind_zplate:ind_bot_asth,1);
        Vs_adavg = mean(Vs_adbt_vec);

        % Calculate the residual between predicted and measured
        Res = (Vs_adavg - Vs_adavg_obs).^2/Vs_adavg_obs;
        Res_Vs_adavg_mat(i_var1,i_var2) = Res;
        Vs_adavg_mat(i_var1,i_var2)     = Vs_adavg;

    end
end



end


% PLOTTING FUNCTIONS

function plot_Q_profiles(layout, Box, fits, seismic_obs,...
    i_fmin, i_fmax)
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %
%  Plots all of the Q profiles and the best fitting profiles   %
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %

axes('position',layout); hold on; box on;
axlims = [1 8; 0 350];

[n_var1, n_var2] = size(Box);

for i_var1 = 1:n_var1
    for i_var2 = 1:n_var2
        Z_km = Box(i_var1, i_var2).run_info.Z_km;
        frame = Box(i_var1, i_var2).Frames(end);
        switch seismic_obs.q_method
            case 'AndradePsP'
                Qs_f_mat = frame.VBR.out.anelastic.AndradePsP.Qa;
            otherwise
                Qs_f_mat = ...
                    frame.VBR.out.anelastic.(seismic_obs.q_method).Q;
        end
        Qs_f_band = Qs_f_mat(:,i_fmin:i_fmax);
        Qs_mnstd = meanstd_fband(Qs_f_band);

        p1 = plot(log10(Qs_mnstd(:,1)),Z_km,'k-','linewidth',1);
        p1.Color(4) = 0.2;
    end
end


% Add standard deviations to these lines
clr = [1 0 0 0.9];
i_best = fits.i_best;
Z_km = Box(i_best).run_info.Z_km;
frame = Box(i_best).Frames(end);
switch seismic_obs.q_method
    case 'AndradePsP'
        Qs_f_mat = frame.VBR.out.anelastic.AndradePsP.Qa;
    otherwise
        Qs_f_mat = frame.VBR.out.anelastic.(seismic_obs.q_method).Q;
end
Qs_f_band = Qs_f_mat(:,i_fmin:i_fmax);
Qs_mnstd = meanstd_fband(Qs_f_band);

plot(log10(Qs_mnstd(:,1)),Z_km,'-','color',clr,'linewidth',2);

% Plot best fitting Z_LAB from best Q profile
Z_LAB_Q_km = fits.LAB(i_best);
fprintf('\nLAB depth from Q:   %.1f km\n', Z_LAB_Q_km);
p2 = plot(axlims(1,:),Z_LAB_Q_km*[1 1],':','color',...
    clr,'linewidth',2);
p2.Color(4) = 0.5;

% Observed plate thickness line
p3 = plot(axlims(1,:), seismic_obs.LAB*[1 1],...
    '-','color',clr,'linewidth',7);
p3.Color(4) = 0.2;


xlim(axlims(1,:)); ylim(axlims(2,:));
title('Fit LAB depth with Q_L_A_B'); xlabel('Log_1_0(Q)');
ylabel('Depth (km)'); axis ij; set(gca,'XAxisLocation','Top')

end

function plot_Vs_profiles(layout, Box, fits, seismic_obs,...
    i_fmin, i_fmax)
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %
%  Plots all of the Vs profiles and the best fitting profiles  %
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %

axes('position',layout); hold on; box on;
axlims = [4.0 5.05; 0 350];

[n_var1, n_var2] = size(Box);
for i_var1 = 1:n_var1
    for i_var2 = 1:n_var2
        Z_km = Box(i_var1, i_var2).run_info.Z_km;
        frame = Box(i_var1, i_var2).Frames(end);
        switch seismic_obs.q_method
            case 'AndradePsP'
                Vs_f_mat = frame.VBR.out.anelastic.AndradePsP.Va;
            otherwise
                Vs_f_mat = ...
                    frame.VBR.out.anelastic.(seismic_obs.q_method).V;
        end
        Vs_f_band = Vs_f_mat(:,i_fmin:i_fmax);
        Vs_mnstd = meanstd_fband(Vs_f_band);

        p1 = plot(Vs_mnstd(:,1).*1e-3,Z_km,'k-','linewidth',1);
        p1.Color(4) = 0.2;
    end
end


% Add standard deviations to these lines
clr = [1 0 0 0.9];
i_best = fits.i_best;
Z_km = Box(i_best).run_info.Z_km;
frame = Box(i_best).Frames(end);
switch seismic_obs.q_method
    case 'AndradePsP'
        Vs_f_mat = frame.VBR.out.anelastic.AndradePsP.Va;
    otherwise
        Vs_f_mat = frame.VBR.out.anelastic.(seismic_obs.q_method).V;
end
Vs_f_band = Vs_f_mat(:,i_fmin:i_fmax);
Vs_mnstd = meanstd_fband(Vs_f_band);

plot(Vs_mnstd(:,1).*1e-3,Z_km,'-','color',clr,'linewidth',2);

% observed Vs lines
p2 = plot(seismic_obs.asth_v*[1 1], axlims(2,:),...
    '-','color',clr,'linewidth',5);
p2.Color(4) = 0.2;

% Calculated Vs lines
p3 = plot(fits.Vs_adavg_mat(i_best)*[1 1],...
    axlims(2,:),':','color',clr,'linewidth',2);
p3.Color(4) = 0.5;


xlim(axlims(1,:)); ylim(axlims(2,:));
title('Fit V_s in asth.'); xlabel('V_s (km/s)');
ylabel('Depth (km)'); axis ij; set(gca,'XAxisLocation','Top')


end

function residPlots(layouts, Tpot_vec, zPlate_vec, fits, clr)
  %  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %
  %  Plots the residual profiles                                 %
  %  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %

  ij_best = fits.ij_best;
  x_lims = zPlate_vec([1 end]); y_lims = Tpot_vec([1 end]);


  % Plot the residuals for ZLAB
  axes('position', layouts.c1); hold on; box on;
  imagesc(x_lims, y_lims,log10(fits.Res_LAB));
  colormap(gray)
  scatter(zPlate_vec(ij_best(2)), Tpot_vec(ij_best(1)), 10, clr,'filled')
  title(['Res: Z_L_A_B']);
  xlabel('Z_p_l_a_t_e (km)'); ylabel('T_p_o_t (\circC)');
  xlim(x_lims); ylim(y_lims); axis ij

  % Plot the residuals for Vs_adavg
  axes('position', layouts.c2); hold on; box on;
  imagesc(x_lims, y_lims, log10(fits.Res_Vs_adavg_mat));
  colormap(gray)
  scatter(zPlate_vec(ij_best(2)), Tpot_vec(ij_best(1)), 10, clr,'filled')
  title('Res: V_s');
  xlabel('Z_p_l_a_t_e (km)');
  set(gca,'YTickLabel', []);
  xlim(x_lims); ylim(y_lims); axis ij;

  % Plot the residuals for Vs_adavg HOT
  axes('position', layouts.c3); hold on; box on;
  imagesc(x_lims, y_lims, fits.P_JOINT);
  colormap(gray)
  scatter(zPlate_vec(ij_best(2)), Tpot_vec(ij_best(1)), 10, clr,'filled')
  title('Joint Prob');
  xlabel('Z_p_l_a_t_e (km)');
  set(gca,'YTickLabel', []);
  xlim(x_lims); ylim(y_lims); axis ij


end

function Layout = buildLayout(num_inpt)

% Depth plots
L1 = 0.06;
B1 = 0.1;
W1 = 0.15;
H1 = 0.7;

Layout.Qz  = [L1, B1, W1, H1];
Layout.Vsz = [L1+W1+0.05, B1, W1, H1];


% Residual plots
hdel = 0.015;
hdX = 10.5*hdel;

if num_inpt == 2
    % top row
    Layout.rh.c1 = [Layout.Vsz(1)+Layout.Vsz(3)+hdel*4, B1+0.5, 0.165,0.275];
    Layout.rh.c2 = [Layout.rh.c1(1)+hdX+hdel, Layout.rh.c1(2:end)];
    Layout.rh.c3 = [Layout.rh.c2(1)+hdX+hdel, Layout.rh.c1(2:end)];

    % bottom row
    Layout.rb.c1 = [Layout.Vsz(1)+Layout.Vsz(3)+hdel*4, B1+0.05, 0.165, 0.275];
    Layout.rb.c2 = [Layout.rb.c1(1)+hdX+hdel, Layout.rb.c1(2:end)];
    Layout.rb.c3 = [Layout.rb.c2(1)+hdX+hdel, Layout.rb.c1(2:end)];

elseif num_inpt == 1
    Layout.rh.c1 = [Layout.Vsz(1)+Layout.Vsz(3)+hdel*4, B1+mean([0.05,0.5]),...
        0.165,0.275];
    Layout.rh.c2 = [Layout.rh.c1(1)+hdX+hdel, Layout.rh.c1(2:end)];
    Layout.rh.c3 = [Layout.rh.c2(1)+hdX+hdel, Layout.rh.c1(2:end)];


end


end
