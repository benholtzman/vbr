function [P_mod] = multivariate_Bayesian_inversion(...
    Work, seismic_obs, zPlate)


% Bayes' theorem:  P(A|B) = ( P(B|A)P(A) ) / P(B)
% In this case, A is [phi, Tp, grain size];  B is observed Vs


% Parameter sweep - define ranges for phi and grain size!
sweep.phi_range = (0.0:0.002:0.03).*1e2; % melt fraction
sweep.phi_units = ' %';
sweep.gs_range = 10.^(0:0.5:4); % grain size
sweep.gs_units = ' microns';

sweepBox = generate_parameter_sweep(Work, zPlate, sweep);
vs_vals = extract_Vs(sweepBox, seismic_obs);
residual = vs_vals - seismic_obs.asth_v;
normalised_residual = residual./max(residual(:));


% P(T, phi, gs | Vs) = __P(_Vs_|_T,_phi,_gs)_*_P(_T,_phi,_gs_)___
%                                        P(Vs)
% (the RHS there is a fraction, btw...)

% Let's calculate some probabilities!
% Remember the PDF for a normal distribution
%    PDF  = 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2))

P_mod = zeros(size(sweepBox));

for k = 1:numel(sweepBox)
    
    values = sweepBox(k).info;
    
    % P(Vs)
    % This is the probability of the observed Vs being correct
    % We'll assume a normal distribution, mean = observed Vs, std = obs error
    sigma = seismic_obs.asth_v_error;
    P_Vs = 1/sqrt(2*pi*sigma^2); % assume mean = observed Vs so (x-mu) = 0
    
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
    sigma = 0.5*abs(diff(values.gs_range([1 end])));
    mu    = mean(values.gs_range);
    x     = values.gs;
    P_gs = 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2));
    
    
    % P(Vs | T, phi, gs)
    % We have a residual between observed and calculated Vs
    % So... probability = 1 - residual^2 ???
    P_Vs_given_mod = 1 - normalised_residual(k).^2;
    
    
    % All together now!
    P_mod(k) = (P_Vs_given_mod * P_T * P_phi * P_gs)./(P_Vs);
    

end




end

function sweepBox = generate_parameter_sweep(Work, zPlate, sweep)


Work.boxname = ['Box_' Work.Box_base_name '_VBR_py.mat'];
Work.boxpath = [Work.hmdir 'Boxes/' Work.Box_base_name '/'];
load([Work.boxpath Work.boxname]);
Work.sweep_boxname = ['Box_' Work.Box_base_name '_sweep_' ...
    num2str(zPlate.zPlate) 'km.mat'];

fprintf('\n\nPlate thickness: %g\n', ...
    Box(1,1).info.var2range(zPlate.zPlate_ind))


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
                'Run %g of %g\n'], n_gs*n_phi*n_Tpot -i_gs + ...
                (i_phi-1)*n_gs + (i_Tpot-1)*n_phi*n_gs, n_gs*n_phi*n_Tpot);
    
            sweep.gs = sweep.gs_range(i_gs);

            % Put values in the VBR and run it
            VBR.in = VBR_in.in;
            VBR.in.SV.phi = sweep.phi;
            VBR.in.SV.dg_um = sweep.gs;
            
            VBR = VBR_spine(VBR);
            

            sweepBox(i_Tpot, i_phi, i_gs).info = orderfields(sweep);
            sweepBox(i_Tpot, i_phi, i_gs).VBR = VBR;
        end
    end
end
tend = cputime;

save([Work.boxpath Work.sweep_boxname], sweepBox);

fprintf(['-----------------------------------------\n' ...
    'Computations complete!\n\tTotal VBR time: %.1f min\n' ...
    'Box saved to %s.'],(tend-tstart)/60, Work.sweep_boxname);

end

function vs_vals = extract_Vs(sweepBox, seismic_obs)

vs_vals = zeros(size(sweepBox));

for k = 1:numel(vs_vals)
    
    z = sweepBox(k).info.Z_km;
    vs  = sweepBox(k).VBR.out.anelastic.(seismic_obs.q_method).Vave;
    
    vs_vals(k) = mean(vs(z >= seismic_obs.depthrange(1) && ...
        z <= seismic_obs.depthrange(2)));
    
end

end