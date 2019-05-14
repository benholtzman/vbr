function Fits = fit_LAB_Tp(VbrBoxFile, seismic_obs, settings)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Fit the measured LAB and Vs to find reasonable models
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *    *    *    *
%
% Does:
%   1. a joint fit between potential temperature and seismic LAB depth
%   2. a fit for just Plate Thickness (zPlate) at fixed potential temp.
%
% If there is melt present (i.e. seismic velocities are low), Joint inversion
% will not return good constraints on the potential temperature but zPlate
% estimates are much more strongly dependent on seismic LAB depth than on
% temperature estimates.
%
% Input:
%   VbrBoxFile : filename for box with vbr results
%   seismic_obs : structure with seismic observations
%   settings : structure settings for fitting. Following field names.
%     .set_Tp : return best zPlate for this potential temp.
%     .q_method : Q method to use for fitting (must be in VBR box), defaults to
%                 first available
%     .per_bw_max : max period to use for fitting (s), default 30
%     .per_bw_min : min period to use for fitting (s), default 10
%     .dz_adi_km : depth range to average over to get astheno Q (default 40)
%     .interp_pts : number of points to use for interpolation of Z, Q_z
%                   (default 1000)
%     .q_LAB_method : method for LAB search (see find_LAB_Q), default 'Q_factor'
%     .q_LAB_value : value to use in LAB search, default 20 (i.e., Q_factor=20)
%
% Output:
%  zPlate_fit : structure with following fitting results
%     .resids : residual structure
%     .bestJoint : best fitting from joint minimization
%     .predictions : predicted observations (zLAB from Q, average astheno Vs)
%     .fixedTp : fitting results for fixed Tp
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp("\n\n\nFitting LAB and asthenosphere velocity")
  load(VbrBoxFile); % loads VBR structure
  settings=checkFittingSettings(settings,VBR); % process input args

  % predictions (from prior models):
  %   preds.zLAB_Q: LAB depth as measured by Q
  %   preds.meanVs average asthenosphere velocity
  preds =  predictObs(VBR,settings);

  % calculate residuals
  [resids,best]=jointResiduals(...
      'zLAB', preds.zLAB_Q, seismic_obs.LAB, seismic_obs.LAB_error, ...
      'Vs', preds.meanVs, seismic_obs.asth_v, seismic_obs.asth_v_error);

  % find best fitting zPlate given input Tpot, set_Tp.
  if isfield(settings,'set_Tp')
    Fits.fixed_Tp = findBest_Zplate_at_Tp(VBR,preds,resids,settings.set_Tp);
  end

  % store in final structure
  Fits.resids=resids;
  Fits.bestJoint=best;
  Fits.predictions=preds;
end

function settings = checkFittingSettings(settings,VBRbox)
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % define defaults, add to settings if field does not exist, calc freq range.
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % define defaults
  Defaults.per_bw_max=30;
  Defaults.per_bw_min=10;
  available_Q=fieldnames(VBRbox(1).out.anelastic);
  Defaults.q_method=available_Q{1};
  Defaults.dz_adi_km = 40;
  Defaults.interp_pts=1000;
  Defaults.q_LAB_method = 'Q_factor'; % method for LAB search (see find_LAB_Q)
  Defaults.q_LAB_value = 20; % value to use in LAB search

  % check input settings for fields, add defaults if they dont exist
  available_opts=fieldnames(Defaults);
  input_names=fieldnames(settings);
  for i = 1:numel(available_opts)
    this_opt=available_opts{i};
    if ~isfield(this_opt,settings)
      settings.(this_opt)=Defaults.(this_opt);
    end
  end

  % calculate frequency range
  settings.freq_range=sort([1/settings.per_bw_max 1/settings.per_bw_min]);
end

function predictions =  predictObs(Box,settings)
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  predicts Z_LAB based on Q for a Box
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % allocate zLAB, meanVs for every Box
  Z_LAB_Q = zeros(size(Box));
  meanVs= zeros(size(Box));
%   meanQs= zeros(size(Box)); % ** Currently unused **

  % Get frequency axis (same for all boxes)
  freq=Box(1).in.SV.f;

  % for interpolating Q, Vs
  nn_pts=settings.interp_pts;
  dz_adi_km=settings.dz_adi_km; % average Vs within dz_adi_km of LAB for adiabatic

  q_method=settings.q_method; % Q method for VBR
  q_LAB_method = settings.q_LAB_method; % method for LAB search
  q_LAB_value = settings.q_LAB_value; % value to use in LAB search
  
  for iBox = 1:numel(Box)
    Z_km = Box(iBox).Z_km;
    Qs_fz = Box(iBox).out.anelastic.(q_method).Q;
    Vs_fz = Box(iBox).out.anelastic.(q_method).V/1000; % m/s to km/s

    % 2D interplations of Q(z,f) and Vs(z,f) over depth and frequency
    [ Vs_zf, ~, ~ ] = interp_FreqZ(Vs_fz,freq,nn_pts,...
                                         Z_km,nn_pts);
    [ Qs_zf, freq_interp, Z_km_interp ] = interp_FreqZ(Qs_fz,freq,nn_pts,...
                                                             Z_km,nn_pts);
                                                         
    % Mask Vs and Q in frequency and depth
    f_mask = (freq_interp>=settings.freq_range(1) & freq_interp<=settings.freq_range(2));
    % index adiabatic velocity within dZ of LAB
    zplate=Box(iBox).BoxParams.var2val;
    z_mask = (Z_km_interp>=zplate & Z_km_interp<=zplate+dz_adi_km);
    Vs_zf_mask = Vs_zf(z_mask,f_mask);
    Qs_zf_mask = Qs_zf(z_mask,f_mask);
    
    % Calculate mean Vs and Q of mask
    meanVs(iBox)=mean(Vs_zf_mask(:));
%     meanQs(iBox)=mean(Qs_zf_mask(:)); % ** Currently unused **
    
    % find seismic LAB z
    Qs_z = mean(Qs_zf,2); % mean within frequency mask
    Z_LAB_Q(iBox) = find_LAB_Q(Qs_z,Z_km_interp,'method',q_LAB_method, ...
                                        'value',q_LAB_value,'z_min_km',zplate);
  end

  predictions.zLAB_Q=Z_LAB_Q;
  predictions.meanVs=meanVs;
end

function zPlate_at_Tp = findBest_Zplate_at_Tp(VBR,preds,resids,set_Tp)
  % finds best zPlate for a give Tpot

  % pull out the parameter vectors
  Tpot_vec = VBR(1).BoxParams.var1range;
  zPlate_vec = VBR(1).BoxParams.var2range;

  % find closest Tpot in Tpot_vec
  [~, zPlate_at_Tp.Tpot_ind] = min(abs(Tpot_vec - set_Tp));
  zPlate_at_Tp.Tpot = Tpot_vec(zPlate_at_Tp.Tpot_ind);

  % find best fitting zPlate at that Tpot and store indeces
  P_zPlate=resids.P_zLAB;
  [~, zPlate_at_Tp.zPlate_ind] = max(P_zPlate(zPlate_at_Tp.Tpot_ind,:));
  zPlate_at_Tp.zPlate = zPlate_vec(zPlate_at_Tp.zPlate_ind);
  zPlate_at_Tp.ij_best = [zPlate_at_Tp.Tpot_ind, zPlate_at_Tp.zPlate_ind];

end
