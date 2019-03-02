function predictions =  predictObs(Box,q_method,freq_range)
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %
%  predicts Z_LAB base on Q for a Box
%
%  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   %

  % allocate zLAB, meanVs for every Box
  Z_LAB_Q = zeros(size(Box));
  meanVs= zeros(size(Box));

  % build frequency mask (same for all boxes)
  freq=Box(1).in.SV.f;
  f_mask=(freq>=freq_range(1)&freq<=freq_range(2));

  % for interpolating Q, Vs
  nn_pts=1000;
  dz_adi_km=40; % average Vs within dz_adi_km of LAB for adiabatic

  for iBox = 1:numel(Box)
    % get Q(z), V(z) (average over frequency range)
    Qs_fz = Box(iBox).out.anelastic.(q_method).Q(:,f_mask);
    Vs_fz = Box(iBox).out.anelastic.(q_method).V(:,f_mask)/1000; % m/s to km/s
    Qs_z=mean(Qs_fz,2);
    Vs_z=mean(Vs_fz,2);

    % Interpolate Qs_z, Vs_z to higher resolution on Z_km
    Z_km = Box(iBox).Z_km;
    Z_km_interp = linspace(Z_km(1),Z_km(end),nn_pts);
    Qs_z = interp1(Z_km,Qs_z,Z_km_interp);
    Vs_z = interp1(Z_km,Vs_z,Z_km_interp);

    % find seismic LAB z
    zplate=Box(iBox).BoxParams.var2val;
    Z_LAB_Q(iBox) = find_LAB_Q(Qs_z,Z_km_interp,'method','Q_factor','value',20,'z_min_km',zplate);

    % find average adiabatic velocity within dZ of LAB
    z_mask=(Z_km_interp>=zplate & Z_km_interp<=zplate+dz_adi_km);
    meanVs(iBox)=mean(Vs_z(z_mask));
  end

  predictions.zLAB_Q=Z_LAB_Q;
  predictions.meanVs=meanVs;
end
