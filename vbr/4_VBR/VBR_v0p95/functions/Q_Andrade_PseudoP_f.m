function [VBR]=Q_Andrade_PseudoP_f(VBR)
%disp('doing Andrade PseudoPeriod calc now')
%% ================================
%% read in variables and parameters
%% ================================

%  (BAD RENAMING OF VARIABLES !G==M, sorta... M,J general, G is shear !)
   if isfield(VBR.in.elastic,'poro_Takei')
     Mu_in = VBR.out.elastic.poro_Takei.Gu ;
     %disp('using poroTakei Gu')
   elseif isfield(VBR.in.elastic,'anharmonic')
     Mu_in = VBR.out.elastic.anharmonic.Gu ;
     %disp('using anharmonic Gu')
   end
   rho_in = VBR.in.SV.rho ;
   T_K_mat = VBR.in.SV.T_K ;
   P_Pa_mat = VBR.in.SV.P_GPa.*1e9 ; % convert pressure GPa to Pa = GPa*1e9
   d_mat = VBR.in.SV.dg_um ; % microns grain size
   phi = VBR.in.SV.phi ; % should be called _mat for consistency-- matrix w same dims as other SVs.
   f_vec = VBR.in.SV.f;  % frequency [this, however, is not a mat -- is a 1D vector]

%  Andrade parameters, set in params file
   Andrade_params=VBR.in.anelastic.AndradePsP;
   n = Andrade_params.n  ; % andrade exponent
   Tau_MR = Andrade_params.Tau_MR ; %??
   Beta = Andrade_params.Beta ;
   TR = Andrade_params.TR;% Kelvins
   PR = Andrade_params.PR *1e9 ; %2e8 ; % convert pressure GPa to Pa = GPa*1e9
   dR = Andrade_params.dR ; % 3.1 microns grain size

   E = Andrade_params.E ; % J/mol
   R = Andrade_params.R ;
   Vstar = Andrade_params.Vstar ; % m^3/mol (Activation Volume? or molar volume?)
   m = Andrade_params.m ;

%  Elastic-GBS relaxation peak parameters, set in params file
   Te = Andrade_params.Te ; % = 0.1 ; % not sure what this is
   Tgbs = Andrade_params.Tgbs ; % = 0.1 ;% sec
   Delta = Andrade_params.Delta ;% = 0.43 ; % Relaxation strength

%  EFFECT OF MELT, hypothetical:
%  sharper response than creep at phi_c, but less sensitive at higher phi (for HTB only)
   alpha = Andrade_params.melt_alpha ;
   phi_c = Andrade_params.phi_c ;
   x_phi_c = Andrade_params.x_phi_c ;


%% ========================================
%% simple calculations from input paramters
%% ========================================

   w_vec = 2*pi.*f_vec ; % angular frequency
   Tau0_vec = 1./f_vec ; % period
   Ju_in = 1./Mu_in ; % compliance

%  calculate X_tilde:
%  X_tilde in paper == scale_mat (dimensionless) (like strain rate, not viscosity)

%  truly melt free
  Xtilde = ((d_mat./dR).^-m).*exp((-E/R).*(1./T_K_mat-1/TR)) ...
                           .*exp(-(Vstar/R).*(P_Pa_mat./T_K_mat-PR/TR));
  if VBR.in.GlobalSettings.melt_enhacement==0
    x_phi_c=1;
  else
    Xtilde = Xtilde / x_phi_c ;
  end

%  melt enhancement
   [Xtilde_prime] = sr_melt_enhancement(phi,alpha,x_phi_c,phi_c) ;
   Xtilde = Xtilde_prime.*Xtilde ;

%% ===========================
%% allocation of Qstruct and V
%% ===========================

   n_freq = numel(f_vec);
   sz = size(Mu_in);

%  frequency dependent vars
   J1 = proc_add_freq_indeces(zeros(sz),n_freq);
   J2 = J1; Qa = J1; Qinv = J1; Ma = J1; Va = J1;
   J1_gbs = J1; J2_gbs = J1; Q_gbs = J1; M_gbs = J1;
   J1_comp = J1; J2_comp = J1; Q_comp = J1; M_comp = J1; Va_comp = J1;

%  vectorized rho and Vave
   n_th = numel(Ju_in); % total elements
   rho_vec = reshape(rho_in,size(Ma(1:n_th)));
   Vave=reshape(zeros(sz),size(Ma(1:n_th)));


%% =============================
%% calculate material properties
%% =============================

% scalar values
  param1 = Beta*gamma(1+n)*cos(n*pi/2);
  param2 = Beta.*gamma(1+n)*sin(n*pi/2);

% loop over frequency
  for f = 1:n_freq
%     get linear index of J1, J2, etc.
      ig1 = 1+(f - 1) * n_th; % the first linear index in current frequency
      ig2 = (ig1-1)+ n_th; % the last linear index in current frequency

%     pseudoperiod master variable
      wX_mat = 2*pi./(Tau0_vec(f).*Xtilde) ;
      w = w_vec(f);

%     andrade model
      J1(ig1:ig2) = Ju_in.*(1 + param1 * (wX_mat.^-n)) ;
      J2(ig1:ig2) = Ju_in.*(param2 * (wX_mat.^-n) + Xtilde./(Tau_MR.*w));

      Qa(ig1:ig2) = J1(ig1:ig2)./J2(ig1:ig2) ;
      Qinv(ig1:ig2) = 1./Qa(ig1:ig2);
      Ma(ig1:ig2) = (J1(ig1:ig2).^2 + J2(ig1:ig2).^2).^(-1/2) ;

%     gbs relaxation bump
      J1_gbs(ig1:ig2) = Ju_in.*Delta./(1+Tgbs.^2.*w.^2) ;
      J2_gbs(ig1:ig2) = Ju_in.*Delta.*(w.*Te)./(1+Tgbs.^2.*w.^2) ;
      Q_gbs(ig1:ig2) = J1_gbs(ig1:ig2)./J2_gbs(ig1:ig2) ;
      M_gbs(ig1:ig2) = (J1_gbs(ig1:ig2).^2 + J2_gbs(ig1:ig2).^2).^(-1/2) ;

%     composite
      J1_comp(ig1:ig2) = J1(ig1:ig2) + J1_gbs(ig1:ig2) ;
      J2_comp(ig1:ig2) = J2(ig1:ig2) + J2_gbs(ig1:ig2) ;
      Q_comp(ig1:ig2) = J1_comp(ig1:ig2)./J2_comp(ig1:ig2) ;
      M_comp(ig1:ig2) = (J1_comp(ig1:ig2).^2 + J2_comp(ig1:ig2).^2).^(-1/2);

%     velocities [m/s]
%      andrade
       Va(ig1:ig2) = sqrt(Ma(ig1:ig2)./rho_vec) ;
%      average (divided by nfreq outside loop)
       Vave = Vave + Va(ig1:ig2);
%      composite
       Va_comp(ig1:ig2) = sqrt(M_comp(ig1:ig2)./rho_vec) ;
  end

%% BOOKKEEPING
   VBR.out.anelastic.AndradePsP.Vave=reshape(Vave/n_freq,sz);
   VBR.out.anelastic.AndradePsP.J1 = J1;
   VBR.out.anelastic.AndradePsP.J2 = J2;
   VBR.out.anelastic.AndradePsP.Q = Qa;
   VBR.out.anelastic.AndradePsP.Qinv = Qinv;
   VBR.out.anelastic.AndradePsP.M=Ma;
   VBR.out.anelastic.AndradePsP.V=Va;
   VBR.out.anelastic.AndradePsP.J1_gbs = J1_gbs;
   VBR.out.anelastic.AndradePsP.J2_gbs = J2_gbs;
   VBR.out.anelastic.AndradePsP.Q_gbs = Q_gbs;
   VBR.out.anelastic.AndradePsP.M_gbs=M_gbs;
   VBR.out.anelastic.AndradePsP.J1_comp = J1_comp;
   VBR.out.anelastic.AndradePsP.J2_comp = J2_comp;
   VBR.out.anelastic.AndradePsP.Q_comp = Q_comp;
   VBR.out.anelastic.AndradePsP.M_comp=M_comp;
   VBR.out.anelastic.AndradePsP.Va_comp=Va_comp;
   %disp('VBR shoulda added And...')

end
