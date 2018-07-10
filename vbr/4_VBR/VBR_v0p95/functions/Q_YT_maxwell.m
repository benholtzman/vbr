function[VBR]=Q_Andrade_PseudoP_f(VBR)
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

  % state VARIABLES
  rho_in = VBR.in.SV.rho ;
  T_K = VBR.in.SV.T_K ;
  P_Pa = VBR.in.SV.P_GPa.*1e9 ; % convert pressure GPa to Pa = GPa*1e9
  d = VBR.in.SV.dg_um ; % microns grain size
  phi = VBR.in.SV.phi ;

  % Frequency  =========================================
  f_vec = VBR.in.SV.f;  % frequency
  omega_vec = f_vec.*(2*pi) ;
  tau_vec = 1./omega_vec
  period_vec = 1./f_vec
%  Maxwell scaling parameters, set in params file
%   YT_maxwell_params=VBR.in.anelastic.YT_maxwell;
% for the moment, define them here:

% The scaling function:  ===============================
  eta_diff_yt = VBR.out.viscous.LH2012.diff.eta ;
  % Tau_M_yt =
  tau_Mxw_flolaw = eta_diff_yt./ Mu_in ; % not the model reference viscosity !

  tau_n = tau_vec./ tau_Mxw_flolaw ;


% The fitting function: ==============================
% The relaxation spectrum X(tau)
  Beta1 = 0.32 ;
  alpha1 = 0.39 - 0.28./(1+2.6*(tau_n.^0.1)) ;

  Beta2 = 1853.0
  alpha2 = 0.5

% ====================================================
% integration to get J1 and J2:
% ====================================================

%% ===========================
%% allocation of Qstruct and V
%% ===========================

   n_freq = numel(f_vec);
   sz = size(Mu_in);

%  frequency dependent vars
   J1 = proc_add_freq_indeces(zeros(sz),n_freq);
   J2 = J1; Qa = J1; Qinv = J1; Ma = J1; Va = J1;
   %J1_gbs = J1; J2_gbs = J1; Q_gbs = J1; M_gbs = J1;
   %J1_comp = J1; J2_comp = J1; Q_comp = J1; M_comp = J1; Va_comp = J1;

%  vectorized rho and Vave
   n_th = numel(Ju_in); % total elements
   rho_vec = reshape(rho_in,size(Ma(1:n_th)));
   Vave=reshape(zeros(sz),size(Ma(1:n_th)));


   %% =============================
   %% calculate material properties
   %% =============================

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
