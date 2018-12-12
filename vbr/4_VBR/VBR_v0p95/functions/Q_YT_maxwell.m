function[VBR]=Q_YT_maxwell(VBR)
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
  Ju_in  = 1./Mu_in ;

  % state VARIABLES
  rho_in = VBR.in.SV.rho ;
  T_K = VBR.in.SV.T_K ;
  P_Pa = VBR.in.SV.P_GPa.*1e9 ; % convert pressure GPa to Pa = GPa*1e9
  d = VBR.in.SV.dg_um ; % microns grain size
  phi = VBR.in.SV.phi ;

  % Frequency  =========================================
  f_vec = VBR.in.SV.f;  % frequency
  period_vec = 1./f_vec ;
  omega_vec = f_vec.*(2*pi) ;
  tau_vec = 1./omega_vec ;

% The scaling function:  ===============================
  visc_method=VBR.in.viscous.methods_list{1};
  eta_diff = VBR.out.viscous.(visc_method).diff.eta ; % viscosity for maxwell relaxation time
  % Tau_M_yt =
  tau_mxw_x = eta_diff./ Mu_in ; % maxwell relaxtion time
  %tau_norm = tau_vec./ tau_mxw ;

function[X_tau] = X_func(tau_norm_vec,params)
  %% the relaxation spectrum function
  Beta  = params.beta1 .* ones(size(tau_norm_vec));
  Alpha = params.Alpha_a - params.Alpha_b./(1+params.Alpha_c*(tau_norm_vec.^params.Alpha_taun));

  Beta(tau_norm_vec<1e-11)=params.beta2;
  Alpha(tau_norm_vec<1e-11)=params.alpha2;
  X_tau = Beta .* tau_norm_vec.^Alpha;
end

%% ===========================
%% allocation of new matrixes
%% ===========================

   n_freq = numel(f_vec);
   sz = size(Mu_in);
   n_th = numel(Mu_in); % total elements

%  frequency dependent vars
   J1 = proc_add_freq_indeces(zeros(sz),n_freq);
   J2 = J1; Q = J1; Qinv = J1; M1 = J1; M2 = J1; M = J1; V = J1;
   f_norm_glob=J1; tau_norm_glob=J1;

%  vectorized rho and Vave
   rho_vec = reshape(rho_in,size(Mu_in(1:n_th)));
   Vave=reshape(zeros(sz),size(Mu_in(1:n_th)));

% ====================================================
% LOOP over the DOMAIN of the state variables
% ====================================================
% use linear indexing!! will loop over n-dimensions of Ju_mat.
for x1 = 1:n_th  % loop using linear index!

%   pull out variables at current index
  tau_mxw = tau_mxw_x(x1);
  Ju = Ju_in(x1) ;
  rho = rho_in(x1) ;
  tau_norm = tau_vec ./ tau_mxw ; % vector ./ scalar

%   loop over frequency
  for i=1:n_freq
      i_glob = x1 + (i - 1) * n_th; % the linear index of the arrays with
                                   % a frequency index

      freq=f_vec(i); % current frequency
      f_norm=tau_mxw*freq; % normalized frequency
      max_tau_norm=1./(2*pi*f_norm); % maximum normalized tau

      tau_norm_f = max_tau_norm;%tau_norm(i) ;
      tau_norm_vec_local = linspace(0,tau_norm_f,100) ;
      %display(size(tau_norm_vec_local))
      %% QUESTION: in the theory, what is ln(tau) ??  enter ln(tau) into function instead of tau ?
      X_tau = X_func(tau_norm_vec_local,VBR.in.anelastic.YT_maxwell) ;

      %FINT1 = trapz(X_tau) ;  %@(taup) (X_tau, taup
      %int1 = Tau_fac.*quad(FINT1, 0, tau_norm_i);
      int1 = trapz(tau_norm_vec_local,X_tau./tau_norm_f) ;

      J1(i_glob) = Ju.*(1 + int1);
      %J2(i_glob) = Ju.*((pi/2)*X_tau(end) + 1/(2*pi*tau_mxw)); %this was wrong
      J2(i_glob) = Ju.*((pi/2)*X_tau(end) + tau_norm(i));

      % See McCarthy et al, 2011, Appendix B, Eqns B6 !
      J2_J1_frac=(1+sqrt(1+(J2(i_glob)./J1(i_glob)).^2))/2;
      Qinv(i_glob) = J2(i_glob)./J1(i_glob).*(J2_J1_frac.^-1);
      Q(i_glob) = 1./Qinv(i_glob);

      % THIS IS WRONG !! BH probably put this here !
      % (though may be approximate in some conditions!)
      % M1(i_glob) = 1./J1(i_glob) ;
      % M2(i_glob) = 1./J2(i_glob) ;
      % M1(i_glob) =
      M(i_glob) = 1./sqrt(J1(i_glob).^2+J2(i_glob).^2);%(M1(i_glob).^2 + M2(i_glob).^2).^(0.5) ;
      V(i_glob) = sqrt(1./(J1(i_glob)*rho)).*(J2_J1_frac.^(-1/2));

      f_norm_glob(i_glob)=f_norm;
      tau_norm_glob(i_glob)=tau_norm_f;
      Vave(x1) = Vave(x1) + V(i_glob); % add them all, divide by nfreq later
  % end loop over frequency
  end
% end the loop(s) over spatial dimension(s)
end

%% WRITE VBR
 VBR.out.anelastic.YT_maxwell.J1 = J1;
 VBR.out.anelastic.YT_maxwell.J2 = J2;
 %VBR.out.anelastic.YT_maxwell.M1 = M1;
 %VBR.out.anelastic.YT_maxwell.M2 = M2;
 VBR.out.anelastic.YT_maxwell.Q = Q;
 VBR.out.anelastic.YT_maxwell.Qinv = Qinv;
 VBR.out.anelastic.YT_maxwell.M=M;
 VBR.out.anelastic.YT_maxwell.V=V;
 VBR.out.anelastic.YT_maxwell.f_norm=f_norm_glob;
 VBR.out.anelastic.YT_maxwell.tau_norm=tau_norm_glob;
 VBR.out.anelastic.YT_maxwell.Vave = Vave./n_freq;

end
