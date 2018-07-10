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

%  Maxwell scaling parameters, set in params file
%   YT_maxwell_params=VBR.in.anelastic.YT_maxwell;
% for the moment, define them here:

% The scaling function:  ===============================
  eta_diff = VBR.out.viscous.LH2012.diff.eta ;
  % display('the size of eta_diff:')
  % display(size(eta_diff))
  % Tau_M_yt =
  tau_mxw_x = eta_diff./ Mu_in ; % not the model reference viscosity !
  %tau_norm = tau_vec./ tau_mxw ;

  %% =============================
  %% the relaxation spectrum function
  %% =============================



function[X_tau] = X_func(tau_norm_vec)
    beta1 = 0.32 ;
    beta2 = 1853.0 ;
    alpha2 = 0.5 ;

    X_tau = zeros(1,length(tau_norm_vec)) ;
    for ff = 1:length(tau_norm_vec)

      tau_norm_f = tau_norm_vec(ff) ;
      alpha1 = 0.39 - 0.28./(1+2.6*(tau_norm_f.^0.1)) ; %

      if tau_norm_f>= 1e-11
        X_tau(ff) = beta1.*tau_norm_f.^alpha1;
      elseif tau_norm_f < 1e-11
        X_tau(ff) = beta2.*tau_norm_f.^alpha2;
      end

    end
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
   %J1_gbs = J1; J2_gbs = J1; Q_gbs = J1; M_gbs = J1;
   %J1_comp = J1; J2_comp = J1; Q_comp = J1; M_comp = J1; Va_comp = J1;

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
      tau_norm_f = tau_norm(i) ;
      tau_norm_vec_local = linspace(0,tau_norm_f,10) ;
      %display(size(tau_norm_vec_local))
      %% QUESTION: in the theory, what is ln(tau) ??  enter ln(tau) into function instead of tau ?
      X_tau = X_func(tau_norm_vec_local) ;

      %FINT1 = trapz(X_tau) ;  %@(taup) (X_tau, taup
      %int1 = Tau_fac.*quad(FINT1, 0, tau_norm_i);
      int1 = trapz(X_tau) ;

      J1(i_glob) = Ju.*(1 + int1);
      %J2(i_glob) = Ju.*((pi/2)*X_tau(end) + 1/(2*pi*tau_mxw));
      J2(i_glob) = Ju.*((pi/2)*X_tau(end) + tau_norm(i));

      % if method==0
      % %% MY METHOD--
      %     D_vec = (alf.*Tau_X_vec.^(alf-1))./(Tau_H^alf - Tau_L^alf) ;
      %
      %     int_J1 = trapz(Tau_X_vec,(D_vec./(1+w^2.*Tau_X_vec.^2))) ;
      %     J1(i_glob) = Ju.*(1+Delta.*int_J1) ;
      %
      %     int_J2 = trapz(Tau_X_vec,((Tau_X_vec.*D_vec)./(1+w^2.*Tau_X_vec.^2))) ;
      %     J2(i_glob) = Ju.*(w*Delta*int_J2 + 1/(w*Tau_M)) ;

      % elseif method==1
      %% GEOFF's METHOD-- works !
      % GA: try iterative quadrature
      % (still seems to be a problem with high freq, w*tauH > 1 say)

          %FINT1 = @(x) (x.^(alf-1))./(1+(w.*x).^2);
          %int1 = Tau_fac.*quad(FINT1, 0, tau_norm_i);


          %J1(i_glob) = Ju.*(1 + int1);
          %J2(i_glob) = Ju.*(int2 + 1./(w.*Tau_M));
      % end

      Q(i_glob) = J1(i_glob)./J2(i_glob) ;
      Qinv(i_glob) = 1./Q(i_glob) ; % J2 / J1
      M1(i_glob) = 1./J1(i_glob) ;
      M2(i_glob) = 1./J2(i_glob) ;
      M(i_glob) = (M1(i_glob).^2 + M2(i_glob).^2).^(0.5) ;
      V(i_glob) = sqrt(M(i_glob)./rho) ;

      %Vave(x1) = Vave(x1) + V(i_glob); % add them all, divide by nfreq later
  % end loop over frequency
  end
% end the loop(s) over spatial dimension(s)
end

%% WRITE VBR
 VBR.out.anelastic.YT_maxwell.J1 = J1;
 VBR.out.anelastic.YT_maxwell.J2 = J2;
 VBR.out.anelastic.YT_maxwell.M1 = M1;
 VBR.out.anelastic.YT_maxwell.M2 = M2;
 VBR.out.anelastic.YT_maxwell.Q = Q;
 VBR.out.anelastic.YT_maxwell.Qinv = Qinv;
 VBR.out.anelastic.YT_maxwell.M=M;
 VBR.out.anelastic.YT_maxwell.V=V;
 %VBR.out.anelastic.YT_maxwell.Vave = Vave./nfreq;

end
