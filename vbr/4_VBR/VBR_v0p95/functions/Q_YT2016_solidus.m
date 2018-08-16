function [VBR] = Q_YT2016_solidus(VBR)
%% Hatsuki Yamauchi and Yasuko Takei, JGR 2016, "Polycrystal anelasticity at
%% near-solidus temperatures,"
%% Eqs. 13,14,15
%%
%% Requires the solidus, VBR.in.SV.Tsolidus_K, as an additional state variable

  if isfield(VBR.in.SV,'Tsolidus_K')
    has_solidus=1;
  else
    has_solidus=0;
    disp('To use Q_YT2016_solidus, you must provide VBR.in.SV.Tsolidus_K')
  end

  if has_solidus
    % ================================
    % read in variables and parameters
    % ================================

    % elastic modulus
      if isfield(VBR.in.elastic,'poro_Takei')
       Gu_in = VBR.out.elastic.poro_Takei.Gu;
      elseif isfield(VBR.in.elastic,'anharmonic')
       Gu_in = VBR.out.elastic.anharmonic.Gu;
      end
      Ju_in  = 1./Gu_in ;

    % state VARIABLES
      rho = VBR.in.SV.rho ;
      T_K = VBR.in.SV.T_K ;
      P_Pa = VBR.in.SV.P_GPa.*1e9 ; % convert pressure GPa to Pa = GPa*1e9
      d = VBR.in.SV.dg_um ; % microns grain size
      phi = VBR.in.SV.phi ;
      Tsol_K=VBR.in.SV.Tsolidus_K;
      Tn=T_K./Tsol_K ; % solidus-normalized temperature

    % ==================================
    % frequency independent calculations
    % ==================================

    % viscoscity & maxwell time
      [eta]=visc_calc_YT2016_solidus(T_K,P_Pa,Tsol_K,d,phi);
      tau_m=eta./Gu_in;

    % calculate the Tn-dependent coefficients, A_p and sig_p
      [A_p,sig_p]=calcApSigp(Tn,phi);

    % set other constants
      alpha_B=0.38;
      A_B=0.664;
      tau_pp=6*1e-5;

    % ===========================
    % set up frequency dependence
    % ===========================
      period_vec = 1./VBR.in.SV.f ;
      pifac=sqrt(2*pi)/2;

    % frequency dependent vars
      n_freq = numel(VBR.in.SV.f);
      sz = size(Gu_in);
      J1 = proc_add_freq_indeces(zeros(sz),n_freq);
      J2 = J1; Q = J1; Qinv = J1; M1 = J1; M2 = J1; M = J1; V = J1;

    % loop over points & frequencies
      n_SVs=numel(T_K); % total elements in state variables
      for i_SV=1:n_SVs % the state variable index
        for i=1:n_freq % the frequency index
          i_glob=i_SV + (i - 1) * n_SVs; % the global linear index in J1, J2, Q, etc.

          p_p=period_vec(i)./(2*pi*tau_m(i_SV));
          ABppa=A_B*(p_p^alpha_B);
          lntaupp=log(tau_pp./p_p);
          Ju0=Ju_in(i_SV);
          J1_i=Ju0 .* (1 + ABppa/alpha_B+ ...
               pifac*A_p(i_SV).*sig_p(i_SV).*(1-erf(lntaupp./(sqrt(2).*sig_p(i_SV)))));
          J2_i=Ju0*pi/2.*(ABppa+A_p(i_SV).*(exp(-(lntaupp^2)./(2*sig_p(i_SV).^2)))) + ...
               Ju0*p_p;

          % store in global
          J1(i_glob)=J1_i;
          J2(i_glob)=J2_i;
          J2_J1_frac=(1+sqrt(1+(J2_i./J1_i).^2))/2;
          Qinv(i_glob) = J2_i./J1_i.*(J2_J1_frac.^-1);
          Q(i_glob) = 1./Qinv(i_glob);
          M1(i_glob) = 1./J1_i ;
          M2(i_glob) = 1./J2_i ;
          M(i_glob) = 1./sqrt(J1_i.^2+J2_i.^2);%(M1(i_glob).^2 + M2(i_glob).^2).^(0.5) ;
          V(i_glob) = sqrt(1./(J1_i*rho(i_SV))).*(J2_J1_frac.^(-1/2));
        end
      end

    % ===========================
    % WRITE VBR
    % ===========================
      VBR.out.anelastic.YT2016_solidus.J1 = J1;
      VBR.out.anelastic.YT2016_solidus.J2 = J2;
      VBR.out.anelastic.YT2016_solidus.M1 = M1;
      VBR.out.anelastic.YT2016_solidus.M2 = M2;
      VBR.out.anelastic.YT2016_solidus.Q = Q;
      VBR.out.anelastic.YT2016_solidus.Qinv = Qinv;
      VBR.out.anelastic.YT2016_solidus.M=M;
      VBR.out.anelastic.YT2016_solidus.V=V;

  end
end

function [A_p,sig_p] = calcApSigp(Tn,phi);
  % Tn-dependent coefficients, A_p and sig_p
    Beta=0; % ?
    A_p=zeros(size(Tn));
    A_p(Tn < 0.91)=0.01 ;
    msk=(Tn >= 0.91) & (Tn < 0.96);
    A_p(msk)= (0.01 +0.4*(Tn(msk)-0.91));
    A_p((Tn(msk) >= 0.96) & (Tn(msk) < 1))= 0.03;
    msk=(Tn>=1);
    A_p(msk)= 0.03+Beta*phi(msk);

    sig_p=zeros(size(Tn));
    sig_p(Tn < 0.92)=4;
    msk=(Tn >= 0.92) & (Tn < 1.0);
    sig_p(msk)=4 + 37.5.*(Tn(msk)-0.92);
    sig_p(Tn >= 1.0 )=7;
end
