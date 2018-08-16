function eta = visc_calc_YT2016_solidus(T_K,P_Pa,T_sol_K,d_um,phi)
% calculates the viscosity from Hatsuki Yamauchi and Yasuko Takei, JGR 2016,
% "Polycrystal anelasticity at near-solidus temperatures,"

% reference values
  Tr=1200+273;
  Pr=1.5*1e9;
  etar=6.22*1e21;
  dr=d_um; % grain size independent beyond what is in etar

% constants
  H=462.5*1e3; % activation energy [J/mol]
  Vol=7.913*1e-6; % activation vol [m3/mol]
  R=8.314; % gas constant [J/mol/K]
  m=3; % grain size exponent -- but this does not matter since dr = d.

% calculate solidus-depdendent activation energy factor
  A_n=calcA_n(T_K./T_sol_K,phi);

% calculate the viscoscity
  eta=etar.*(d_um./dr).^m  .* exp(Vol/R.*(P_Pa./T_K-Pr./Tr)) .* ...
        A_n .* exp(H/R.*(1./T_K-1./Tr));

end

function A_n = calcA_n(Tn,phi)
  T_eta=0.94;
  gamma=5;
  lambda=0;
  A_n=zeros(size(Tn));

  A_n(Tn<T_eta)=1;

  msk=(Tn >= T_eta) & (Tn < 1);
  A_n(msk)=exp(-(Tn(msk)-T_eta)./(Tn(msk)-Tn(msk)*T_eta)*log(gamma));

  msk=(Tn > 1);
  A_n(msk)=exp(-lambda*phi(msk))/gamma;
end
