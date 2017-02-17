function eta = get_VBR_visc(Vark)

  addpath ../../4_VBR/VBR_v0p93/functions/

  T_K = Vark.T+273 ;
  P_MPa = Vark.P/1e6 ;
  sig_MPa = Vark.sig_MPa ;
  d = Vark.dg_um ;
  phi = Vark.phi ;
  ch2o=Vark.Cs; 
   
  
  % set FlowLawOpts   
  FlowLawOpts.ch2o_o = 50; % reference water content [ppm] ("dry" below this value)
  FlowLawOpts.n_mech=3;
  FlowLawOpts.P_dep_calc='no'; 
  FlowLawOpts.StudyNameDry='LH12'; % 'HK03' or 'LH12' 
          % !!! HK03 will fail on gbs right now !!! keep at LH12 for now!!!
  FlowLawOpts.StudyNameWet='HK03'; % 'HK03' only
  FlowLawOpts.phi_c = 0.0001;
  FlowLawOpts.speedfac = 5;
  
% calculate water fugacity  
  fH2O=water_fugacity(ch2o,FlowLawOpts.ch2o_o,P_MPa,T_K,FlowLawOpts); % [MPa]

% calculate strain rate
  [sr,sr_tot] = flolaw_olv_3fl_Opts_f(T_K,P_MPa,sig_MPa,d,phi,fH2O,FlowLawOpts);

  eta = (Vark.sig_MPa)*1e6./sr_tot ;
  eta(eta>1e26)=1e26;
end
