%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates strain rates and viscosities for input state variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VBR= sr_visc_calc(VBR)

%% extract state variables and parameters
   T_K = VBR.in.SV.T_K ; % [K]
   P_MPa = 1e3.*(VBR.in.SV.P_GPa) ; % [GPa] to [MPa]
   sig = VBR.in.SV.sig_MPa; % deviatoric stress [MPa]
   d = VBR.in.SV.dg_um ; % [um]
   phi = VBR.in.SV.phi ;
   Ch2o=VBR.in.SV.Ch2o; % in [PPM]!
   ch2o=VBR.params_viscous.ch2o_o; % effectively dry below this [PPM]
   params=VBR.params_viscous;
   
%% calculate water fugacity (store in VBR structure)    
   fH2O=water_fugacity(Ch2o,ch2o,P_MPa,T_K,params); % [MPa]
   VBR.in.SV.Fh2o=fH2O;
  
%% calculate strain rate [1/s]
   [sr,sr_tot] = flolaw_olv_3fl_Opts_f(T_K,P_MPa,sig,d,phi,fH2O,params);
  
%% store strain rates, calculate viscosity  
   VBR.out.viscous.sr_tot=sr_tot; % total strain rate 
   VBR.out.viscous.sr_diff=sr(1).sr; % diffusion (coble) creep
   VBR.out.viscous.sr_disl=sr(2).sr; % dislocation creep
      
   VBR.out.viscous.eta_total = sig*1e6./sr_tot ; % total viscosity
   VBR.out.viscous.eta_diff = sig*1e6./VBR.out.viscous.sr_diff; % diffusion creep visc
   VBR.out.viscous.eta_disl = sig*1e6./VBR.out.viscous.sr_disl; % dislocation creep visc
   
   if params.n_mech == 3; % check if we're gbs'ing
       VBR.out.viscous.sr_gbs=sr(3).sr;
       VBR.out.viscous.eta_gbs = sig*1e6./VBR.out.viscous.sr_gbs;
   end
  
end