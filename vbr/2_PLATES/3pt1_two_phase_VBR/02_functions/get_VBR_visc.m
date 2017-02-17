function eta = get_VBR_visc(Vark)

  addpath ../../4_VBR/VBR_v0p95/functions/
  addpath ../../4_VBR/VBR_v0p95/params/
  addpath ../../4_VBR/VBR_v0p95/

  VBR.in.viscous.methods_list={'LH2012'}; 
  
  VBR.in.SV.P_GPa = Vark.P./1e9 ;
  VBR.in.SV.T_K = Vark.T +273;  
  VBR.in.SV.rho = Vark.rho ;
  VBR.in.SV.phi =  Vark.phi ;
  VBR.in.SV.dg_um = Vark.dg_um ;
  VBR.in.SV.sig_MPa = Vark.sig_MPa ;
  VBR.in.SV.chi = Vark.comp;
  VBR.in.SV.Ch2o = Vark.Cs;
  
  [VBR] = VBR_spine(VBR) ;
  
  eta = VBR.out.viscous.LH2012.eta_total;
  eta(eta>1e26)=1e26;
end
