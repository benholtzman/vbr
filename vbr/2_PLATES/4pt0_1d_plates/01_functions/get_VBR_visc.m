function eta = get_VBR_visc(Vark)

% %% make sure VBR is in path
%    VBR_version = 'VBR_v0p95';
%    addpath(['../../4_VBR/'  VBR_version ],...
%            ['../../4_VBR/'  VBR_version '/functions'],...
%            ['../../4_VBR/'  VBR_version '/params'])

%% write VBR methods lists (these are the things to calculate)
   visc_method='HK2003';
   VBR.in.viscous.methods_list={visc_method};

%% set relevant state variables
   VBR.in.SV.P_GPa = Vark.P./1e9 ;
   VBR.in.SV.T_K = Vark.T +273;
   VBR.in.SV.phi =  Vark.phi ;
   VBR.in.SV.dg_um = Vark.dg_um ;
   VBR.in.SV.sig_MPa = Vark.sig_MPa ;
   VBR.in.SV.Ch2o = Vark.Cs_H2O;

%% calculate viscosity
   [VBR] = VBR_spine(VBR) ;

%% return viscosity
   eta=VBR.out.viscous.(visc_method).eta_total;
   eta(eta>1e26)=1e26;

end
