   close all; clear
  
%% add required paths  
   addpath ../functions
   addpath ../params
   addpath ../
    
%% construct state variable fields          
   z = linspace(0,300,150)*1e3; z= z';     
   to = 110 * (1e6 * 365 * 24 * 3600); % Myrs * ( Myrs_to_seconds)   
   VBR.in.SV.T_K = 1400*erf(z./(2*sqrt(1e-6.*to))) + 273 + z/1e3 * 0.4;
   VBR.in.SV.P_GPa = z * 3300 * 9.8 /1e9;
   VBR.in.SV.sig_MPa = 0.1*ones(size(z));
   VBR.in.SV.dg_um = 0.01 * 1e6 * ones(size(z)); 
   VBR.in.SV.phi = 0.001 * (z >= 80000);
   VBR.in.SV.Ch2o = zeros(size(z)); % in PPM!   
   VBR.in.SV.rho = 3300 * ones(size(z)); % [Pa]
   VBR.in.SV.chi = ones(size(z)); 
   VBR.in.SV.f = logspace(-2.2,-1.3,10);
      
%% write method list (these are the things to calculate)
   VBR.in.elastic.methods_list={'anharmonic';'poro_Takei';'Stixrude_LithgowBertelloni'};    
   VBR.in.viscous.methods_list={'HK2003'; 'LH2012'};   
   VBR.in.anelastic.methods_list={'eBurgers';'AndradePsP'};    
   
%% load in settings that you might want to overwrite (optional)
%  (each will be called internally if you don't call them here)   
%    VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity 
%    VBR.in.elastic.poro_Takei=Params_Elastic('poro_Takei'); % unrelaxed poro-elasticity    
%    VBR.in.viscous.HK2003 = Params_Viscous('HK2003'); % viscous parameters 
%    VBR.in.viscous.LH2012 = Params_Viscous('LH2012'); % viscous parameters 
%    VBR.in.anelastic.eBurgers = Params_Anelastic('eBurgers'); % anelasticity   
%    VBR.in.anelastic.AndradePsP = Params_Anelastic('AndradePsP'); % anelasticity
   
%% call VBR    
   [VBR] = VBR_spine(VBR);  
     
%% check time in each
   disp(VBR.out.computation_time)
   
