function [params] = Params_Elastic(method)

%% ========================================================================
%% Elastic Properties =====================================================
%% ========================================================================
  
if strcmp(method,'anharmonic')
  params.anharm_scale_mthd = 'Isaak' ; % Isaak or Cammarano ; 
  
  params.T_K_ref = 300 ;% room temp [K] (THIS WAS AT 1173!!)
  params.P_Pa_ref = 1e5;% 1 atm [Pa]
  params.Gu_0_ol = 81; % olivine reference shear modulus [GPa]
  params.Gu_0_crust = 30; % effective reference shear modulus for crust [GPa]
  
  if strcmp(params.anharm_scale_mthd,'Isaak')
%    Isaak, 1992
%    NOTE (bh, feb 17, 2015): the reference values for G, T, P should be self consistent 
%    (i.e. the T and P are where G was determined, but should not depend on the experiment!)
%    when using G_0 = 65 GPa (from Faul + Jackson) 
%	VBR.Gu_0_ref = 65e9 ; % in Pa	
	params.dG_dT =  -13.6*1e6 ; % Pa/K (so 13.6 is in MPa/K)    
	params.dG_dP = 1.8 ; % unitless ; % Pa/Pa 
  elseif strcmp(params.anharm_scale_mthd,'Cammarano')
%    Cammarano et al 2003: 
%    Cammarano has G_0 = 81 - 31*X_Fe, for X = 0.9 (Mg# 91), G_0 = 78.2
	params.dG_dT = -14.0e6 ; % Cammarano 2003) Pa/K (so 13.6 is in MPa/K)
	params.dG_dP = 1.4 ; %(Cammarano 2003)  unitless ; % Pa/Pa 
  end
    
  params.nu = 0.25 ; % poisson's ratio

end


if strcmp(method,'poro_Takei')  
  
%% parameters for poro-elastic melt effect  
   params.Melt_A  = 1.6 ; % 1:2.3 depending upon the wetting angle (see Yoshino). 
   params.Melt_Km = 30e9; % melt bulk modulus [Pa], Takei 2002, Table 2    
   
end
  
end
