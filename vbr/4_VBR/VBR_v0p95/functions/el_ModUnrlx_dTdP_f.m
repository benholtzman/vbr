%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ModUnrlx_dTdP_f
% calculates the effects of pressure and temperature on the unrelaxed
% shear modulus Gu (shear modulus?) 
% Gu in Pa 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ VBR ] = el_ModUnrlx_dTdP_f( VBR )

% read in temperature, pressure and density
  T_K = VBR.in.SV.T_K ; 
  P_Pa = (VBR.in.SV.P_GPa).*1e9 ;
  rho = VBR.in.SV.rho; 
  
% read in elastic parameters
  ela = VBR.in.elastic.anharmonic;
  nu = ela.nu ;
  dG_dT0 = ela.dG_dT ; % Pa/K
  dG_dP0 = ela.dG_dP  ; % dimensionless
  T_K_ref = ela.T_K_ref ;
  P_Pa_ref = ela.P_Pa_ref ;
  Gu_0=VBR.out.elastic.Gu_0;  
  dT = (T_K-T_K_ref); 
  dP = (P_Pa-P_Pa_ref); 
  
% calculate shear modulus at T,P of interest  
  Gu_TP = calc_Gu(Gu_0,dT,dP,dG_dT0,dG_dP0); 
  
% calculate bulk modulus 
  Ku_TP = calc_Ku(Gu_TP,nu);

% calculate velocities   
  [Vp,Vs] = el_VpVs_unrelaxed(Ku_TP,Gu_TP,rho);  
  
% store in VBR structure    
  anharmonic.Gu = Gu_TP ;   
  anharmonic.Ku = Ku_TP;
  anharmonic.Vpu = Vp;
  anharmonic.Vsu = Vs;   
  VBR.out.elastic.anharmonic = anharmonic; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates unrelaxed modulus at temperature, pressure above or below the
% reference values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Gu_TP = calc_Gu(Gu_0,dT,dP,dG_dT,dG_dP)  
  Gu_TP = Gu_0 + dT.*dG_dT + dP.*dG_dP;  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates bulk modulus from shear modulus and poisson ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ku = calc_Ku(Gu,nu)
 Ku = 2/3 * Gu .* (1+nu)./(1-2*nu); 
end
