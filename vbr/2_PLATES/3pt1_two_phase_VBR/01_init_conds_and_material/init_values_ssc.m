function [Info] = init_values_ssc(settings,Vars,Set0,Info0)
 
% read in depth arrays
   z = settings.Zinfo.z_m; %  new one
   zold = Set0.Zinfo.z_m; %   old one
% pull out last solution and relevant init conds from before. 
% interpolate to new mesh
   Info.init.T = interp1(zold,Vars.T(:,end),z);
   Info.init.Gdot = interp1(zold,Vars.Gdot(:,end),z);   
   Info.init.grain = interp1(zold,Info0.init.grain,z);
   Info.init.sig_MPa = interp1(zold,Info0.init.sig_MPa,z);
   Info.init.Rho_0 =interp1(zold,Info0.init.Rho_0,z);
   Info.init.Cp_0 =interp1(zold,Info0.init.Cp_0,z);
   Info.init.Kc_0 =interp1(zold,Info0.init.Kc_0,z);
   Info.init.Gu_0_GPa =interp1(zold,Info0.init.Gu_0_GPa,z);
   Info.init.comp =interp1(zold,Info0.init.comp,z);
   CH2OS = interp1(zold,Info0.init.Cs0,z);
   CCO2S= interp1(zold,Info0.init.CCO2,z);
   
% reset melt-specific
%  background velocity - convert to m/s      
   Info.init.Vbg=-settings.Vbg / (1e2 * 365 * 24 * 3600);     
      
%  porosity (initial_porosity uses the other initial conditions to calculate
%  an initial solidus, so do this last unless initial_porsotiy is changed)
   phi0=settings.phi_init*ones(size(z)); 
   [phi0,Gdot0]=initial_porosity(Info,phi0,CH2OS,CCO2S,settings,z); 
   Info.init.phi=phi0; 
   Info.init.Cs0=CH2OS;
   Info.init.CCO2=CCO2S; 
   Info.init.Gdot = Info.init.Gdot .* Gdot0; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% porosity initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi,Gdot0] = initial_porosity(Info,phi0,CH2OS,CCO2S,settings,z)

 phi = phi0; 
%  CH2OS = settings.Cs0 * ones(size(z));  
%  CCO2S = settings.Cs0_CO2 * ones(size(z)); 
% calculate Rho, Cp and Kc at T,P (in order to get P)             
  [rho,cp,Kc,P_hyd] = ...
     MaterialProperties(Info.init.Rho_0,Info.init.Kc_0,Info.init.Cp_0,...
                       Info.init.T+273,z,settings.P0,settings.dTdz_ad,...
                        settings.Flags.PropType);  
% calculate solidus           
  CH2O = CH2OS./settings.kd_H2O;
  CCO2 = CCO2S./settings.kd_CO2;
  Solidus = SoLiquidus(P_hyd,CH2O,CCO2,'katz');
  Tsol=Solidus.Tsol; 
  
  iz = 1; nz = numel(phi); 
  while iz < nz; 
      if Info.init.T(iz)<Tsol(iz)
          phi(iz) = settings.phimin;          
          iz = iz + 1; 
      elseif Info.init.T(iz) > Tsol(iz) 
          iz = nz * 2;
      end    
  end
  phi(Info.init.T(:)<Tsol(:))=settings.phimin;
%   Gdot0 = Info.init.T(:) > Tsol(:); 
  Gdot0 = 1; 
end
