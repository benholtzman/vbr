function [Info] = init_values(settings)
 
% read in depth array
   z = settings.Zinfo.z_m; 

%  temperature [C]
%    halfspace cooling
     age = settings.age0*1e6*365*24*3600; 
     Tpot = settings.Tpot; 
     dTdz_ad = settings.dTdz_ad; 
     Diff0 = settings.Kc_olv / settings.Cp_olv / settings.rhos;  
     if strcmp(settings.T_init,'oceanic')
       Info.init.T = Tpot * erf( z / sqrt(age * Diff0)) + z*dTdz_ad;
     elseif strcmp(settings.T_init,'adiabatic')
       Info.init.T = Tpot + z*dTdz_ad;
     elseif strcmp(settings.T_init,'continental')
       zPlate=settings.T_init_Zlab*1e3; 
       [val,id]=min(abs(z-zPlate)); 
       Info.init.T = Tpot * z/zPlate;
       Info.init.T(id:end) = Info.init.T(id) + (z(id:end)-zPlate)*dTdz_ad;
     end
   
%  background velocity - convert to m/s      
   Info.init.Vbg=-settings.Vbg / (1e2 * 365 * 24 * 3600); 
  
%  uniform (for now?) scalar quantities
   Info.init.grain = settings.grain0*ones(size(z)); % grain size [m]
   Info.init.sig_MPa = 0.1*ones(size(z)); % deviatoric stress [MPa] 
   Info.init.Cs0 = settings.Cs0 * ones(size(z)); % H2O concentration [wt %]

%  composition-dependent properties
%     calculate the weighting function
      Zm=settings.Z_moho_km*1e3; 
      dzMoho = settings.Moho_thickness_km * 1e3; 
      Grade = 0.5 * (1 + erf( (z - Zm) / dzMoho));
%     initialize properties with crustal values
      Rho_0 = settings.rhos_crust * ones(size(z)); % density [kg/m3]
      Cp_0 = settings.Cp_crust * ones(size(z)); % specific heat [J/kg/K]  
      Kc_0 = settings.Kc_crust * ones(size(z)); % conductiviy [W/m/K] 
      Gu_0 = settings.Gu_crust * ones(size(z)); % conductiviy [W/m/K] 
%     add on the difference via weighting
      Info.init.Rho_0 = Rho_0 + (settings.rhos - settings.rhos_crust)*Grade;
      Info.init.Cp_0 = Cp_0 + (settings.Cp_olv - settings.Cp_crust)*Grade;
      Info.init.Kc_0 = Kc_0 + (settings.Kc_olv - settings.Kc_crust)*Grade;
      Info.init.Gu_0_GPa = Gu_0 + (settings.Gu_olv - settings.Gu_crust)*Grade; 
      Info.init.comp = Grade;     
      Info.init.Gdot = zeros(size(z)); 
      
%  porosity (initial_porosity uses the other initial conditions to calculate
%  an initial solidus, so do this last unless initial_porsotiy is changed)
   phi0=settings.phi_init*ones(size(z)); 
   [phi0,CH2O,CCO2]=initial_porosity(Info,phi0,settings,z); 
   Info.init.phi=phi0; 
   Info.init.Cs0=CH2O;
   Info.init.CCO2=CCO2;         
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% porosity initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi,CH2OS,CCO2S] = initial_porosity(Info,phi0,settings,z)

 phi = phi0; 
 CH2OS = settings.Cs0 * ones(size(z));  
 CCO2S = settings.Cs0_CO2 * ones(size(z)); 
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
%           CH2OS(iz) = settings.Cs0./100;
%           CCO2S(iz) = settings.Cs0_CO2./100;
          iz = iz + 1; 
      elseif Info.init.T(iz) > Tsol(iz) 
          iz = nz * 2;
      end    
  end
  phi(Info.init.T(:)<Tsol(:))=settings.phimin;
  
%   for iz = 1:nz
%       if Info.init.T(iz) > Tsol(iz)
%           CH2OS(iz) = CH2OS(iz) ./ (phi(iz) / settings.kd_H2O + 1-phi(iz));
%       end
%   end
  
end
