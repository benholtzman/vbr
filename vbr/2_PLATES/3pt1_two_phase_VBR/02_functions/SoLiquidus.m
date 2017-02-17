function [Solidus] = SoLiquidus(P,H2O,CO2,solfit)   
%  get solidus info
   P_GPa = P * 1e-9;    
   
   [dTH2O,dTdH2O] = depression_katz(H2O,P_GPa);  % H2O in melt phase [wt%]
   dTCO2 = depression_dasgupta(CO2); % CO2 in melt phase [wt%]      
   
   if strcmp(solfit,'katz')
       Solidus = solidus_katz(P_GPa);      
       Solidus.Tsol = Solidus.Tsol_dry - dTH2O -dTCO2;
       Solidus.Tliq = Solidus.Tliq_dry - dTH2O;%-dTCO2;
       Solidus.Tlherz = Solidus.Tlherz_dry - dTH2O - dTCO2;       
       Solidus.dTdPsol = Solidus.dTdPsol.*(1 - 0 * (CO2>0)) / 1e9; % artificial reduction in productivity with CO2?
       Solidus.dTdPlherz=Solidus.dTdPlherz/1e9;
       Solidus.dTdPliq = Solidus.dTdPliq /1e9; % true liquidus
       Solidus.dTdH2O = -dTdH2O; % [C / wt%]
   elseif strcmp(solfit,'hirschmann')
       Solidus = solidus_hirschmann(P_GPa);
       Solidus.Tsol = Solidus.Tsol_dry - dTH2O -dTCO2;
   end

   
end
      
%__________________________________________________________________________
%
%                     [Sols]=solidus_katz(P,T)
%__________________________________________________________________________
% 
% Dry peridotie solidus from Katz et al, A new parametrization of hydrous
% mantle melting, G3, 2003, DOI: 10.1029/2002GC000433
%
% input
%    T     	temperature [C]
%    P    	pressure [GPa]
%
% output
%    Sols.           structure with solidi temperatures
%            .Tsol_dry       dry solidus [C]
%            .Tliq_dry        dry liquidus [C]
%            .Tlherz_dry   lherzolite solidus [C]
%            .dTdPliq
%            .dTdPsol
%            .dTdPlherz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sols] = solidus_katz(P)
% Parameterization Constants
    % for solidus
      A1 = 1085.7; % [C]
      A2 = 132.9; % [C/GPa]
      A3 = -5.1; % [C/GPa^2]
    % for lherzolite liquidus
      B1 = 1475; % [C]
      B2 = 80; % [C/GPa]
      B3 = -3.2; % [C/GPa^2] 
    % for true liquidus
      C1 = 1780; % [C]
      C2 = 45; % [C/GPa]
      C3 = -2.0; % [C/GPa^2]  

% Calculate Solidii
      Sols.Tsol_dry = A1+A2*P+A3*P.^2; % solidus
      Sols.Tlherz_dry=B1+B2*P+B3*P.^2; % lherzolite liquidus
      Sols.Tliq_dry = C1+C2*P+C3*P.^2; % true liquidus

% Calculate Solidii Derivatives w.r.t. P
      Sols.dTdPsol = A2+2*A3*P; % solidus
      Sols.dTdPlherz=B2+2*B3*P; % lherzolite liquidus
      Sols.dTdPliq = C2+2*C3*P; % true liquidus
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%

%__________________________________________________________________________
%
%                     [Sols]=solidus_hirschmann(P)
%__________________________________________________________________________
% 
% Mantle solidus: Experimental constraints and the effects of peridotite
%  composition, G3, 2000 
%
% input
%    T     	temperature [C]
%    P    	pressure [GPa]
%
% output
%    Sols.           structure with solidi temperatures
%            .Tsol_dry       dry solidus [C]
%            .Tliq_dry        dry liquidus [C]
%            .Tlherz_dry   lherzolite solidus [C]
%            .dTdPliq
%            .dTdPsol
%            .dTdPlherz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sols] = solidus_hirschmann(P)
% Parameterization Constants
    % for solidus
      A1 = 1108.08; % [C]
      A2 = 139.44; % [C/GPa]
      A3 = -5.904; % [C/GPa^2]
    
% Calculate Solidii
      Sols.Tsol_dry = A1+A2*P+A3*P.^2; % solidus     
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%



%__________________________________________________________________________
%
%                   dT = depression_katz()
%__________________________________________________________________________
%
% depression of peridotite solidus due to water from  Katz et al, A new 
% parametrization of hydrous mantle melting, G3, 2003, 
% DOI: 10.1029/2002GC000433
%
% P in GPa
% H2O IN MELT (Katz sets bulk H2O, calculates H2O in melt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dT,dTdH2O] = depression_katz(H2O,P)
%    constants     
     gamma = 0.75; % temperature depression exponent
     K = 43; % [C/wt%/gamma]
%    H2O Saturation
     H2Osat = 12*P.^0.6 + P;
     H2O = H2O.*(H2O<=H2Osat)+H2Osat.*(H2O>H2Osat);
%    calculate freezing point depression                  
     dT = K * (H2O).^gamma;  
     dTdH2O = gamma * K * (H2O.^(gamma - 1)); % [C / wt%]
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%

%__________________________________________________________________________
%
%                   dT = depression_dasgupta()
%__________________________________________________________________________
%
% Water follows carbon: CO 2 incites deep silicate melting and dehydration 
% beneath mid-ocean ridges, Geology, 2007
% DOI: 10.1130/G22856A
%
% P in GPa
% CO2 IN MELT in wt %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dTz] = depression_dasgupta(CO2_z)
     dTz = zeros(size(CO2_z)); 
     for iz = 1:numel(CO2_z)
         CO2 = CO2_z(iz);
         if CO2 <= 25
             dT = 27.04 * CO2 +  1490.75 * log((100-1.18 * CO2)/100);
         elseif CO2 > 25 && CO2 < 37;
             dTmax = 27.04 * 25 +  1490.75 * log((100-1.18 * 25)/100);
             dT = dTmax +  (160 - dTmax)/(37 - 25) * (CO2-25);
         elseif CO2 > 37
             dTmax = 27.04 * 25 +  1490.75 * log((100-1.18 * 25)/100);
             dTmax = dTmax +  (160 - dTmax);
             dT = dTmax + 150;
         end
         dTz(iz)=dT;
     end
%      dT = dT * (dT < 300) + 300 * (dT > 300); 
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%