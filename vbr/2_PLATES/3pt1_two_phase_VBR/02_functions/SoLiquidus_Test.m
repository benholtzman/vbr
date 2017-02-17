function SoLiquidus_Test
z = linspace(0,200,100)*1e3; 
P = 3300*9.8*z;



T = 1325  + 0.4 *1e-3 * z; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot type 1: constant F, plot Tsol vs P (or z) for varying concentrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = 0.0001; 
% no carbon
  Cs_H2O = [0 100 200 300 400 500] * 1e-4; % weight percent
  Cs_CO2 = 0 * 1e-4; % weight percent   
  Tsol = zeros(numel(P),numel(Cs_H2O)); 
  for iH2O = 1:numel(Cs_H2O) 
    Solidus = get_Solidus(Cs_H2O(iH2O),Cs_CO2,P,F);
    Tsol(:,iH2O) = Solidus.Tsol; 
  end
  figure
  subplot(1,3,1)
  plot(Tsol,z/1e3); set(gca,'ydir','rev'); 
  legendCell=strcat(strtrim(cellstr(num2str(Cs_H2O(:)*1e4))),' ppm H_2O');
  legend(legendCell)
  hold on
  plot(T,z/1e3,'--k')
  
% no water
  Cs_H2O = 0* 1e-4; % weight percent
  Cs_CO2 = [0 25 50 75 100 150] * 1e-4; % weight percent   
  Tsol = zeros(numel(P),numel(Cs_CO2)); 
  for iCO2 = 1:numel(Cs_CO2) 
    Solidus = get_Solidus(Cs_H2O,Cs_CO2(iCO2),P,F);
    Tsol(:,iCO2) = Solidus.Tsol; 
  end
  subplot(1,3,2)
  plot(Tsol,z/1e3); set(gca,'ydir','rev'); 
  legendCell=strcat(strtrim(cellstr(num2str(Cs_CO2(:)*1e4))),' ppm CO_2');
  legend(legendCell) 
  hold on
  plot(T,z/1e3,'--k')
  
% both
  Cs_H2O = [0 100 0 100]* 1e-4; % weight percent
  Cs_CO2 = [0 0 100 100] * 1e-4; % weight percent   
  Tsol = zeros(numel(P),numel(Cs_CO2)); 
  for iCO2 = 1:numel(Cs_CO2) 
    Solidus = get_Solidus(Cs_H2O(iCO2),Cs_CO2(iCO2),P,F);
    Tsol(:,iCO2) = Solidus.Tsol; 
  end
    
  subplot(1,3,3)
  plot(Tsol,z/1e3); set(gca,'ydir','rev'); 
  legendCell=strcat(strtrim(cellstr(num2str(Cs_CO2(:)*1e4))),' ppm CO_2,',...
              strtrim(cellstr(num2str(Cs_H2O(:)*1e4))),' ppm H_2O');
  legend(legendCell) 
  hold on
  plot(T,z/1e3,'--k')
end


function Solidus = get_Solidus(Cs_H2O,Cs_CO2,P,F)


  kd_H2O = 1e-2; 
  kd_CO2 = 1e-4; 
  

% F = linspace(0.00001,1,1000); 
Cf_H2O = Cs_H2O / (kd_H2O + F * (1-kd_H2O)); 
Cf_CO2 = Cs_CO2 ./ (kd_CO2+ F * (1-kd_CO2)); 



%F + D * (1-F)

[Solidus] = SoLiquidus(P,Cf_H2O,Cf_CO2,'hirschmann');
end