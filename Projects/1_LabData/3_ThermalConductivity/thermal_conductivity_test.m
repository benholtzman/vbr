function thermal_conductivity_test  
  close all
  addpath ../../2_PLATES/4pt0_prem_forward/01_functions/
  
% load grabit points 
  load Xu_fig9

  T = linspace(200,1400,100);  
  
% minimize for Kc_o
  [Kc_o,KcMisfit,KcTest] = Minimize_Kco;

% calculate conductivity   
  [Rho,Cp,Kc_10,P] = MaterialProperties(3300,Kc_o,1000,T,1e5,10*1e9,0.5,'Prescribed_P');
  [Rho,Cp,Kc_7,P] = MaterialProperties(3300,Kc_o,1000,T,1e5,7*1e9,0.5,'Prescribed_P');
  [Rho,Cp,Kc_4,P] = MaterialProperties(3300,Kc_o,1000,T,1e5,4*1e9,0.5,'Prescribed_P');

figure('color',[1 1 1])
subplot(1,2,1)
plot(k_10GPa_T,k_10GPa,'.r','displayname','10 GPa','markersize',16)
hold on 
plot(k_7GPa_T,k_7GPa,'.b','displayname','7 GPa','markersize',16)
plot(k_4GPa_T,k_4GPa,'.k','displayname','4 GPa','markersize',16)
plot(T,Kc_10,'r','displayname','calculated')
plot(T,Kc_7,'b','displayname','calculated')
plot(T,Kc_4,'k','displayname','calculated')
xlim([200 1400])

xlabel('T [^oC]')
ylabel('k [W m^{-1} K^{-1}]')
legend('location','northeast')

subplot(1,2,2)
semilogy(KcTest,KcMisfit,'k')
hold on 
semilogy(Kc_o,KcMisfit(KcTest==Kc_o),'.k','markersize',16)
xlabel('test k_o [W m^{-1} K^{-1}')
ylabel('residual')
title(['Best fit Kc_o: ' num2str(Kc_o) ' W m^{-1} K^{-1}'])
end

function [Kc_o,KcMisfit,KcTest] = Minimize_Kco
%  simple minimzation
   KcTest=linspace(2,6,1000); 
   KcMisfit = zeros(size(KcTest));
   
   load Xu_fig9

   for iKco = 1:numel(KcTest)
       Kc_o=KcTest(iKco);
       [Rho,Cp,Kc_10,P] = MaterialProperties(3300,Kc_o,1000,k_10GPa_T,1e5,10*1e9,0.5,'Prescribed_P');
       [Rho,Cp,Kc_7,P] = MaterialProperties(3300,Kc_o,1000,k_7GPa_T,1e5,7*1e9,0.5,'Prescribed_P');
       [Rho,Cp,Kc_4,P] = MaterialProperties(3300,Kc_o,1000,k_4GPa_T,1e5,4*1e9,0.5,'Prescribed_P');
       
       err_10=sum(abs(Kc_10 - k_10GPa)./k_10GPa);
       err_7=sum(abs(Kc_7 - k_7GPa)./k_7GPa);
       err_4=sum(abs(Kc_4 - k_4GPa)./k_4GPa);
       
       KcMisfit(iKco)=err_10+err_7+err_4; 
   end

   [minval,id]=min(KcMisfit);    
   Kc_o=KcTest(id); 
end