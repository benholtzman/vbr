function thermal_conductivity_test
  close all

  % put VBR in the path
  path_to_top_level_vbr='../../../';
  addpath(path_to_top_level_vbr)
  vbr_init

  % data if it works
  data=tryLoad();

  % calculate a basic fit 
  T = linspace(200,1400,100);
  if isfield(data,'k_4GPa')
    % minimize for Kc_o
    [Kc_o,KcMisfit,KcTest] = Minimize_Kco(data);
  else
    Kc_0=4.17;
  end

  % calculate conductivity
  [Rho,Cp,Kc_10,P] = MaterialProperties(3300,Kc_o,1000,T,1e5,10*1e9,0.5,'Prescribed_P');
  [Rho,Cp,Kc_7,P] = MaterialProperties(3300,Kc_o,1000,T,1e5,7*1e9,0.5,'Prescribed_P');
  [Rho,Cp,Kc_4,P] = MaterialProperties(3300,Kc_o,1000,T,1e5,4*1e9,0.5,'Prescribed_P');

  figure('color',[1 1 1])
  if isfield(data,'k_4GPa')
    subplot(1,2,1)
  end
  plot(T,Kc_10,'r','displayname','calculated')
  hold on
  plot(T,Kc_7,'b','displayname','calculated')
  plot(T,Kc_4,'k','displayname','calculated')
  if isfield(data,'k_4GPa')
    plot(data.k_10GPa_T,data.k_10GPa,'.r','displayname','10 GPa','markersize',16)
    plot(data.k_7GPa_T,data.k_7GPa,'.b','displayname','7 GPa','markersize',16)
    plot(data.k_4GPa_T,data.k_4GPa,'.k','displayname','4 GPa','markersize',16)
  end
  xlim([200 1400])
  xlabel('T [^oC]')
  ylabel('k [W m^{-1} K^{-1}]')
  legend('location','northeast')

  if isfield(data,'k_4GPa')
    subplot(1,2,2)
    semilogy(KcTest,KcMisfit,'k')
    hold on
    semilogy(Kc_o,KcMisfit(KcTest==Kc_o),'.k','markersize',16)
    xlabel('test k_o [W m^{-1} K^{-1}')
    ylabel('residual')
    title(['Best fit Kc_o: ' num2str(Kc_o) ' W m^{-1} K^{-1}'])
  end
end

function [Kc_o,KcMisfit,KcTest] = Minimize_Kco(data)
  %  simple minimzation
  KcTest=linspace(2,6,1000);
  KcMisfit = zeros(size(KcTest));

  for iKco = 1:numel(KcTest)
    Kc_o=KcTest(iKco);
    [Rho,Cp,Kc_10,P] = MaterialProperties(3300,Kc_o,1000,data.k_10GPa_T,1e5,10*1e9,0.5,'Prescribed_P');
    [Rho,Cp,Kc_7,P] = MaterialProperties(3300,Kc_o,1000,data.k_7GPa_T,1e5,7*1e9,0.5,'Prescribed_P');
    [Rho,Cp,Kc_4,P] = MaterialProperties(3300,Kc_o,1000,data.k_4GPa_T,1e5,4*1e9,0.5,'Prescribed_P');

    err_10=sum(abs(Kc_10 - data.k_10GPa)./data.k_10GPa);
    err_7=sum(abs(Kc_7 - data.k_7GPa)./data.k_7GPa);
    err_4=sum(abs(Kc_4 - data.k_4GPa)./data.k_4GPa);

    KcMisfit(iKco)=err_10+err_7+err_4;
  end

  [minval,id]=min(KcMisfit);
  Kc_o=KcTest(id);
end

function data=tryLoad()
  dataDir='../../../../vbrWork/expt_data/thermalconductivity/';
  data=struct();
  if exist([dataDir,'Xu_fig9.mat'],'file')
    load([dataDir,'Xu_fig9.mat'])
    data.k_10GPa=k_10GPa;
    data.k_10GPa_T=k_10GPa_T;
    data.k_4GPa=k_4GPa;
    data.k_4GPa_T=k_4GPa_T;
    data.k_7GPa=k_7GPa;
    data.k_7GPa_T=k_7GPa_T;
  end
end
