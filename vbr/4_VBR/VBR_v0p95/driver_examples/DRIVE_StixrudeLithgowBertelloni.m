   close all; clear
  
%% add required paths  
   addpath ../functions
   addpath ../params
   addpath ../
    
%% construct state variable fields          
   z = linspace(0,450,500)*1e3; z= z';     
   to = 110 * (1e6 * 365 * 24 * 3600); % Myrs * ( Myrs_to_seconds)   
   VBR.in.SV.T_K = 1350*erf(z./(2*sqrt(1e-6.*to))) + 273 + z/1e3 * 0.4;         
   VBR.in.SV.sig_MPa = 0.1*ones(size(z));
   VBR.in.SV.dg_um = 0.01 * 1e6 * ones(size(z)); 
   VBR.in.SV.phi = 0.000001 * (z >= 100000) .* (z<=300000);
   VBR.in.SV.Ch2o = zeros(size(z)); % in PPM!      
   VBR.in.SV.chi = 0.5*(1+ erf((z/1e3-6)/5));
   VBR.in.SV.f = [1 35]/1e3; % 35 mHz
                  
   
%% write method list (these are the things to calculate)
   VBR.in.elastic.methods_list={'anharmonic';'poro_Takei';'SLB2005'};    
   VBR.in.anelastic.methods_list={'eBurgers';'AndradePsP'};        
   
%% load anharmonic parameters, adjust Gu_0_ol   
   VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity 
   VBR.in.elastic.anharmonic.Gu_0_ol = 75.5; % olivine reference shear modulus [GPa]
   
   
%% calculate density and pressure
   addpath ../../../2_PLATES/4pt0_prem_forward/01_functions/
   Rho_o = (3300-2800) * VBR.in.SV.chi +2800;   
   [Rho,Cp,Kc,P] = MaterialProperties(Rho_o,ones(size(z)),ones(size(z)),VBR.in.SV.T_K,z,0,0.4/1e3,'PT_dep');
   VBR.in.SV.rho=Rho; 
   VBR.in.SV.P_GPa=P/1e9; 

%% call VBR    
   [VBR] = VBR_spine(VBR);  
     
%% load reference models, remove water layers
   VeloModDir='../../../6_FitVobs/velocity_models/';
   PA5=load([VeloModDir 'PA5/PA5.mat']);
   PAC=load([VeloModDir 'PAC/PAC.mat']);
      
   PA5.Vs_ave = (PA5.Vsh_kms + PA5.Vsv_kms)/2;  % average shear wave
   PA5.z_km=PA5.z_km(PA5.Vs_ave>0); % remove water layer
   PA5.z_km = PA5.z_km - PA5.z_km(1); 
   PA5.Vs_ave=PA5.Vs_ave(PA5.Vs_ave>0); % remove water layer
   
   PAC.z_km=PAC.z_km(PAC.Vs_kms>0); % remove water layer
   PAC.z_km = PAC.z_km - PAC.z_km(1); 
   PAC.Vs_kms=PAC.Vs_kms(PAC.Vs_kms>0); % remove water layer
      
%% plot
   close all; figure('color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
   subplot(1,3,1)
   plot(PAC.Vs_kms,PAC.z_km,'k','displayname','PAC')
   hold on
   plot(PA5.Vs_ave,PA5.z_km,'--k','displayname','PA5')   
   plot(VBR.out.elastic.anharmonic.Vsu/1e3,z/1e3,'m','displayname','Anharmonic')  
   
   for ifreq = 1:numel(VBR.in.SV.f)
       scl=1-(ifreq-1)/(numel(VBR.in.SV.f)-1); 
       clr = [0 scl (1-scl)]; 
       dn = ['eBurg, f=' num2str(VBR.in.SV.f(ifreq)*1000) ' mHz'];
       Vs=squeeze(VBR.out.anelastic.eBurgers.V(:,ifreq)/1e3);
       plot(Vs,z/1e3,'color',clr,'displayname',dn)
   end
   title('VBR calculator anharmonic and extended burgers')
   
   subplot(1,3,2)
   plot(PAC.Vs_kms,PAC.z_km,'k','displayname','PAC')
   hold on
   plot(PA5.Vs_ave,PA5.z_km,'--k','displayname','PA5')
   plot(VBR.out.elastic.SLB2005.Vs,z/1e3,'r','displayname','SLB2005 - LVZ')   
   title('Stixrude and Lithgow-Bertelloni 2005 LVZ parametrization')
   
   subplot(1,3,3)
   plot(PAC.Vs_kms,PAC.z_km,'k','displayname','PAC')
   hold on
   plot(PA5.Vs_ave,PA5.z_km,'--k','displayname','PA5')      
   for ifreq = 1:numel(VBR.in.SV.f)
       scl=1-(ifreq-1)/(numel(VBR.in.SV.f)-1); 
       clr = [0 scl (1-scl)]; 
       dn = ['eBurg, f=' num2str(VBR.in.SV.f(ifreq)*1000) ' mHz'];
       Vs=squeeze(VBR.out.anelastic.eBurgers.V(:,ifreq)/1e3);
       plot(Vs,z/1e3,'color',clr,'displayname',dn)
   end
   title('comparison')
   
   
   plot(VBR.out.elastic.SLB2005.Vs,z/1e3,'r','displayname','SLB2005 - LVZ')   
   

   
   for ip=1:3
       subplot(1,3,ip)
       legend('location','southwest')
       set(gca,'ydir','rev')
       xlim([4.2 4.8])
       ylim([0 max(z/1e3)])
       xlabel('V_s [km/s]')
   end
   ylabel('depth [km]')
