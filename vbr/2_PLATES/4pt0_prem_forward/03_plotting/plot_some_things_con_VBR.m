clear; close all
% paths and VBR
  addpath ../02_box_functions
  VBR_version = 'VBR_v0p93';
  addpath(['../../../4_VBR/',VBR_version ])
  addpath(['../../../4_VBR/', VBR_version, '/functions'])
  addpath(['../../../4_VBR/', VBR_version, '/params'])

% choose a thing
  savedir0='../../../../../0_BOXES'; % the closet to store the box in 
  savebase = '2015-12-4-vary_dT'; % boxes end up named ['Box_' savebase]

% choose some other things
%   example plot
    Target=1550; % target for plotting the time series
    DepVar='Vbg'; % for title display -- will extract this setting and put it in the title
    DepUnits='cm/yr';
     
%     Target = 2; 
%     DepVar='Tpot_excess';
%     DepUnits='^oC';
    

    t_skip = 10; 
    ylims=[0 120];

%   plot configs
    x_width_LAB=3 ;y_width_LAB=2.5;
    x_width_Example=8 ;y_width_Example=3.5;
    x_width_3a=5 ;y_width_3a=2;
    x_width_3b=8 ;y_width_3b=2;
    plotfigs=[1 1 0];
    savefigs=[0 1 0];
  
% load a thing
  loadname=['Box_' savebase '_VBR.mat']; 
  loaddir=[savedir0 '/' savebase '/'];
  addpath ../02_box_functions
  load([loaddir loadname])
  Box1 = Box;
  
% load another thing
  loadname=['Box_' savebase '_VBR_nomelt.mat']; 
  loaddir=[savedir0 '/' savebase '/'];
  addpath ../02_box_functions
  load([loaddir loadname])
  Box2 = Box;
%   
  nv1 = numel(Box1(:,1)); 
  nv2 = numel(Box1(1,:)); 
  Val = Box1(1,1).Movie.info.settings.(DepVar); 
    titname=[num2str(Val) ' ' DepUnits]; 
  
%% ------------------------------------------------------------------------ 
%% PLOT 1  
%% ------------------------------------------------------------------------
  if plotfigs(1)==1
  figure('color',[1 1 1])
  iv2 = 1; 
  for iv1 = 1:nv1
      tMyr = Box1(iv1,iv2).Movie.info.timesteps_myrs;
      nt = numel(tMyr); 
      zLAB = Box1(iv1,iv2).Movie.info.zLAB(1:nt)/1e3; 
      zSOL = Box1(iv1,iv2).Movie.info.zSOL(1:nt)/1e3; 
      
      if max(zSOL) >= Box1(iv1,iv2).Movie.info.Z_km(end-10);
          zSOL = ones(1,nt) * 80;
          clr = [0 0 0];
          lsty='--';
      else
          R = iv1/nv1;
          B = 1-iv1/nv1;
          clr = [R 0 B];
          lsty='-';
      end
      dnm = [num2str(Box1(iv1,iv2).info.var1val) Box1(iv1,iv2).info.var1units];
      dnm1 = [dnm ', LAB'];dnm2 = [dnm ', SOL'];
      
%       plot(tMyr,zLAB,'color',clr,'linestyle','--','displayname',dnm1)
%       hold on
      plot((tMyr),zSOL,'color',clr,'linestyle',lsty,'displayname',dnm)
      hold on
      set(gca,'ydir','rev')
      xlabel('t [Myrs]')
      ylabel('z [km]')
      legend('location','northwest')
      xlim([0 10])
      ylim([30 90])
      title(titname)
  end
  if savefigs(1) == 1; 
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperPosition', [0 0 x_width_LAB y_width_LAB]); %
      saveas(gcf,[loaddir 'FIG_zSOL_t.eps'],'epsc')
  end
  end
  
%% ------------------------------------------------------------------------
%% PLOT 2
%% ------------------------------------------------------------------------
  if plotfigs(2)==1
   figure('color',[1 1 1])   

%  now find the min  
   var1range=Box1(1,1).info.var1range;  
   [minval,iv1]=min(abs(var1range-Target)); 
   tit2name=[titname ', ' num2str(Box1(iv1,1).info.var1val) Box1(iv1,1).info.var1units];
%  plot the time series
   tMyrs=Box1(iv1,1).Movie.info.timesteps_myrs;
   Z = Box1(iv1,1).Movie.info.Z_km;
   nt = numel(tMyrs); 
   it_1=2; 
   nts=numel(it_1:t_skip:nt); ip = 1;  
   
   for it = it_1:t_skip:nt   
     F = Box1(iv1,1).Movie.Frames(it); 
     clr=[ip/nts 0 ip/nts]; ip = ip + 1;          
     
%    T,Tsol
     subplot(1,4,1)
     if it == it_1; plot(F.Tsol,Z,'--k'); end
     hold all; 
     plot(F.T,Z,'color',clr)
     hold off; 
     xlabel('T [^oC]'); ylabel('z [km]')
     xlim([0 1550]) 
     set(gca,'xtick',[0 500 1000 1500])
%    phi
     subplot(1,4,2)
     if it > it_1; hold all; end
     plot(F.phi,Z,'color',clr)
     if it > it_1; hold off; end
     xlabel('\phi')
     title(tit2name); 
%    visc
     subplot(1,4,3)
     if it > it_1; hold all; end
     semilogx(F.eta,Z,'color',clr)
     if it > it_1; hold off; end
     xlabel('\eta [Pa s]')
     xlim([1e18 1e26])
     set(gca,'xtick',[10.^(18:26)])
     
%    velocity 
%    first get depth varying Vs, need to average over the frequency range 
     Vs_full = Box1(iv1,iv2).Movie.Frames(it).VBR.AndradePsP.Va/1e3;        
     Vs_z = zeros(size(F.phi));
     for iz = 1:numel(F.phi)
         Vs_z(iz) = mean(Vs_full(:,iz));
     end
     
     subplot(1,4,4)
     if it > it_1; hold all; end
     plot(Vs_z,Z,'color',clr)
     if it > it_1; hold off; end
     xlabel('V_s [km/s]')     
     
   end

   for ip=1:4
       subplot(1,4,ip)
       set(gca,'ydir','rev')
       ylim(ylims)
   end

   if savefigs(2) == 1; 
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperPosition', [0 0 x_width_Example y_width_Example]); %
      saveas(gcf,[loaddir 'FIG_Profiles_t_' num2str(Target) '.eps'],'epsc')
   end
   end
  
%% ------------------------------------------------------------------------ 
%% PLOT 3
%% ------------------------------------------------------------------------
  if plotfigs(3)==1
  f1=figure('color',[1 1 1]);
  f2=figure('color',[1 1 1]);
  iv2 = 1; 
  stg=@(f) (f(2:end)+f(1:end-1))/2;

  dZcal = [2 3.5 5 10 15 30]; % average over this lengthscale about the boundary
  ndZcal=numel(dZcal);
% average values for each run
  maxphi=zeros(1,nv1); 
  meanphi=zeros(ndZcal,nv1); 
  dV_phi=zeros(ndZcal,nv1);
  dVdz_Tphi=zeros(ndZcal,nv1);
  meandTdz=zeros(ndZcal,nv1); 
  dVdz_T=zeros(ndZcal,nv1);

% loop over each run
  for iv1 = 1:nv1
      tMyr = Box1(iv1,iv2).Movie.info.timesteps_myrs;
      SOL = Box1(iv1,iv2).Movie.info.zSOL/1e3;
      Z = Box1(iv1,iv2).Movie.info.Z_km;
      nt = numel(tMyr);       

      
%     initialize average values for time steps
      dV_ave = zeros(ndZcal,nt);       
      dTdzmax = zeros(1,nt);
      dTdzmean = zeros(ndZcal,nt);
      dVdzmean = zeros(ndZcal,nt);
      dVdz1mean = dVdzmean;
      phimax = zeros(1,nt); 
      phimean = zeros(ndZcal,nt);

%     get average values for each time step
      for it = 1:nt; 

%        extract melt fraction      
         phi=Box1(iv1,iv2).Movie.Frames(it).phi; 
         ZLAB=Box1(iv1,iv2).Movie.info.zLAB(it)/1e3;
         ZSOL=Box1(iv1,iv2).Movie.info.zSOL(it)/1e3;
         
%        extract the temperature in the thermal boundary layer
         T=Box1(iv1,iv2).Movie.Frames(it).T; 
         T0=Box1(iv1,iv2).Movie.Frames(1).T;
         Tsol=Box1(iv1,iv2).Movie.Frames(1).Tsol; 
         nsol = sum(T0>Tsol); 

%        get VBR Vs with and without melt
         Vs_full = Box1(iv1,iv2).Movie.Frames(it).VBR.AndradePsP.Va/1e3;
         Vs_NoMelt_full = Box2(iv1,iv2).Movie.Frames(it).VBR.AndradePsP.Va/1e3;                          

%        construct depth varying Vs with/without melt, need to average over
%        the frequency range calculated for each 
         Vs_Melt = zeros(size(phi));
         Vs_NoMelt = zeros(size(phi));
         for iz = 1:numel(phi)
           Vs_Melt(iz) = mean(Vs_full(:,iz));
           Vs_NoMelt(iz) = mean(Vs_NoMelt_full(:,iz));
         end

%        calculate average drop due to melt (should be identical
%        where there isn't any melt)
         dV = (Vs_NoMelt - Vs_Melt)./Vs_NoMelt; 
         
         
%        max/mean phi where partially molten	 
         if nsol > 0
             phimax(it)=max(phi);
             for izcal = 1:ndZcal         
                 phitoave = phi.*(Z - ZSOL <= 2*dZcal(izcal)).*(T>=Tsol); 
                 phimean(izcal,it)=mean(phitoave(phitoave~=0));
                 dVtoave = dV.*(Z - ZSOL <= 2*dZcal(izcal)).*(T>=Tsol); 
                 dV_ave(izcal,it) = mean(dVtoave(dVtoave~=0)); 
             end
         end            
                          
%        calculate T-related gradients 
         dTdz0=(T(2:end) - T(1:end-1))./(Z(2:end)-Z(1:end-1));
         dVdz0=(Vs_NoMelt(2:end)-Vs_NoMelt(1:end-1))./(Z(2:end)-Z(1:end-1));
         dVdz10=(Vs_Melt(2:end)-Vs_Melt(1:end-1))./(Z(2:end)-Z(1:end-1));
         Zs = stg(Z); 
         
%        find some way to average... 
         for izcal = 1:ndZcal
             if nsol > 0
                 dTdz = dTdz0 .* (abs(Zs - ZSOL) <= dZcal(izcal));
                 dVdz = dVdz0 .* (abs(Zs - ZSOL) <= dZcal(izcal));
                 dVdz1 = dVdz10 .* (abs(Zs - ZSOL) <= dZcal(izcal));
             else
                 dTdz = dTdz0 .* (abs(Zs - ZLAB) <= dZcal(izcal));
                 dVdz = dVdz0 .* (abs(Zs - ZLAB) <= dZcal(izcal));
                 dVdz1 = dVdz10 .* (abs(Zs - ZLAB) <= dZcal(izcal));
             end             
             dVdzmean(izcal,it)=mean(dVdz(dVdz~=0));             
             dTdzmean(izcal,it)=mean(dTdz(dTdz~=0));             
             dVdz1mean(izcal,it)=mean(dVdz1(dVdz1~=0));   
         end
         
         
         
         
      end
      
%     average over the time steps       
      maxphi(iv1)=max(phimax); 
             
      for izcal = 1:ndZcal
          meandTdz(izcal,iv1)=mean(dTdzmean(izcal,:));
          dVdz_T(izcal,iv1)=mean(dVdzmean(izcal,:));
          dVdz_Tphi(izcal,iv1)=mean(dVdz1mean(izcal,:));
          meanphi(izcal,iv1)=mean(phimean(izcal,:)); 
          dV_phi(izcal,iv1)=mean(dV_ave(izcal,:));
      end
  end
  var1range=Box1(1,1).info.var1range;
  
  figure(f1)
  subplot(1,2,1)
  plot(var1range,maxphi,'k','marker','.','displayname','max');
  
  for izcal = 1:ndZcal
  clr = [(izcal-1)/ndZcal 0 (ndZcal-izcal)/ndZcal];
  dnm= num2str(dZcal(izcal)*2);    
  
  figure(f1)
  subplot(1,2,1)  
    hold on
    plot(var1range,meanphi(izcal,:),'color',clr,'marker','.','displayname',dnm);
    hold off
    xlabel([Box1(1,1).info.var1name '[' Box1(1,1).info.var1units ']'])
    ylabel('\phi')
    legend('location','northwest')
    xlim([var1range(1) var1range(end)])
  
    subplot(1,2,2)  
    if izcal> 1; hold on; end
    plot(var1range,meandTdz(izcal,:),'color',clr,'marker','.');
    xlabel([Box1(1,1).info.var1name '[' Box1(1,1).info.var1units ']'])
    ylabel('mean TBL dT/dz [^oC/km]')
    xlim([var1range(1) var1range(end)])
    ylim([min(min(meandTdz)) max(max(meandTdz))])
    
  figure(f2)
  subplot(1,3,1) 
    if izcal> 1; hold on; end
    plot(var1range,dVdz_Tphi(izcal,:),'color',clr,'marker','.','displayname',dnm);    
    ylabel('mean dV/dz, with melt [(km/s)/km]')
    xlabel([Box1(1,1).info.var1name '[' Box1(1,1).info.var1units ']'])
    xlim([var1range(1) var1range(end)])
    ylim([-.08 0.01])
         
  subplot(1,3,2)
    if izcal> 1; hold on; end
    plot(var1range,dVdz_T(izcal,:),'color',clr,'marker','.','displayname',dnm);
    xlabel([Box1(1,1).info.var1name '[' Box1(1,1).info.var1units ']'])
    ylabel('mean dV/dz, no melt [(km/s)/km]')
    xlim([var1range(1) var1range(end)])
    ylim([-.08 0.01])
    legend('location','southwest')
    
  subplot(1,3,3) 
    if izcal> 1; hold on; end
    plot(var1range,dVdz_Tphi(izcal,:)./dVdz_T(izcal,:),'color',clr,'marker','.','displayname',dnm);    
    ylabel('with melt / without melt')
    xlabel([Box1(1,1).info.var1name '[' Box1(1,1).info.var1units ']'])
    xlim([var1range(1) var1range(end)])     
  end
  
  
  if savefigs(3) == 1;
      set(f1, 'PaperUnits', 'inches');
      set(f1, 'PaperPosition', [0 0 x_width_3a y_width_3a]); %
      saveas(f1,[loaddir 'FIG_characterize_a.eps'],'epsc')
      
      set(f2, 'PaperUnits', 'inches');
      set(f2, 'PaperPosition', [0 0 x_width_3b y_width_3b]); %
      saveas(f2,[loaddir 'FIG_characterize_b.eps'],'epsc')
  end
  end
