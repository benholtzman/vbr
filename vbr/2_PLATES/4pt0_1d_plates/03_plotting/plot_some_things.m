clear; close all

% choose a thing
  savedir0='../../../../../0_BOXES'; % the closet to store the box in 
  savebase = '2015-12-4-vary_dT'; % boxes end up named ['Box_' savebase]

% choose some other things
%   example plot
%     Target=1550;
%     DepVar='Vbg';
%     DepUnits='cm/yr';
     
    Target = 2; 
    DepVar='Tpot_excess';
    DepUnits='^oC';
    

    t_skip = 10; 
    ylims=[0 120];

%   plot configs
    x_width_LAB=3 ;y_width_LAB=2.5;
    x_width_Example=8 ;y_width_Example=3.5;
    x_width_3=6 ;y_width_3=2;
    plotfigs=[0 0 1];
    savefigs=0*[1 1 1];
  
% load a thing
  loadname=['Box_' savebase '_VBR.mat']; 
  loaddir=[savedir0 '/' savebase '/'];
  addpath ../02_box_functions
  load([loaddir loadname])
  
%   
  nv1 = numel(Box(:,1)); 
  nv2 = numel(Box(1,:)); 
  Val = Box(1,1).Movie.info.settings.(DepVar); 
    titname=[num2str(Val) ' ' DepUnits]; 
  
%% ------------------------------------------------------------------------ 
%% PLOT 1  
%% ------------------------------------------------------------------------
  if plotfigs(1)==1
  figure('color',[1 1 1])
  iv2 = 1; 
  for iv1 = 1:nv1
      tMyr = Box(iv1,iv2).Movie.info.timesteps_myrs;
      nt = numel(tMyr); 
      zLAB = Box(iv1,iv2).Movie.info.zLAB(1:nt)/1e3; 
      zSOL = Box(iv1,iv2).Movie.info.zSOL(1:nt)/1e3; 
      
      if max(zSOL) >= Box(iv1,iv2).Movie.info.Z_km(end-10);
          zSOL = ones(1,nt) * 80;
          clr = [0 0 0];
          lsty='--';
      else
          R = iv1/nv1;
          B = 1-iv1/nv1;
          clr = [R 0 B];
          lsty='-';
      end
      dnm = [num2str(Box(iv1,iv2).info.var1val) Box(iv1,iv2).info.var1units];
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
   var1range=Box(1,1).info.var1range;  
   [minval,iv1]=min(abs(var1range-Target)); 
   tit2name=[titname ', ' num2str(Box(iv1,1).info.var1val) Box(iv1,1).info.var1units];
%  plot the time series
   tMyrs=Box(iv1,1).Movie.info.timesteps_myrs;
   Z = Box(iv1,1).Movie.info.Z_km;
   nt = numel(tMyrs); 
   it_1=2; 
   nts=numel(it_1:t_skip:nt); ip = 1;  
   
   for it = it_1:t_skip:nt   
     F = Box(iv1,1).Movie.Frames(it); 
     clr=[ip/nts 0 ip/nts]; ip = ip + 1; 
     
%    T,Tsol
     subplot(1,3,1)
     if it == it_1; plot(F.Tsol,Z,'--k'); end
     hold all; 
     plot(F.T,Z,'color',clr)
     hold off; 
     xlabel('T [^oC]'); ylabel('z [km]')
     xlim([0 1550]) 
     set(gca,'xtick',[0 500 1000 1500])
%    phi
     subplot(1,3,2)
     if it > it_1; hold all; end
     plot(F.phi,Z,'color',clr)
     if it > it_1; hold off; end
     xlabel('\phi')
     title(tit2name); 
%    visc
     subplot(1,3,3)
     if it > it_1; hold all; end
     semilogx(F.eta,Z,'color',clr)
     if it > it_1; hold off; end
     xlabel('\eta [Pa s]')
     xlim([1e18 1e26])
     set(gca,'xtick',[10.^(18:26)])
     
   end

   for ip=1:3
       subplot(1,3,ip)
       set(gca,'ydir','rev')
       ylim(ylims)
   end

   if savefigs(2) == 1; 
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperPosition', [0 0 x_width_Example y_width_Example]); %
      saveas(gcf,[loaddir 'FIG_Profiles_t.eps'],'epsc')
   end
   end
  
  
 %% ------------------------------------------------------------------------ 
%% PLOT 3
%% ------------------------------------------------------------------------
  if plotfigs(3)==1
  figure('color',[1 1 1])
  iv2 = 1; maxphi=zeros(1,nv1); maxdTdz=zeros(1,nv1); meanphi=zeros(1,nv1); 
  stg=@(f) (f(2:end)+f(1:end-1))/2; 
  for iv1 = 1:nv1
      tMyr = Box(iv1,iv2).Movie.info.timesteps_myrs;
      SOL = Box(iv1,iv2).Movie.info.zSOL/1e3;
      Z = Box(iv1,iv2).Movie.info.Z_km;
      nt = numel(tMyr);       
      
      phimax=0; dTdzmax=0; phimean=0; 
      for it = 1:nt; 
         phi=Box(iv1,iv2).Movie.Frames(it).phi; 
         phimax=max([phimax max(phi)]); 
         
         if max(phi>1e-7)
         phimean(it)=mean(phi(phi>1e-7)); 
         end
         
         T=Box(iv1,iv2).Movie.Frames(it).T; 
         T0=Box(iv1,iv2).Movie.Frames(1).T;
         T=T(Z<(SOL(it))); T0=T0(Z<(SOL(it))); 
         Zz=Z(Z<(SOL(it))); 
         dTdz=(T(2:end) - T(1:end-1))./(Zz(2:end)-Zz(1:end-1));
%          dTdzmax=max([dTdzmax max(dTdz)]);
dTval=5;
         if sum(T>T0+dTval)>0
         dTdzmax(it)=mean(dTdz(stg(T)>(stg(T0)+dTval)));
         end
      end
      dTdzmax=mean(dTdzmax); 
      
      maxphi(iv1)=phimax; 
      meanphi(iv1)=mean(phimean); 
      maxdTdz(iv1)=dTdzmax; 
      
  end
  var1range=Box(1,1).info.var1range; 
  
  subplot(1,2,1)
  plot(var1range,maxphi,'k','marker','.','displayname','max'); 
  hold on
  plot(var1range,meanphi,'r','marker','.','displayname','mean'); 
  hold off
  xlabel([Box(1,1).info.var1name '[' Box(1,1).info.var1units ']'])
  ylabel('\phi')
  legend('location','northwest')
  xlim([var1range(1) var1range(end)])
    subplot(1,2,2)
  plot(var1range,maxdTdz,'k','marker','.'); 
  xlabel([Box(1,1).info.var1name '[' Box(1,1).info.var1units ']'])
  ylabel('mean TBL dT/dz [^oC/km]')
  xlim([var1range(1) var1range(end)])
  if savefigs(3) == 1;
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperPosition', [0 0 x_width_3 y_width_3]); %
      saveas(gcf,[loaddir 'FIG_characterize.eps'],'epsc')
  end
  end
