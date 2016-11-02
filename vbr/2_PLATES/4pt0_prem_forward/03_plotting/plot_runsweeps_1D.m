    
    loadnew='yes';

    
    if strcmp(loadnew,'yes')
        clear;
        loadnew='yes';
    end
    close all

% choose a thing
%  savedir0='../../../../../0_BOXES'; % the closet to store the box in 
  savedir0='/home/chris/Dropbox/Research/0_Boxes'; % the closet to store the box in 
  savebase = '2016-03-24-Sd_dz_25e-2'; % boxes end up named ['Box_' savebase]
    titname='Sd';

% choose some other things
%   example plot
%     Target=1550;
%     DepVar='Vbg';
%     DepUnits='cm/yr';

    Targetv1 = 1; 
    DepVar='Sd';
    DepUnits=' ';
    
    Targett=1; 
    

    t_skip = 2; 
    ylims=[0 150];
    imagesc_clim=[0 35];
    imagesc_yrange=[1400 1550];
%   plot configs
    x_width_LAB=3 ;y_width_LAB=2.5;
    x_width_Example=8.5 ;y_width_Example=3.;
    x_width_3=6 ;y_width_3=2;
    x_width_Cont = 3; y_width_Cont=2; 
    savefigs=[1 1 1];
  
    if strcmp(loadnew,'yes')
        loadname=['Box_' savebase '.mat'];
        loaddir=[savedir0 '/' savebase '/'];
        addpath ../02_box_functions
        load([loaddir loadname])
    end
%   
  nv1 = numel(Box(:,1)); 
  nv2 = numel(Box(1,:)); 
  
  
%% ------------------------------------------------------------------------ 
%% PLOT 1  
%% ------------------------------------------------------------------------
  figure('color',[1 1 1]) 
  for iv1 = 1:nv1
    for iv2 = 1:nv2; 
      tMyr = Box(iv1,iv2).Movie.info.timesteps_myrs;
      nt = numel(tMyr); 
      zLAB = Box(iv1,iv2).Movie.info.zLAB(1:nt)/1e3; 
      zSOL = Box(iv1,iv2).Movie.info.zSOL(1:nt)/1e3; 
      
      [val,id]=min(abs(tMyr-Targett)); 
      id = id * (id > 1) + (id == 1) * 2; 

      Sol_at_Target(iv1,iv2)=zSOL(1)-zSOL(id); 
      
      zkm=Box(iv1,iv2).Movie.info.Z_km;
      dZ=zkm(2)-zkm(1);
      sets=Box(iv1,iv2).Movie.info.settings;
      a=sets.grain0;
      muf=sets.mufo;
      nxi=sets.nxi;
      np=sets.np;
      C=sets.C;
      
      phi=0;phidbl=0;
      for it = 1:numel(tMyr)
          phiz=Box(iv1,iv2).Movie.Frames(it).phi;
          etaz=Box(iv1,iv2).Movie.Frames(it).eta;
          T=Box(iv1,iv2).Movie.Frames(it).T;
          Tsol=Box(iv1,iv2).Movie.Frames(it).Tsol;
          k=phiz.^np.*a*a/C;
          xi = etaz.*(phiz.^nxi);
          dc=sqrt(xi.*k./muf)/1e3;
          
          MOz=max(zkm(T>Tsol));
          [valz,MOid]=min(abs(zkm-MOz));
          
          phi(it)=mean(phiz(T>Tsol));
          
          
          [zval,DBLid]=min(abs(zkm-zSOL(it)));
          dc0=mean(dc(MOid-round(20/dZ):MOid));
         
          
          phidbl(it)=mean(phiz(DBLid:DBLid+round(dc0/dZ)));
          phimax(it)=max(phiz);
      end
      PHI_Mean(iv1,iv2)=mean(phi(phi>1e-5));%(Box(iv1,iv2).Movie.Frames(id).phi);
      PHI_DBL(iv1,iv2)=mean(phidbl);
      PHI_Max(iv1,iv2)=mean(phimax); 
      plotit='yes';
      if max(zSOL) >= Box(iv1,iv2).Movie.info.Z_km(end-10);
          %zSOL = ones(1,nt) * 80;
	  %zSOL=zLAB;
          clr = [0 0 0];
          plotit='no';
          lsty='--';
      else
          R = 1-(iv1-1)/(nv1-1);
          B = (iv1-1)/(nv1-1);
          clr = [R 0 B];
          lsty='-';
      end
      dnm = [num2str(Box(iv1,iv2).info.var1val) Box(iv1,iv2).info.var1units];
      dnmb = [num2str(Box(iv1,iv2).info.var2val) Box(iv1,iv2).info.var2units];
      dnm =[dnm ',' dnmb];
      dnm1 = [dnm ', LAB'];dnm2 = [dnm ', SOL'];
      
%       plot(tMyr,zLAB,'color',clr,'linestyle','--','displayname',dnm1)
%       hold on
      if strcmp(plotit,'yes')
      plot(sqrt(tMyr),zSOL,'color',clr,'linestyle',lsty,'displayname',dnm)
      hold on
      set(gca,'ydir','rev')
      xlabel('sqrt(t) [Myrs^{1/2}]')
      ylabel('z [km]')
      legend('location','southeast')
      xlim([0 sqrt(2)])
      ylim([20 140])
      title(titname)
      end
    end
  end
  if savefigs(1) == 1; 
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperPosition', [0 0 x_width_LAB y_width_LAB]); %
      saveas(gcf,[loaddir 'FIG_zSOL_t.eps'],'epsc')
  end
  
%% ------------------------------------------------------------------------
%% PLOT 2
%% ------------------------------------------------------------------------
   figure('color',[1 1 1])   

%  now find the min  
   var1range=Box(1,1).info.var1range;  
   
   [minval,iv1]=min(abs(var1range-Targetv1));    
   iv2=1;
   tit2name=[titname ', ' num2str(Box(iv1,iv2).info.var1val) Box(iv1,iv2).info.var1units];
%  plot the time series
   tMyrs=Box(iv1,iv2).Movie.info.timesteps_myrs;
   Z = Box(iv1,iv2).Movie.info.Z_km;
   nt = numel(tMyrs); 
   it_1=1; 
   nts=numel(it_1:t_skip:nt); ip = 1;  
   
   for it = it_1:t_skip:nt   
     F = Box(iv1,iv2).Movie.Frames(it); 
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
     xlim([0 0.03])
%    visc
     subplot(1,4,3)
     if it > it_1; hold all; end
     semilogx(F.eta,Z,'color',clr)
     if it > it_1; hold off; end
     xlabel('\eta [Pa s]')
     xlim([1e18 4e26])
     set(gca,'xtick',[10.^(18:26)])
     
     subplot(1,4,4)
     dz =(Z(2)-Z(1))*1e3;
     Kc = (F.Kc(2:end)+F.Kc(1:end-1))/2;
     q=(F.T(2:end)-F.T(1:end-1))/dz .* Kc; 
     zs=(Z(2:end)+Z(1:end-1))/2;
     if it > it_1; hold all; end
     plot(q,zs,'color',clr)
     xlim([0 0.5])
     set(gca,'xtick',[0 0.1 0.2 0.3 0.4 0.5])
    
     if it > it_1; hold off; end
     xlabel('q [W m^2]')     
   end

   for ip=1:4
       subplot(1,4,ip)
       set(gca,'ydir','rev')
       ylim(ylims)
   end
% 
   addpath ../02_functions
   settings=Box(1,1).Movie.info.settings;
   fxc=calc_XtalFactor(Z*1e3,settings);
   subplot(1,4,4)
   hold on
   plot(fxc * 0.4,Z,'linestyle','--','color',[0 0.8 0]); 
   hold off
   
   
 if savefigs(2) == 1; 
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperPosition', [0 0 x_width_Example y_width_Example]); %
      saveas(gcf,[loaddir 'FIG_Profiles_t.eps'],'epsc') 
 end

 figure('color',[1 1 1])
 
 semilogx(var1range(1:nv1),PHI_Mean,'b','marker','.','displayname','in melting column')
 hold on
 semilogx(var1range(1:nv1),PHI_DBL,'k','marker','.','displayname','in DBL')
 semilogx(var1range(1:nv1),PHI_Max,'r','marker','.','displayname','max')
 hold off 
 xlabel([Box(1,1).info.var1name Box(1,1).info.var1units])
 ylabel('mean \phi') 
%  set(gca,'xtick',var1range(1:nv1))
 legend('location','northeast')
%  ylim(imagesc_yrange)
 if savefigs(3) == 1; 
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperPosition', [0 0 x_width_Cont y_width_Cont]); %
      saveas(gcf,[loaddir 'FIG_sweep.eps'],'epsc') 
 end
%   
%   
%  %% ------------------------------------------------------------------------ 
% %% PLOT 3
% %% ------------------------------------------------------------------------
%   figure('color',[1 1 1])
%   iv2 = 1; maxphi=zeros(1,nv1); maxdTdz=zeros(1,nv1); meanphi=zeros(1,nv1); 
%   stg=@(f) (f(2:end)+f(1:end-1))/2; 
%   for iv1 = 1:nv1
%       tMyr = Box(iv1,iv2).Movie.info.timesteps_myrs;
%       SOL = Box(iv1,iv2).Movie.info.zSOL/1e3;
%       Z = Box(iv1,iv2).Movie.info.Z_km;
%       nt = numel(tMyr);       
%       
%       phimax=0; dTdzmax=0; phimean=0; 
%       for it = 1:nt; 
%          phi=Box(iv1,iv2).Movie.Frames(it).phi; 
%          phimax=max([phimax max(phi)]); 
%          
%          if max(phi>1e-7)
%          phimean(it)=mean(phi(phi>1e-7)); 
%          end
%          
%          T=Box(iv1,iv2).Movie.Frames(it).T; 
%          T0=Box(iv1,iv2).Movie.Frames(1).T;
%          T=T(Z<(SOL(it))); T0=T0(Z<(SOL(it))); 
%          Zz=Z(Z<(SOL(it))); 
%          dTdz=(T(2:end) - T(1:end-1))./(Zz(2:end)-Zz(1:end-1));
% %          dTdzmax=max([dTdzmax max(dTdz)]);
% dTval=5;
%          if sum(T>T0+dTval)>0
%          dTdzmax(it)=mean(dTdz(stg(T)>(stg(T0)+dTval)));
%          end
%       end
%       dTdzmax=mean(dTdzmax); 
%       
%       maxphi(iv1)=phimax; 
%       meanphi(iv1)=mean(phimean); 
%       maxdTdz(iv1)=dTdzmax; 
%       
%   end
%   var1range=Box(1,1).info.var1range; 
%   
%   subplot(1,2,1)
%   plot(var1range,maxphi,'k','marker','.','displayname','max'); 
%   hold on
%   plot(var1range,meanphi,'r','marker','.','displayname','mean'); 
%   hold off
%   xlabel([Box(1,1).info.var1name '[' Box(1,1).info.var1units ']'])
%   ylabel('\phi')
%   legend('location','northwest')
%   xlim([var1range(1) var1range(end)])
%     subplot(1,2,2)
%   plot(var1range,maxdTdz,'k','marker','.'); 
%   xlabel([Box(1,1).info.var1name '[' Box(1,1).info.var1units ']'])
%   ylabel('mean TBL dT/dz [^oC/km]')
%   xlim([var1range(1) var1range(end)])
%   if savefigs(3) == 1;
%       set(gcf, 'PaperUnits', 'inches');
%       set(gcf, 'PaperPosition', [0 0 x_width_3 y_width_3]); %
%       saveas(gcf,[loaddir 'FIG_characterize.eps'],'epsc')
%   end
