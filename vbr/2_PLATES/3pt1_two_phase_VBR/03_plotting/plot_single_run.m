%  load (or reload) a box?     
   loadnew='yes'; % 'yes' will load     
   if strcmp(loadnew,'yes');clear;loadnew='yes'; end
   close all

% choose a thing
  cwd=pwd; cd ~; hmdir=pwd; cd(cwd);
  loaddir0=[hmdir '/Desktop']; % the closet to store the box in 
  loadbase = '2016-02-16-TEST'; % boxes end up named ['Box_' savebase]
  addpath ./postproc/  
  
% choose another thing
  savedir='/home/chris/Desktop/'; % where to save the figures
   
% (plot 1) zLAB, zSOL vs t
   Fig(1).box_indeces='both'; % 'both', 'var1', 'var2'
   Fig(1).box_indeces_v1_skip=1; % increment for box index
   Fig(1).box_indeces_v2_skip=1; % increment for box index
   Fig(1).paper_x_width = 3; % figure width for save (inches)
   Fig(1).paper_y_width = 2.5; % figure height for save (inches) 
   Fig(1).savefig='yes'; % 'yes' or 'no'
   Fig(1).savedir=savedir; % filename (saved in box location)
   Fig(1).savename='FIG_zSOL_t.eps'; % filename (saved in box location)   
   
% (plot 2) time evolution of single run (with heat flow)
   Targetv1 = 80; % variable 1 value
   Targetv2 = 1525; % variable 2 value (will find the closest box index)
   t_skip = 1; % time step skip for time evolution plot
   ylims=[0 120];
   x_width_Example=8.5 ;y_width_Example=3.;   
   Fig(2).paper_x_width = 8.5; % figure width for save (inches)
   Fig(2).paper_y_width = 3.0; % figure height for save (inches)
   Fig(2).savefig='yes'; % 'yes' or 'no'
   Fig(2).savedir=savedir; % filename (saved in box location)
   Fig(2).savename='FIG_Profiles_t.eps'; % filename (saved in box location)

% (plot 5) time evolution of single run with upwelling velo instead of heat flow 
   Fig(3).paper_x_width = 8.5; % figure width for save (inches)
   Fig(3).paper_y_width = 3.0; % figure height for save (inches)   
   Fig(3).savefig='yes'; % 'yes' or 'no'
   Fig(3).savedir='/home/chris/Desktop/'; % filename (saved in box location)
   Fig(3).savename='FIG_Profiles_t_Vbg.eps'; % filename (saved in box location)   
   
   
%% ------------------------------------------------------------------------  
%% Load the box
%% ------------------------------------------------------------------------  
    if strcmp(loadnew,'yes')
        loadname=['Box_' loadbase '.mat'];
        loaddir=[loaddir0 '/' loadbase '/'];
        addpath ../02_box_functions
        load([loaddir loadname])
    end
    
    Sets = Box(1,1).run_info.settings; 
    titlebase = 'Figure';
           
%% ------------------------------------------------------------------------
%% PLOT 1 - zLAB vs t 
%% ------------------------------------------------------------------------   
   nv1 = numel(Box(:,1)); 
   nv2 = numel(Box(1,:)); 
   if strcmp(Fig(1).box_indeces,'both')
       var1range=1:Fig(1).box_indeces_v1_skip:nv1;
       var2range=1:Fig(1).box_indeces_v2_skip:nv2;
   elseif strcmp(Fig(1).box_indeces,'var1')
       var1range=1:Fig(1).box_indeces_v1_skip:nv1;
       var2range=Fig(1).box_indeces_v2_skip;
   elseif strcmp(Fig(1).box_indeces,'var2')
       var1range=Fig(1).box_indeces_v1_skip;
       var2range=1:Fig(1).box_indeces_v2_skip:nv2;
   end
   
   Fig(1).fig=figure('color',[1 1 1]); iplt = 1; 
   for iv1 = var1range
       for iv2 = var2range
           tVars = get_q_phi_dVdz_vs_t(Box,iv1,iv2,'calc_melt');

           dnm = num2str([Box(iv1,iv2).info.var1val Box(iv1,iv2).info.var2val],2);
                         
           if iplt > 1; hold all; end
           plot(tVars.tMyr,tVars.zSOL,'displayname',dnm)
           if iplt > 1; hold off; end
           iplt = iplt+1; 
       end      
   end
   
   legend('location','southeast')
   set(gca,'ydir','rev')
   xlabel('t [Myrs]')
   ylabel('z_{SOL} [km]')    
  
%% ------------------------------------------------------------------------
%% PLOT 2 - time evolution
%% ------------------------------------------------------------------------
   Fig(2).fig=figure('color',[1 1 1]);

%  now find the min  
   var1range=Box(1,1).info.var1range;  
   var2range=Box(1,1).info.var2range;  
   [minval,iv1]=min(abs(var1range-Targetv1)); 
   [minval,iv2]=min(abs(var2range-Targetv2));
   if isempty(iv2); iv2=1; end
   tit2name=[titlebase ', ' num2str(Box(iv1,iv2).info.var1val) ' ' Box(iv1,iv2).info.var1units ...
             ',' num2str(Box(iv1,iv2).info.var2val) ' ' Box(iv1,iv2).info.var2units];
%  plot the time series
   tMyrs=Box(iv1,iv2).run_info.timesteps_myrs;
   Z = Box(iv1,iv2).run_info.Z_km;
   nt = numel(tMyrs); 
   it_1=1; 
   nts=numel(it_1:t_skip:nt); ip = 1;  
   
   for it = it_1:t_skip:nt   
     F = Box(iv1,iv2).Frames(it); 
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
   settings=Box(1,1).run_info.settings;
   fxc=calc_XtalFactor(Z*1e3,settings);
   subplot(1,4,4)
   hold on
   plot(fxc * 0.4,Z,'linestyle','--','color',[0 0.8 0]); 
   hold off      
 
   
%% ------------------------------------------------------------------------
%% PLOT 3 - time evolution, slightly different vars plotted
%% ------------------------------------------------------------------------
   Fig(3).fig=figure('color',[1 1 1]);

%  now find the min  
   var1range=Box(1,1).info.var1range;  
   var2range=Box(1,1).info.var2range;  
   [minval,iv1]=min(abs(var1range-Targetv1)); 
   [minval,iv2]=min(abs(var2range-Targetv2));
   if isempty(iv2); iv2=1; end
   tit2name=[titlebase ', ' num2str(Box(iv1,iv2).info.var1val) ' ' Box(iv1,iv2).info.var1units ...
             ',' num2str(Box(iv1,iv2).info.var2val) ' ' Box(iv1,iv2).info.var2units];
%  plot the time series
   tMyrs=Box(iv1,iv2).run_info.timesteps_myrs;
   Z = Box(iv1,iv2).run_info.Z_km;
   nt = numel(tMyrs); 
   it_1=1; 
   nts=numel(it_1:t_skip:nt); ip = 1;  
   
   for it = it_1:t_skip:nt   
     F = Box(iv1,iv2).Frames(it); 
     clr=[ip/nts 0 ip/nts]; ip = ip + 1; 
     
%    bg velocity
     subplot(1,4,1)     
     if it > it_1; hold all; end
     plot(-F.Vbgz*3600*24*365*100,Z,'color',clr)
      ylabel('z [km]')
     xlabel('V_{bg} [cm/yr]')
     
%    T,Tsol
     subplot(1,4,2)
     if it == it_1; plot(F.Tsol,Z,'--k'); end
     hold all; 
     plot(F.T,Z,'color',clr)
     hold off; 
     xlabel('T [^oC]');
     xlim([0 1550]) 
     set(gca,'xtick',[0 500 1000 1500])
     
%    phi
     subplot(1,4,3)
     if it > it_1; hold all; end
     plot(F.phi,Z,'color',clr)
     if it > it_1; hold off; end
     xlabel('\phi')
     title(tit2name); 
     xlim([0 0.03])
%    visc
     subplot(1,4,4)
     if it > it_1; hold all; end
     semilogx(F.eta,Z,'color',clr)
     if it > it_1; hold off; end
     xlabel('\eta [Pa s]')
     xlim([1e18 4e26])
     set(gca,'xtick',[10.^(18:26)])
          
           
   end

   for ip=1:4
       subplot(1,4,ip)
       set(gca,'ydir','rev')
       ylim(ylims)
   end
  
   


%% ------------------------------------------------------------------------
%% Done plotting, save the figures
%% ------------------------------------------------------------------------

 for ifig = 1:numel(Fig)
 if strcmp(Fig(ifig).savefig,'yes') 
     figure(Fig(ifig).fig); 
     set(gcf, 'PaperUnits', 'inches');
     set(gcf, 'PaperPosition', [0 0 Fig(ifig).paper_x_width Fig(ifig).paper_y_width]); %
     saveas(gcf,[Fig(ifig).savedir Fig(ifig).savename],'epsc')
 end
 end
 
 