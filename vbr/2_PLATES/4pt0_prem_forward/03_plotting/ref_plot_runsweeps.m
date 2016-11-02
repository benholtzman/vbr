%   load (or reload) a box?     
    loadnew='yes';  
    recalc='yes'; 
    if strcmp(loadnew,'yes')
        clear;
        loadnew='yes'; recalc = 'yes'; 
    end
    close all

% choose a thing
  cwd=pwd; cd ~; hmdir=pwd; cd(cwd);
  savedir0=[hmdir '/Dropbox/Research/0_Boxes']; % the closet to store the box in 
  savebase = '2016-04-04-Tpot_vs_Vbg'; % boxes end up named ['Box_' savebase]
  addpath ../../../7_Analysis/functions/  

%  set save flag for each figure
   savefigs=[1 1 1 1];
  
% (plot 1) zLAB, zSOL vs t
   Fig(1).box_indeces='both'; % 'both', 'var1', 'var2'
   Fig(1).box_indeces_v1_skip=1; % increment for box index
   Fig(1).box_indeces_v2_skip=1; % increment for box index
   Fig(1).paper_x_width = 3; % figure width for save (inches)
   Fig(1).paper_y_width = 2.5; % figure height for save (inches)   
   Fig(1).savename='FIG_zSOL_t.eps'; % filename (saved in box location)   
   
% (plot 2) time evolution of single run
   Targetv1 = 1550; % variable 1 value
   Targetv2 = 5; % variable 2 value (will find the closest box index)
   t_skip = 2; % time step skip for time evolution plot
   ylims=[0 150];
   x_width_Example=8.5 ;y_width_Example=3.;   
   Fig(2).paper_x_width = 8.5; % figure width for save (inches)
   Fig(2).paper_y_width = 3.0; % figure height for save (inches)   
   Fig(2).savename='FIG_Profiles_t.eps'; % filename (saved in box location)
   
   
% (plot 3) contour plot of thinning amount
   Fig(3).Targett=4;
   Fig(3).ylab='T_{pot} [^oC]';
   Fig(3).xlab='V_{bg} [cm/yr]';
   Fig(3).iv1 = 1;
   Fig(3).iv2 = 2;
   Fig(3).ylimits = [1400 1575];
   Fig(3).xlimits = [.9 10];
   Fig(3).titlename=['\Delta{z_{SOL}} [km] at ' num2str(Fig(3).Targett) ' Myr'];
   Fig(3).paper_x_width=5;
   Fig(3).paper_y_width=3;
   Fig(3).savename=['FIG_dZLAB_' num2str(Fig(3).Targett) '_Myr.eps'];
      
% (plot 4) movie of contour plot. Uses labels and x/y limits from Fig(3)
%  outputs png frames, use ImageMagick on command line to concantenate to 
%  a movie. 
   Fig(4).movietime='popcorn'; % 'popcorn' to run the film. 
   Fig(4).timevec=linspace(1,15,10); % the target times to animate
   Fig(4).paper_x_width=5; 
   Fig(4).paper_y_width=4;
   Fig(4).savename='FIG_dZLAB';
   
%% ------------------------------------------------------------------------  
%% Load the box
%% ------------------------------------------------------------------------  
    if strcmp(loadnew,'yes')
        loadname=['Box_' savebase '.mat'];
        loaddir=[savedir0 '/' savebase '/'];
        addpath ../02_box_functions
        load([loaddir loadname])
    end
    
    Sets = Box(1,1).Movie.info.settings; 
    titlebase = ['V_{bg}=' num2str(Sets.Vbg) ' cm/yr,' ...
                 'T_{pot}^{asth}=' num2str(Sets.Tpot_excess) ' ^oC'];
             
%% ------------------------------------------------------------------------  
%% Postprocessing Calculations  
%% ------------------------------------------------------------------------ 
   if strcmp(recalc,'yes')
       Vals=calc_run_analyses(Box,Fig(3).Targett); % single-valued measures of a run
   end
   
    
    
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
                      
           R = 1-(iv1 - 1) ./ (nv1 - 1); 
           G = 0; 
           B = (iv1 - 1) ./ (nv1 - 1); 
           
           dnm = num2str([Box(iv1,iv2).info.var1val Box(iv1,iv2).info.var2val],2);
                         
           if iplt > 1; hold on; end
           plot(tVars.tMyr,tVars.zSOL,'color',[R G B],'displayname',dnm)
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
 
%% ------------------------------------------------------------------------
%% PLOT 3 - box contour
%% ------------------------------------------------------------------------
 Fig(3).fig=figure('color',[1 1 1]);
  
 if nv2 > 1
     pcolorvar=Vals.zSol_at_Target(Fig(3).iv1:nv1,Fig(3).iv2:nv2); 
     pcolorvar=pcolorvar.* (pcolorvar > 0); 
     imagesc(var2range(Fig(3).iv2:nv2),var1range(Fig(3).iv1:nv1),pcolorvar)     
     colorbar          
 else
    plot(var1range(1:nv1),zSol_at_Target)     
 end
 ylabel(Fig(3).ylab)
 xlabel(Fig(3).xlab)
 ylim(Fig(3).ylimits)
 xlim(Fig(3).xlimits)
 title([Fig(3).titlename ', ' titlebase]); 

 for ifig = 1:3
 if savefigs(ifig) == 1; 
     figure(Fig(ifig).fig); 
     set(gcf, 'PaperUnits', 'inches');
     set(gcf, 'PaperPosition', [0 0 Fig(ifig).paper_x_width Fig(ifig).paper_y_width]); %
     saveas(gcf,[loaddir Fig(ifig).savename],'epsc')
 end
 end
 
%% ------------------------------------------------------------------------
%% PLOT 4 - movie of box contour
%% ------------------------------------------------------------------------
if strcmp(Fig(4).movietime,'popcorn')
 Fig(4).fig=figure('color',[1 1 1]);

%first, calculate frames for movie
if strcmp(recalc,'yes')
    maxzSol = 0;
    for it = 1:numel(Fig(4).timevec)
        Vals=calc_run_analyses(Box,Fig(4).timevec(it)); % single-valued measures of a run
        CumulativeVals(it).Vals=Vals;
        zSols=Vals.zSol_at_Target(Fig(3).iv1:nv1,Fig(3).iv2:nv2);
        maxzSol=max([maxzSol max(max(zSols))]);
        disp(it)
    end
end
% now, run through timevec again (need maxzSol for colorbar clim)
 if exist([loaddir 'movie_files/'],'dir')==0
     mkdir([loaddir 'movie_files/'])
 end
 for it = 1:numel(Fig(4).timevec) 
     Vals=CumulativeVals(it).Vals; 
     pcolorvar=Vals.zSol_at_Target(Fig(3).iv1:nv1,Fig(3).iv2:nv2); 
     pcolorvar=pcolorvar.* (pcolorvar > 0); 
     imagesc(var2range(Fig(3).iv2:nv2),var1range(Fig(3).iv1:nv1),pcolorvar)     
     caxis([0 maxzSol])
     colorbar          
%      colormap('gray')
     ylabel(Fig(3).ylab)
     xlabel(Fig(3).xlab)
     ylim(Fig(3).ylimits)
     xlim(Fig(3).xlimits)
     title([num2str(Fig(4).timevec(it)) ' Myr'])
      
     
     set(Fig(4).fig, 'PaperUnits', 'inches');
     set(Fig(4).fig, 'PaperPosition', [0 0 Fig(4).paper_x_width Fig(4).paper_y_width]); %
     
     pause(0.5)
     figid = num2str(it); 
     ndigs=numel(figid);
     while ndigs < 3
         figid = ['0' figid];
         ndigs=numel(figid);
     end
     MovieBase=[loaddir 'movie_files/' Fig(4).savename]; 
     saveas(Fig(4).fig,[MovieBase '_' figid],'epsc')


 end
 pause(1)
 
%  close(Fig(4).vidobject);
end

%% ------------------------------------------------------------------------
%% Done plotting, save the figures
%% ------------------------------------------------------------------------

 for ifig = 1:3
 if savefigs(ifig) == 1; 
     figure(Fig(ifig).fig); 
     set(gcf, 'PaperUnits', 'inches');
     set(gcf, 'PaperPosition', [0 0 Fig(ifig).paper_x_width Fig(ifig).paper_y_width]); %
     saveas(gcf,[loaddir Fig(ifig).savename],'epsc')
 end
 end
 
 