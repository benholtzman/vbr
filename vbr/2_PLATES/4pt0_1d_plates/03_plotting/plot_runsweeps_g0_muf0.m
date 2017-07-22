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
  savebase = '2016-04-21-g0_vs_muf_con'; % boxes end up named ['Box_' savebase]
  addpath ../../../7_Analysis/functions/  

%  set save flag for each figure
   savefigs=[1 1 1 0];
  
% (plot 1) zLAB, zSOL vs t
   Fig(1).box_indeces='both'; % 'both', 'var1', 'var2'
   Fig(1).box_indeces_v1_skip=1; % increment for box index
   Fig(1).box_indeces_v2_skip=1; % increment for box index
   Fig(1).paper_x_width = 3; % figure width for save (inches)
   Fig(1).paper_y_width = 2.5; % figure height for save (inches)   
   Fig(1).savename='FIG_zSOL_t.eps'; % filename (saved in box location)   
   
% (plot 2) time evolution of single run
   Targetv1 = .01; % variable 1 value
   Targetv2 = 1; % variable 2 value (will find the closest box index)
   t_skip = 2; % time step skip for time evolution plot
   ylims=[0 150];
   x_width_Example=8.5 ;y_width_Example=3.;   
   Fig(2).paper_x_width = 8.5; % figure width for save (inches)
   Fig(2).paper_y_width = 3.0; % figure height for save (inches)   
   Fig(2).savename='FIG_Profiles_t.eps'; % filename (saved in box location)
   
   
% (plot 3) contour plot of thinning amount
   Fig(3).Targett=2;
   Fig(3).ylab='g [m]';
   Fig(3).xlab='\eta_f [Pa s]';
   Fig(3).iv1 = 1;
   Fig(3).iv2 = 1;
   Fig(3).ylimits = [1e-3 .03];
   Fig(3).xlimits = [.01 10];
   Fig(3).titlename=['mean \phi at ' num2str(Fig(3).Targett) ' Myr'];
   Fig(3).paper_x_width=5;
   Fig(3).paper_y_width=3;
   Fig(3).savename=['FIG_phimean_' num2str(Fig(3).Targett) '_Myr.eps'];
      
% (plot 4) movie of contour plot. Uses labels and x/y limits from Fig(3)
%  outputs png frames, use ImageMagick on command line to concantenate to 
%  a movie. 
   Fig(4).movietime='nopecorn'; % 'popcorn' to run the film. 
   Fig(4).timevec=linspace(1,15,10); % the target times to animate
   Fig(4).paper_x_width=5; 
   Fig(4).paper_y_width=4;
   Fig(4).savename='FIG_dZLAB';
   
%% ------------------------------------------------------------------------  
%% Load the box
%% ------------------------------------------------------------------------  
if strcmp(loadnew,'yes')
    addpath ../02_box_functions
    loadname=['Box_' savebase '.mat'];
    loaddir=[savedir0 '/' savebase '/'];
    load([loaddir loadname])
    loadname=['Box_' savebase '_VBR.mat'];
    load([loaddir loadname])
end

Sets = Box(1,1).Movie.info.settings;
titlebase = '';%['V_{bg}=' num2str(Sets.Vbg) ' cm/yr,' ...
%'T_{pot}^{asth}=' num2str(Sets.Tpot_excess) ' ^oC'];

%% ------------------------------------------------------------------------
%% Postprocessing Calculations
%% ------------------------------------------------------------------------
nv1 = numel(Box(:,1));
nv2 = numel(Box(1,:));
if strcmp(recalc,'yes')
    Vals=calc_run_analyses(Box,Fig(3).Targett); % single-valued measures of a run
    
    period_s=10;
    for iv1 = 1:nv1
        for iv2 = 1:nv2
            
            for it = 1:numel(Box(iv1,iv2).Movie.info.timesteps);
                Box(iv1,iv2).Movie.Frames(it).VBR.AndradePsP.Va= ...
                    VBRBox(iv1,iv2).Movie.Frames(it).VBR.Vave_AndradePsP;
            end
            
            [Box(iv1,iv2).dVdz_t,Box(iv1,iv2).melt,Box(iv1,iv2).Vlab]=calc_dVdz_t(Box,iv1,iv2,period_s);
             tMyr = Box(iv1,iv2).Movie.info.timesteps_myrs;
            [val,id]=min(abs(tMyr-Fig(3).Targett)); 
            dVdz(iv1,iv2)=Box(iv1,iv2).dVdz_t(id); 
            VLAB(iv1,iv2)=Box(iv1,iv2).Vlab(id); 
        end
    end
    
end

    
    
%% ------------------------------------------------------------------------
%% PLOT 1 - zLAB vs t 
%% ------------------------------------------------------------------------   
 
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
   
   subplot(1,3,1)
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
   
%    legend('location','southeast')
   set(gca,'ydir','rev')
   xlabel('t [Myrs]')
   ylabel('z_{SOL} [km]')    

%  now find the min  
   var1range=Box(1,1).info.var1range;  
   var2range=Box(1,1).info.var2range;  
 
%% ------------------------------------------------------------------------
%% PLOT 3 - box contour
%% ------------------------------------------------------------------------
 figure(Fig(1).fig);

 pcolorvars(1).var=Vals.PHI_Max(Fig(3).iv1:nv1,Fig(3).iv2:nv2); 
 pcolorvars(2).var=Vals.PHI_DBL(Fig(3).iv1:nv1,Fig(3).iv2:nv2); 
 pcolorvars(1).titname='Max \phi';
 pcolorvars(2).titname='Mean \phi in DBL';
 minmin=min([min(min(pcolorvars(1).var)) min(min(pcolorvars(2).var))]);
 maxmax=max([max(max(pcolorvars(1).var)) max(max(pcolorvars(2).var))]);
 
 for ipcol =1:2
    subplot(1,3,1+ipcol) 
    pcolorvar=pcolorvars(ipcol).var;
    pcolorvar=pcolorvar.* (pcolorvar > 0);
    imagesc(var2range(Fig(3).iv2:nv2),var1range(Fig(3).iv1:nv1),pcolorvar)
    
    caxis([minmin maxmax])
    ylabel(Fig(3).ylab)
    xlabel(Fig(3).xlab)
    ylim(Fig(3).ylimits)
    xlim(Fig(3).xlimits)
    title(pcolorvars(ipcol).titname);
 end
 colorbar
 
%  for ifig = 1:3
%  if savefigs(ifig) == 1; 
%      figure(Fig(ifig).fig); 
%      set(gcf, 'PaperUnits', 'inches');
%      set(gcf, 'PaperPosition', [0 0 Fig(ifig).paper_x_width Fig(ifig).paper_y_width]); %
%      saveas(gcf,[loaddir Fig(ifig).savename],'epsc')
%  end
%  end
 





% Box(iv1,iv2).dVdz_t,Box(iv1,iv2).melt,Fig(3).Targett
figure('color',[1 1 1])

 pcolorvars(1).var=(VLAB); 
 pcolorvars(1).titname='V_s about LAB [km/s]';
 pcolorvars(2).var=abs(dVdz); 
 pcolorvars(2).titname='dV/dz about LAB [1/s]';

 for ipcol =1:2
    subplot(1,2,ipcol) 
    pcolorvar=pcolorvars(ipcol).var;
    pcolorvar=pcolorvar.* (pcolorvar > 0);
    imagesc(var2range(Fig(3).iv2:nv2),var1range(Fig(3).iv1:nv1),pcolorvar)
    
    ylabel(Fig(3).ylab)
    xlabel(Fig(3).xlab)
    ylim(Fig(3).ylimits)
    xlim(Fig(3).xlimits)
    title(pcolorvars(ipcol).titname);
    colorbar
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
        zSols=Vals.PHI_DBL(Fig(3).iv1:nv1,Fig(3).iv2:nv2);
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
     pcolorvar=Vals.PHI_DBL(Fig(3).iv1:nv1,Fig(3).iv2:nv2); 
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

figure
iv1 = 5; 
for iv2=1:6; 
phi=Box(iv1,iv2).Movie.Frames(5).phi.*(Box(iv1,iv2).Movie.Frames(5).T>Box(iv1,iv2).Movie.Frames(5).Tsol);
phi = 10.^(log10(phi)); 
z=Box(iv1,iv2).Movie.info.Z_km; 
hold on 
hold all   
plot(phi,z,'marker','.','displayname',num2str(Box(iv1,iv2).info.var2val))
end



%% ------------------------------------------------------------------------
%% Done plotting, save the figures
%% ------------------------------------------------------------------------

%  for ifig = 1:3
%  if savefigs(ifig) == 1; 
%      figure(Fig(ifig).fig); 
%      set(gcf, 'PaperUnits', 'inches');
%      set(gcf, 'PaperPosition', [0 0 Fig(ifig).paper_x_width Fig(ifig).paper_y_width]); %
%      saveas(gcf,[loaddir Fig(ifig).savename],'epsc')
%  end
%  end
%  
%  