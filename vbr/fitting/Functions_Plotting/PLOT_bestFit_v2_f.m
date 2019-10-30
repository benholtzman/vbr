%function PLOT_bestFit_v2_f(Box,VarInfo,Vel_ExptSet,tsnap,...
 %   depthrange,weighted,obsVs,obsdepth,ylimits,i_Var1_BF,j_Var2_BF)

clf; 

dur = [ num2str(round(Box(1).Movie.info.timesteps_myrs(tsnap))) ' Myrs' ] ; 

          
%% ================================
% PLOTTING  setting
%  ================================

H = 0.6 ; 

% figure size
  figsize_inches = [0.1 0.6 8 6];
% left bottom width height
  plot_A = [0.1 0.3 0.25 H ] ; 
  plot_B = [0.39 0.3 0.15 H ] ;
  plot_C = [0.57 0.3 0.03 H];
  plot_D = [0.7 0.3 0.25 0.28 ] ;  
  plot_E = [0.7 0.78 0.25 0.12 ] ; 
  LBLFNT = 12 ;

%% color setting 
col1 = [0.0 0.4 0.1] ; % Andrade
col2 = [0.0 0.6 0.1] ; % Andrade.COMPOSITE
col3 = [0.0 0.6 0.1] ; % AndradePsP
col4 = [0.0 0.0 0.9] ; % AndradePsP.COMPOSITE
col5 = [0.8 0.0 0.2] ; % eBurgers

heat = [0.9 0.5 0.1] ;



%% start to get the V_VBR from Box and plot Vs(depth)

Z_km = Box(1,1).Movie.info.Z_km ; 

Vel_ExptSet = 3
clear Movie obj_vel 
%Movie=Box(indBox).Movie;
Movie=Box(i_Var1_BF,j_Var2_BF).Movie;
switch Vel_ExptSet
    case 1
        V_VBR = mean(Movie.Frames(tsnap).VBR.Andrade.Va ./1e3);
        Std_V = std(Movie.Frames(tsnap).VBR.Andrade.Va ./1e3);
    case 2    
        V_VBR = mean(Movie.Frames(tsnap).VBR.Andrade.Va_comp ./1e3) ; 
        Std_V = std(Movie.Frames(tsnap).VBR.Andrade.Va_comp ./1e3) ;
    case 3
        V_VBR = mean(Movie.Frames(tsnap).VBR.AndradePsP.Va ./1e3); 
        Std_V = std(Movie.Frames(tsnap).VBR.AndradePsP.Va ./1e3);
    case 4
        V_VBR = mean(Movie.Frames(tsnap).VBR.AndradePsP.Va_comp ./1e3);
        Std_V = std(Movie.Frames(tsnap).VBR.AndradePsP.Va_comp ./1e3);
    case 5
        V_VBR = mean(Movie.Frames(tsnap).VBR.eBurgers.V ./1e3);
        Std_V = std(Movie.Frames(tsnap).VBR.eBurgers.V ./1e3);
end



V_VBR_color = col3;

%% =================
% VELOCITY !
% =================
%figure('color',[1 1 1])
set(gcf, 'Units', 'inches');
set(gcf, 'Position', figsize_inches);
axes('position', plot_A,'box','on'); hold on; 
set(gca,'YDir','reverse') 


%h=fill([meanVs+Std_V fliplr(meanVs-Std_V)],[Z_km fliplr(Z_km)],V_VBR_color,'FaceAlpha', 1);
%h=fill([meanVs+Std_V fliplr(meanVs-Std_V)],[Z_km fliplr(Z_km)],V_VBR_color,'FaceAlpha', 1);

%set(h,'EdgeColor',V_VBR_color);
patchline(V_VBR, Z_km,'linestyle','-','edgecolor',heat,'linewidth',2.0,'edgealpha',1); hold on;

forwardVs= interp1(Z_km,V_VBR,obsdepth);

ind=find(W ~= 0);
plot(forwardVs(ind),obsdepth(ind),'x','color',V_VBR_color,'LineWidth', 2.0); hold on;

% ===========
% plot obs on top ! 
plot(obsVs,obsdepth,'k-','LineWidth', 2);hold on;


xlim([3.5 5.5]);
% ylim([Z_km(1) Z_km(end)]);
ylim(ylimits)
ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
xlabel('Vs (km/s)','fontname','Times New Roman','fontsize',LBLFNT); 


if VarInfo.Var2_n > 1
title_name = ['Best fit: ' VarInfo.Var1_name '=' num2str(VarInfo.Var1_range(i_Var1_BF)) ...
              ' ' VarInfo.Var1_units ', ' VarInfo.Var2_name '=' num2str(VarInfo.Var2_range(j_Var2_BF)) ...
              ' ' VarInfo.Var2_units];
else
title_name = ['Best fit: ' VarInfo.Var1_name '=' num2str(VarInfo.Var1_range(i_Var1_BF)) ...
              ' ' VarInfo.Var1_units];
end
          
title(title_name,...
    'fontname','Times New Roman','fontsize',LBLFNT)

%% ===================
% PLOT difference !
% ===================
axes('position', plot_B,'box','on'); hold on; 
set(gca,'YDir','reverse') 

ind=find(W ~= 0);
plot(100*(obsVs(ind)-forwardVs(ind))./obsVs(ind),obsdepth(ind),'color',V_VBR_color,'LineWidth', 2); hold on;


%xlim([-0.2 0.2]);
plot([0 0], [Z_km(1) Z_km(end)],'k--');
% ylim([Z_km(1) Z_km(end)]);
ylim(ylimits)
set(gca, 'YTick', []);
title('% Difference','fontname','Times New Roman','fontsize',LBLFNT);


%% =================
% PLOT weighted(depth)
% ==================
axes('position', plot_C); hold on; 
set(gca,'YDir','reverse','box','on') 
%plot(W,obsdepth,'.','color',[0.5,0.5,0.5],'LineWidth', 1.0); hold on;
plot(W,obsdepth,'color',[0.2,0.2,0.2],'LineWidth', 3); hold on;

xlim([-0.1 1.2]);
% ylim([Z_km(1) Z_km(end)]);
ylim(ylimits)
set(gca, 'YTick', [],'XTick',[0 1]);
title('W','fontname','Times New Roman','fontsize',LBLFNT)


%% =================
% PLOT misfit !
% =================
axes('position', plot_D,'box','on'); hold on; 


if VarInfo.Var2_n > 1
    % this was the problem: %axes('position', plot_D); hold on;
    % PLOT the matrix of residuals
    x = VarInfo.Var2_range ;
    y = VarInfo.Var1_range ;
    h1 = imagesc(x,y,misfit);%hold on;
    xlabel([VarInfo.Var2_name ' [' VarInfo.Var2_units ']'],'fontname','Times New Roman','fontsize',LBLFNT);
    ylabel([VarInfo.Var1_name ' [' VarInfo.Var1_units ']'],'fontname','Times New Roman','fontsize',LBLFNT);
    title('Misfit','fontname','Times New Roman','fontsize',LBLFNT);
    plot(VarInfo.Var2_range(j_Var2_BF),VarInfo.Var1_range(i_Var1_BF),'or','MarkerSize',16)
    axis tight;
    hold off;
    
    fix_this=0
    if fix_this==1 
    caxis([0 1]);
    % make colorbar, put it on the top (north) and move it higher
    h=colorbar('location','North');
    pos = get(h,'position');
    pos(2) = pos(2) + 0.13; % adjust vertical location (trial and error here...)
    set(h,'position',pos)
    %view(-90,90) 
    end
    
    colormap(gray);

	% ===============================================================
	%%%  Small line plot of min misfit at each Tpot, above 2D plot
    
    for iVar = 1:numel(VarInfo.Var2_range);
       [MinFit(iVar) minTi] = min(misfit(:,iVar));
       MinT(iVar) = VarInfo.Var1_range(minTi);
    end
    E_axis=axes('position', plot_E); hold on;
    plot(VarInfo.Var2_range,MinFit,'k','marker','.')
    xlabel([VarInfo.Var2_name ' (' VarInfo.Var2_units ')'],...
        'fontname','Times New Roman','fontsize',LBLFNT);
    ylabel(['min misfit'])
        
    ax2 = axes('Position',plot_E,'XAxisLocation','top','YAxisLocation',...
     'right','Color','none','xcolor','r','ycolor','r','xticklabel',[]);
    hold on
    plot(VarInfo.Var2_range,MinT,'r','marker','.','parent',ax2)
    
    
elseif VarInfo.Var2_n==1
%    axes('position', plot_D); hold on;
%    plot(VarInfo.Var1_range,misfit);
%    xlabel([VarInfo.Var1_name ' (' VarInfo.Var1_units ')'],...
%        'fontname','Times New Roman','fontsize',LBLFNT);
%    ylabel('Misfit','fontname','Times New Roman','fontsize',LBLFNT);
    
%    axis tight;
    
%    hold on
%    Misfitrange = [0 max(misfit)];
%    Misfitrangex = [VarInfo.Var1_range(i_Var1_BF) VarInfo.Var1_range(i_Var1_BF)];
%    plot(Misfitrangex,Misfitrange,'--k')
%    hold off
end

print(gcf,'-depsc',figName_1)
%end

