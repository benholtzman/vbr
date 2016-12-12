function [F1]=PLOT_bestFit_ch(Box,VarInfo,Fit_Params,Obs,ylimits,xlimits)

i_Var1=Fit_Params.best_Var1;
j_Var2=Fit_Params.best_Var2;
misfit=Fit_Params.misfit;
W=Fit_Params.wtz; 

V_VBR_color =  [0.0 0.6 0.1] ; % AndradePsP
%% =========================
%  get mean V and standard deviation
   
Z_km = Box(i_Var1,j_Var2).run_info.Z_km; 
obsVs=Obs.Vs/1e3;
obsdepth=Obs.depth;
        
iFrame = Fit_Params.Frame_Selection(i_Var1,j_Var2);
VBR = Box(i_Var1,j_Var2).Frames(iFrame).VBR.out.anelastic.(Fit_Params.VBR_anelastic_method);
V_VBR = VBR.Vave/1e3;
        
Std_V = std(VBR.Va./1e3,0,2);

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
  plot_D = [0.68 0.3 0.25 0.25 ] ;  
  plot_E = [0.68 0.75 0.25 0.12 ] ; 
  LBLFNT = 12 ;

heat = [0.9 0.5 0.1] ;


%% =================
% VELOCITY !
% =================
F1=figure('color',[1 1 1]);
set(gcf, 'Units', 'inches');
set(gcf, 'Position', figsize_inches);

axes('position', plot_A,'box','on');  
   
   plot(obsVs,obsdepth,'k-','LineWidth', 1.5);
   hold on;
   set(gca,'YDir','reverse')

   h=fill([(V_VBR+Std_V)' fliplr((V_VBR-Std_V)')],[Z_km' fliplr(Z_km')],...
           V_VBR_color,'FaceAlpha', 1);
   set(h,'EdgeColor',V_VBR_color);
%    plot(V_VBR+2*Std_V,Z_km,'color',[0.2 0.2 0.2])
%    plot(V_VBR-2*Std_V,Z_km,'color',[0.2 0.2 0.2])
   patchline(V_VBR, Z_km,'linestyle','-','edgecolor',heat,'linewidth',1,...
             'edgealpha',1); 

   forwardVs= interp1(Z_km,V_VBR,obsdepth);
   ind=find(W ~= 0);
   plot(forwardVs(ind),obsdepth(ind),'x','color',V_VBR_color,'LineWidth', 0.5); hold on;
   
   xlim([3.5 5.5]);
   ylim(ylimits)
   ylabel('Z (km)','fontname','Times New Roman','fontsize',LBLFNT)
   xlabel('Vs (km/s)','fontname','Times New Roman','fontsize',LBLFNT);
   
   
   if VarInfo.Var2_n > 1
       title_name = ['Best fit: ' VarInfo.Var1_name '=' num2str(VarInfo.Var1_range(i_Var1),'%0.0f') ...
           ' ' VarInfo.Var1_units ', ' VarInfo.Var2_name '=' num2str(VarInfo.Var2_range(j_Var2),'%0.0f') ...
           ' ' VarInfo.Var2_units];
   else
       title_name = ['Best fit: ' VarInfo.Var1_name '=' num2str(VarInfo.Var1_range(i_Var1)) ...
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
    plot(100*(obsVs(ind)-forwardVs(ind))./obsVs(ind),obsdepth(ind),'color',V_VBR_color,'LineWidth', 1.5); hold on;
    
    
    %xlim([-0.2 0.2]);
    plot([0 0], [Z_km(1) Z_km(end)],'k--');
    % ylim([Z_km(1) Z_km(end)]);
    ylim(ylimits)
    set(gca, 'YTick', []);
    title('% Difference','fontname','Times New Roman','fontsize',LBLFNT);


%% =================
% PLOT weighting function
% ==================
axes('position', plot_C); hold on; 
    set(gca,'YDir','reverse','box','on') 
    plot(W,obsdepth,'color',[0.5,0.5,0.5],'LineWidth', 2); hold on;
    
    xlim([-0.1 1.1]);
    ylim(ylimits)
    set(gca, 'YTick', [],'XTick',[0 1]);
    title('W','fontname','Times New Roman','fontsize',LBLFNT)


%% =================
% PLOT misft !
% =================

if VarInfo.Var2_n > 1
    axes('position', plot_D); hold on;
    h1 = imagesc(VarInfo.Var2_range(1:end),VarInfo.Var1_range,log10(misfit(:,1:end)+1e-20));
    xlabel([VarInfo.Var2_name ' (' VarInfo.Var2_units ')'],...
        'fontname','Times New Roman','fontsize',LBLFNT);
    ylabel([VarInfo.Var1_name ' (' VarInfo.Var1_units ')'],...
        'fontname','Times New Roman','fontsize',LBLFNT);
    xlim(xlimits)
    ylim([VarInfo.Var1_range(1) VarInfo.Var1_range(end)]);
    axis tight;

    % make colorbar, put it on the top (north) and move it higher
    h=colorbar('location','North');
    pos = get(h,'position');
    pos(2) = pos(2) + 0.13; % adjust vertical location (trial and error here...)
    set(h,'position',pos)
    colormap(gray);
    title('log_1_0 misfit','fontname','Times New Roman','fontsize',LBLFNT);
    
    hold on
    plot(VarInfo.Var2_range(j_Var2),VarInfo.Var1_range(i_Var1),'or','MarkerSize',16)
    hold off
    
%%%  Small line plot of min misfit at each Tpot, above 2D plot    
    for iVar = 1:numel(VarInfo.Var2_range);
       [MinFit(iVar) minTi] = min(misfit(:,iVar));
       MinT(iVar) = VarInfo.Var1_range(minTi);
    end
    E_axis=axes('position', plot_E); 
      hold on;
      plot(VarInfo.Var2_range,MinFit,'k','marker','.')
      set(E_axis,'Yscale','log')
      xlabel([VarInfo.Var2_name ' (' VarInfo.Var2_units ')'],...
          'fontname','Times New Roman','fontsize',10);
      ylabel(['min misfit'])
      xlim(xlimits)
      
      set(E_axis,'Ytick',10.^[-5 -4 -3 -2 -1 0])
      y0=10.^min(floor(log10(MinFit)));
      y1=10.^max(ceil(log10(MinFit)));
      ylim([y0 y1]);
    
    ax2 = axes('Position',plot_E,'XAxisLocation','top','YAxisLocation',...
     'right','Color','none','xcolor','r','ycolor','r','xticklabel',[]);
    hold on
    plot(VarInfo.Var2_range,MinT,'r','marker','.','parent',ax2)
    xlim(xlimits)
    ylabel([VarInfo.Var1_name ' (' VarInfo.Var1_units ')'],...
          'fontname','Times New Roman','fontsize',10);
      
elseif VarInfo.Var2_n==1
    axes('position', plot_D); hold on;
    plot(VarInfo.Var1_range,misfit);
    xlabel([VarInfo.Var1_name ' (' VarInfo.Var1_units ')'],...
        'fontname','Times New Roman','fontsize',LBLFNT);
    ylabel('Misfit','fontname','Times New Roman','fontsize',LBLFNT);
    
    axis tight;
    
    hold on
    Misfitrange = [0 max(misfit)];
    Misfitrangex = [VarInfo.Var1_range(i_Var1) VarInfo.Var1_range(i_Var1)];
    plot(Misfitrangex,Misfitrange,'--k')
    hold off
end

end

