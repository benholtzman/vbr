% DRIVE fitting a plate thickness, either from receiver functions or other 
% i.e. Dalton EPSL 2016 based on attenuation. 

clear; close all;


wMelt_flag = 0 ; 
T1orS2_flag = 1 ;

if T1orS2_flag == 1
    testDir = 'y161210_TNA_fit'
    z_LAB_OBS_km = 75.0
elseif T1orS2_flag == 2
    testDir = 'y161210_SNA_fit'
    z_LAB_OBS_km = 200.0 
end


closetPath = strcat('../../../VBRcloset/',testDir,'/')
boxName = strcat('Box_',testDir,'.mat') % Box pre VBR
load(strcat(closetPath,boxName))
display(Box(1,1).info)
display(Box(1,1).run_info)

szBox = size(Box) ; 
n_var1 = szBox(1) ; 
n_var2 = szBox(2) ; 

 LBLFNT = 12 ;

if T1orS2_flag == 1
    z_LAB_OBS_km = 70.0
elseif T1orS2_flag == 2
    z_LAB_OBS_km = 200.0 
end


% ==================================================
% loop over Box and calculate residuals 
% between LAB constraint and prediction
% ==================================================
for i_var1 = 1:n_var1
    for i_var2 = 1:n_var2
        
        Z_km = Box(i_var1,i_var2).run_info.Z_km ;
        T_z = Box(i_var1,i_var2).Frames(end).T(:,1); 
        
        Z_LAB_pred_km = Box(i_var1,i_var2).run_info.zLAB(end)/1e3 ; 
        ind_zLAB = find(Z_km > Z_LAB_pred_km,1)-1 ; 
        T_LAB = T_z(ind_zLAB) ; 
        
        Res = ((Z_LAB_pred_km - z_LAB_OBS_km)^2)/z_LAB_OBS_km ;
        
        zLAB_mat(i_var1,i_var2).Z_km = Z_LAB_pred_km ;
        zLAB_mat(i_var1,i_var2).ind = ind_zLAB ; 
        zLAB_mat(i_var1,i_var2).T_LAB = T_LAB ; 
        Res_mat(i_var1,i_var2) = log10(Res) ; 
    end
end


% ==================================
% PLOTTING
% ==================================

LBLFNT = 15 ;
LW = 2 ;

column1 = [0.1 0.1 0.2 0.7] ;
column2 = [0.35 0.1 0.2 0.7] ;
plot1 = [0.65 0.3 0.3 0.35] ;
%plot2 = [0.65 0.1 0.3 0.3] ;

% COLUMN =========================================================
axes('Position', column1); 

colors = colormap(hsv(3*n_var1)) ; 
for i_var1 = 1:n_var1
    for i_var2 = 1:n_var2
        T_vec_C = Box(i_var1,i_var2).Frames(end).T(:) ; 
        Z_km = Box(i_var1,i_var2).run_info.Z_km(:) ; 
        col = colors(i_var1+round(2*n_var1),:) ;
        plot(T_vec_C,Z_km, 'color', col); hold on;
        
        T_lab = zLAB_mat(i_var1,i_var2).T_LAB ; 
        Z_lab = zLAB_mat(i_var1,i_var2).Z_km ;
        plot(T_lab,Z_lab, 'r.', 'color', col, 'markersize', 12 ); hold on; 
    end
end

T_min = 0 ;
T_max = 1800 ;
depth_var = 15 ; 
% lower left corner x,y , width height
w = T_max - T_min ; 
h = 2*depth_var ; 
rectangle('position',[T_min,z_LAB_OBS_km-depth_var,w,h]); 
line([T_min,T_max],[z_LAB_OBS_km,z_LAB_OBS_km]); 
xlim([T_min,T_max])

xlabel('Temperature [C]')
ylabel('Depth [km]')
set(gca,'box','on','xminortick','on','yminortick','on','Ydir','rev',...
            'fontname','Times New Roman','fontsize', LBLFNT)


% COLUMN =========================================================
axes('Position', column2); 


colors = colormap(hsv(3*n_var1)) ; 
for i_var1 = 1:n_var1
    for i_var2 = 1:n_var2
        log10_eta = log10(Box(i_var1,i_var2).Frames(end).eta(:)) ; 
        Z_km = Box(i_var1,i_var2).run_info.Z_km(:) ; 
        col = colors(i_var1+round(2*n_var1),:) ;
        plot(log10_eta,Z_km, 'color', col); hold on;
        
        eta_lab = log10_eta(zLAB_mat(i_var1,i_var2).ind) ; 
        Z_lab = zLAB_mat(i_var1,i_var2).Z_km ;
        plot(eta_lab,Z_lab, 'r.', 'color', col, 'markersize', 12 ); hold on; 
    end
end

eta10_min = 18 ;
eta10_max = 26 ;
depth_var = 15 ; 
% lower left corner x,y , width height
w = eta10_max - eta10_min ; 
h = 2*depth_var ; 
rectangle('position',[eta10_min,z_LAB_OBS_km-depth_var,w,h]); 
line([eta10_min,eta10_max],[z_LAB_OBS_km,z_LAB_OBS_km]); 
xlim([eta10_min,eta10_max])

xlabel('log_{10} \eta, viscosity [Pa.s]')
set(gca,'box','on','xminortick','on','yminortick','on','Ydir','rev',...
            'fontname','Times New Roman','fontsize', LBLFNT)


% COLUMN =========================================================
axes('Position', plot1); 


colormap('gray'); 
var1range = Box(1,1).info.var1range ; 
var2range = Box(1,1).info.var2range ; 

[var1_mesh,var2_mesh] = meshgrid(var1range,var2range) ; 
 %   h1 = imagesc(VarInfo.Var2_range(3:end),VarInfo.Var1_range,log10(misfit(:,3:end)+1e-20));
imagesc(var1range,var2range,Res_mat') ; 
%surf(var1_mesh,var2_mesh,Res_mat'); 

title('Residual $(\log_{10}[(pred-obs)^2/obs])$', 'interpreter', 'latex') ; 
xlabel(Box(1,1).info.var1name) ;
ylabel(Box(1,1).info.var2name) ;
set(gca,'box','on','xminortick','on','yminortick','on',...
            'fontname','Times New Roman','fontsize', LBLFNT)
