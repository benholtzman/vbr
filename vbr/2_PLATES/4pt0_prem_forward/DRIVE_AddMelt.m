% DRIVE_AddMelt
% read in a box, make a new array of boxes and add melt profiles
% (constructed from error functions)
clear all; 
clf; close all; 
% ===================================
% (1) read in box.. and look at it.. 

testDir = 'y161130_test'
closetPath = strcat('../../../VBRcloset/',testDir,'/')
boxName = strcat('Box_',testDir,'.mat') % Box pre VBR
load(strcat(closetPath,boxName))
display(Box(1,1).info)
display(Box(1,1).run_info)

newBoxName = strcat('Box_',testDir,'_wMelt.mat')
newBoxpath = strcat(closetPath,newBoxName)

szBox = size(Box) ; 
n_var1 = szBox(1) ; 
n_var2 = szBox(2) ; 

 LBLFNT = 12 ;


% ===================================
% MAKE THE NEW BOX !@# 
phi_vec = [0.0:0.002:0.02] ; 
% pick, from the previous results, the row to keep: 
n_var2new = length(phi_vec) ; 

% SET THIS AFTER FIRST FITTING EXPERIMENT ! 
i_var2good = 3 ; 

for i_var1 = 1:n_var1
    for i_var2 = 1:n_var2new
        Box_new(i_var1,i_var2).info = Box(i_var1,i_var2good).info ;
        Box_new(i_var1,i_var2).run_info = Box(i_var1,i_var2good).run_info ;
        Box_new(i_var1,i_var2).Frames = Box(i_var1,i_var2good).Frames ;
        
        Box_new(i_var1,i_var2).info.var2range = phi_vec ; 
        Box_new(i_var1,i_var2).info.var2units = 'melt fraction' ;
        Box_new(i_var1,i_var2).info.var2name = 'phi' ; 
        Box_new(i_var1,i_var2).info.var2val = phi_vec(i_var2);
        
        Z_km = Box_new(i_var1,i_var2).run_info.Z_km ;
        T_z = Box_new(i_var1,i_var2).Frames(end).T(:,1); 
        
        % get zLAB, zSOL, zTcrit? and its index
        % (2) find point below which melt fraction exists: 
        Z_phiBumpCenter_km = Box_new(i_var1,i_var2).run_info.zLAB(end)/1e3 ; 
        ind_phiBumpCenter = find(Z_km > Z_phiBumpCenter_km,1)-1 ; 
        T_phiBumpCenter = T_z(ind_phiBumpCenter) ; 
        
        % add the melt fraction ! 
        N = length(Z_km) ; 
        N_z_steps = floor(N/5); 
        i_mid_step = ind_phiBumpCenter ; 
        [step_vec] = make_meltStep(N,N_z_steps,i_mid_step) ; 
        phi = phi_vec(i_var2) ; 
        phi_step = step_vec.*phi ; 

        Box_new(i_var1,i_var2).Frames(end).phi(:) = phi_step ; 
        
    end
end
 
 
 
% ===================================
% (2) find point below which melt fraction exists: 
%zLAB, zSOL, zMO (what is zMO?)

% TEMPERATURE : 
ax1 = subplot(1,2,1); hold on; 
colors = colormap(hsv(3*n_var1)) ; 
for i_var1 = 1:n_var1
    for i_var2 = 1:n_var2new     
        Z_km = Box_new(i_var1,i_var2).run_info.Z_km ;
        T_z = Box_new(i_var1,i_var2).Frames(end).T(:,1); 
        col = colors(i_var1+round(2*n_var1),:) ; 
        plot(T_z,Z_km, 'color', col); hold on; 
        plot(T_phiBumpCenter,Z_phiBumpCenter_km, 'r.', 'markersize', 10); hold on; 
    end
end
set(gca,'box','on','xminortick','on','yminortick','on','Ydir','rev',...
            'fontname','Times New Roman','fontsize', LBLFNT)
        
% MELT FRACTION PLOTS ! 
ax2 = subplot(1,2,2); hold on; 
colors = colormap(summer(n_var2new+2)) ; 
for i_var1 = 1:n_var1
    for i_var2 = 1:n_var2new     
        Z_km = Box_new(i_var1,i_var2).run_info.Z_km ;
        phi_z = Box_new(i_var1,i_var2).Frames(end).phi(:,1);
        col = colors(i_var2,:) ; 
        plot(phi_z,Z_km, 'color', col); hold on; 
        %plot(T_phiBumpCenter,Z_phiBumpCenter_km, 'r.', 'markersize', 10); hold on; 
    end
end

set(gca,'box','on','xminortick','on','yminortick','on','Ydir','rev',...
            'fontname','Times New Roman','fontsize', LBLFNT)
        

% ===================================
% save the box... (then run VBR) 

clear Box
Box = Box_new ; 
save(newBoxpath) ; 
