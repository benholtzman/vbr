function [depth,obsVs] = plot_velprofile(velfile,depth_column,vel_column)
    
    vel=load([velfile,'.txt']);
    depth = vel(:,depth_column);
    obsVs = vel(:,vel_column);
    
    %set(gca,'YDir','reverse');hold on;
    %plot(obsVs,depth,'r:','LineWidth', 2);hold on;
   
end
