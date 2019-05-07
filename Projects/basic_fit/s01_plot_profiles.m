%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a basic fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Vs_file='/path/to/file';
load(Vs_file)

idepth=200;
Depth=Vs_Model.Depth(idepth);
disp(Depth)

lons=Vs_Model.Longitude;
lats=Vs_Model.Latitude;
Vs2d=squeeze(Vs_Model.Vs(:,:,idepth))';

figure

imagesc(lons,lats,Vs2d)
set(gca,'Ydir','normal')
xlabel('lon')
ylabel('lat')
title(['Vs at ',num2str(Depth),' km, Shen and Ritzwoller 2016'])
colormap(hot)
colorbar


[val,idepth]=min(abs(Vs_Model.Depth - 80));

% YELLOWSTONE
Lon=360-110.5885;
Lat=44.4280;
[val,iLon]=min(abs(lons-Lon));
[val,iLat]=min(abs(lats-Lat));
VsY=squeeze(Vs_Model.Vs(iLon,iLat,:));
VsY_depth=VsY(idepth);

% Twin Falls Volc. Field (SRP track, 12-6.5 Ma)
Lon=360-114.4701;
Lat=42.5558;
[val,iLon]=min(abs(lons-Lon));
[val,iLat]=min(abs(lats-Lat));
Vs2=squeeze(Vs_Model.Vs(iLon,iLat,:));
Vs2_depth=Vs2(idepth);


% 44.3833° N, 107.1667° W big horn mountains (wyoming craton)
Lon=360-107.1667;
Lat=44.3833;
[val,iLon]=min(abs(lons-Lon));
[val,iLat]=min(abs(lats-Lat));
Vs3=squeeze(Vs_Model.Vs(iLon,iLat,:));
Vs3_depth=Vs3(idepth);

figure
plot(VsY,Vs_Model.Depth,'k','displayname','Yellowstone','LineWidth',2)
hold on
plot(Vs2,Vs_Model.Depth,'r','displayname','Twin Falls VF','LineWidth',2)
plot(Vs3,Vs_Model.Depth,'--k','displayname','Wyoming Craton','LineWidth',2)
set(gca,'Ydir','rev');
xlabel('Vs')
ylabel('Depth')

legend('Location','SouthWest')
xlim([3,5])
ylim([0,150])

disp(['(YS, TF) at ~80 km:(',num2str(VsY_depth),',',num2str(Vs2_depth),')'])

% (YS, TF) at ~80 km:(4.0598,3.9878)
