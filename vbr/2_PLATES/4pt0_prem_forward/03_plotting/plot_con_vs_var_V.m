target_t = 2; 
target_z = 50; 

x_width_Example=6 ;y_width_Example=2;

F(1).file='../../../../../0_BOXES/2016-03-04-conV/individual_runs/2016-03-04-conV_1.mat';
F(2).file='../../../../../0_BOXES/2016-03-04-varV/individual_runs/2016-03-04-varV_1.mat';
F(1).clr=[0 0.6 0]; F(2).clr=[1 0 0]; 
F(1).dn='con'; F(2).dn='var'; 

xlimits = [0 3.5; ...
           0 1700; ...
           0 1700];
ylimits = [35 80; ...
           20 125; ...
           20 125];
       
labs(1).xticks=[0 1 2 3];
labs(2).xticks=[0 500 1000 1500];
labs(3).xticks=[0 500 1000 1500];
       
labs(1).xlab='sqrt(t) [Myrs^{1/2}]';
labs(3).xlab='T [^oC] at z_{LAB} of 40 km';
labs(2).xlab='T [^oC] at 2 Myr';

labs(1).ylab='z_{LAB} [km]';
labs(2).ylab='depth [km]';
labs(3).ylab='depth [km]';




figure('color',[1 1 1])
load(F(1).file)
for ip = 2:3
    subplot(1,3,ip)    
    plot(Vars.Tsol(:,1),Info.z_km,'--k','displayname','sol'); 
    hold on
    plot(Vars.T(:,1),Info.z_km,'k','displayname','T(0)'); 
    Tadi=Vars.T(end,1) - (Info.z_km(end)-Info.z_km)*settings.dTdz_ad*1e3; 
    plot(Tadi,Info.z_km,'--k','displayname','adi')
    hold off
end


for ifi = 1:2

    load(F(ifi).file)
    F(ifi).zLAB=Info.zLAB/1e3;
    F(ifi).t=Info.tMyrs;
    
    [v1 idz]=min(abs(F(ifi).zLAB-target_z));
    [v1 idt]=min(abs(F(ifi).t-target_t));
    
    F(ifi).Tt=Vars.T(:,idt);
    F(ifi).Tz=Vars.T(:,idz);
    F(ifi).Tz_t=Info.tMyrs(idz);
    F(ifi).z=Info.z_km;

    subplot(1,3,1)
    if ifi>1; hold on; end
    plot(sqrt(F(ifi).t),F(ifi).zLAB,'color',F(ifi).clr,'displayname',F(ifi).dn)
    if ifi>1; hold off; end
    
    subplot(1,3,2)
    hold on
    plot(F(ifi).Tt,F(ifi).z,'color',F(ifi).clr,'displayname',F(ifi).dn)
    hold off
    
    subplot(1,3,3)
    hold on
    plot(F(ifi).Tz,F(ifi).z,'color',F(ifi).clr,'displayname',F(ifi).dn)
    hold off
end


for ip=1:3
   subplot(1,3,ip)
   set(gca,'ydir','rev')
   xlim(xlimits(ip,:))
   ylim(ylimits(ip,:))
   xlabel(labs(ip).xlab)
   ylabel(labs(ip).ylab)
   set(gca,'xtick',labs(ip).xticks)
end

subplot(1,3,2)
legend('location','southwest')

 set(gcf, 'PaperUnits', 'inches');
 set(gcf, 'PaperPosition', [0 0 x_width_Example y_width_Example]); %
 saveas(gcf,['/home/chris/Desktop/blah.eps'],'epsc')
