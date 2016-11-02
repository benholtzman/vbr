% close all;
Boxid1 = 1; 
Boxid2 = 1; 
VBR = Box(Boxid1,Boxid2).Frames(end).VBR; 
Vars = [VBR.in.SV.T_K-273 ...
        VBR.in.SV.P_GPa ...
        VBR.in.SV.rho ...
        VBR.out.elastic.Gu_0/1e9 ...
        VBR.out.elastic.anharmonic.Gu/1e9 ...
        VBR.out.elastic.anharmonic.Vsu/1e3];

Xlabs = {'T [^oC]'; 'P [GPa]'; '\rho [kg m^{-3}]'; ...
         'Gu_o [GPa]'; 'Gu(T,P) [GPa]'; 'V^s_u [km/s]'};
Z =  Box(Boxid1,Boxid2).run_info.Z_km; 

 figure('color',[1 1 1])
for ip = 1:6
   subplot(2,3,ip)
   
   x = Vars(:,ip); 
   plot(x,Z,'k')
   xlabel(char(Xlabs(ip))); 
   set(gca,'ydir','rev')
   ylim([0 300])
end

subplot(2,3,1)
ylabel('depth [km]')
title('1325^oC T_{pot}^{asth}, 50 km Z_{plate}') 
subplot(2,3,4)
ylabel('depth [km]')