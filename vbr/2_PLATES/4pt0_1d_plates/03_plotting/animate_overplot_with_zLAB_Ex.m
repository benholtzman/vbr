
saveframes = 0; 
MovieFolder='Movie_fixed';
if saveframes == 1
    if exist(MovieFolder,'dir')==0
        mkdir(MovieFolder)
    else
        delete([MovieFolder '/*'])
    end
end
% close all

figure('color',[1 1 1]);  

SP(1).variable ='T'; 
SP(1).variable_2 ='Tsol'; 
SP(2).variable ='Kc'; 
SP(3).variable ='phi'; 

SP(1).xlab='T [^oC]'; 
SP(2).xlab='K_c [W/m/K]'; 
SP(3).xlab='\phi'; 

SP(1).ylab='z [km]'; 
SP(2).ylab=' '; 
SP(3).ylab=' '; 

SP(1).xscl = 'linear';
SP(2).xscl = 'linear';
SP(3).xscl = 'log';

SP(1).ylims=[0 10*ceil(Info.z_km(end)/10)]; 
SP(2).ylims=SP(1).ylims;
SP(3).ylims=SP(1).ylims; 

SP(1).xlims=[0 1600]; 
SP(2).xlims=[0 6];
SP(3).xlims=[0.000001 0.01]; 

clr_fade=[0.6 0.6 0.6]; 
clr_new = [1 0 0]; 

z = Box(1,1).Movie.info.Z_km; 

id1=1; 
id2=1; 
nt = numel(Box(id1,id2).Movie.Frames); 
tskip = 1; 
tid_0=tskip*2+1; 
Frame_id = 1; 
for ip = tid_0:tskip:nt
    
     scl = (ip-1)/(nt-1); 
     clr_new = [0 0 0]; 
     clr_fade_2=[scl 1-scl 0]; 
     
     scl = (ip - tskip)/(nt-1); 
     clr_fade=[scl 1-scl 0]; 
     
     
     
     for spid = 1:3
     subplot(2,4,[(spid) spid+4])
       hold on 
       old_step=Box(id1,id2).Movie.Frames(ip-tskip).(SP(spid).variable);
       old_step_2=Box(id1,id2).Movie.Frames(ip-tskip*2).(SP(spid).variable);
       new_step=Box(id1,id2).Movie.Frames(ip).(SP(spid).variable);
       
       if ip > tskip*3
           delete(SP(spid).p3)
           delete(SP(spid).p2)
       end
           
       SP(spid).p1=plot(old_step_2,z,'color',clr_fade);
       SP(spid).p2=plot(old_step,z,'color',clr_fade_2);       
       SP(spid).p3=plot(new_step,z,'color',clr_new,'linewidth',2);
       
       
       if isempty(SP(spid).variable_2)==0
          extra=Box(id1,id2).Movie.Frames(ip).(SP(spid).variable_2);
          plot(extra,z,'--k')
       end
       
       
       xlabel(SP(spid).xlab)
       ylabel(SP(spid).ylab)
       ylim(SP(spid).ylims);
       xlim(SP(spid).xlims); 
       set(gca,'ydir','rev','box','on','xscale',SP(spid).xscl)
     end
     
     Info = Box(id1,id2).Movie.info; 
     Info.tMyrs = Info.timesteps_myrs; 
     it=ip; 
     subplot(2,4,4)
     hold on 
     plot(Info.tMyrs(it-tskip*2),Info.zSOL(it-tskip*2)/1e3,'marker','.','markersize',18,'color',clr_fade)
     plot(Info.tMyrs(it-tskip),Info.zSOL(it-tskip)/1e3,'marker','.','markersize',18,'color',clr_fade_2)
     plot(Info.tMyrs(it),Info.zSOL(it)/1e3,'marker','.','markersize',18,'color',clr_new)     
     xlabel('t [Myrs]'); ylabel('z_{SOL} [km]'); set(gca,'ydir','rev','box','on')
         
     
     subplot(2,4,8)
     hold on 
     plot(Info.tMyrs(it-tskip*2),Info.zLAB(it-tskip*2)/1e3,'marker','.','markersize',18,'color',clr_fade)
     plot(Info.tMyrs(it-tskip),Info.zLAB(it-tskip)/1e3,'marker','.','markersize',18,'color',clr_fade_2)
     plot(Info.tMyrs(it),Info.zLAB(it)/1e3,'marker','.','markersize',18,'color',clr_new)     
     xlabel('t [Myrs]'); ylabel('z_{LAB} [km]'); set(gca,'box','on')
     pause(0.2); 
  
     paper_y_width = 4; 
     paper_x_width = 8; 
     set(gcf, 'PaperUnits', 'inches');     
     set(gcf, 'PaperPosition', [0 0 paper_x_width paper_y_width]); %  
     
     if saveframes == 1
     zerobuff = '_';     
     maxzeros = floor(log10(nt))+1; 
     ibuff0=floor(log10(Frame_id))+1; 
     for ibuff = ibuff0: maxzeros
         zerobuff = [zerobuff '0']; 
     end
     saveas(gcf,[MovieFolder '/Frame' zerobuff num2str(Frame_id)],'epsc')
     Frame_id =Frame_id + 1; 
     end
 end