

figure('color',[1 1 1]);  

SP(1).variable ='T'; 
SP(1).variable_2 ='Tsol'; 
SP(2).variable ='Gdot'; 
SP(3).variable ='phi'; 

SP(1).xlab='T [^oC]'; 
SP(2).xlab='\Gamma [s^{-1}]'; 
SP(3).xlab='\phi'; 

SP(1).ylims=[30 140]; 
SP(2).ylims=SP(1).ylims;
SP(3).ylims=SP(1).ylims; 


clr_fade=[0.6 0.6 0.6]; 
clr_new = [1 0 0]; 

z = Box(1,1).Movie.info.Z_km; 

id1=1; 
id2=1; 
nt = numel(Box(id1,id2).Movie.Frames); 
for ip = 3:nt
    
     scl = (ip-1)/(nt-1); 
     clr_new = [scl 0 0]; 
     clr_fade_2=[scl 0.3 0.3]; 
     
     scl = (ip - 2)/(nt-1); 
     clr_fade=[scl 0.8 0.8]; 
     
     
     
     for spid = 1:3
     subplot(1,3,spid)
       hold on 
       old_step=Box(id1,id2).Movie.Frames(ip-1).(SP(spid).variable);
       old_step_2=Box(id1,id2).Movie.Frames(ip-2).(SP(spid).variable);
       new_step=Box(id1,id2).Movie.Frames(ip).(SP(spid).variable);
       plot(old_step_2,z,'color',clr_fade);
       plot(old_step,z,'color',clr_fade_2);       
       plot(new_step,z,'color',clr_new);
       
       
       if isempty(SP(spid).variable_2)==0
          extra=Box(id1,id2).Movie.Frames(ip).(SP(spid).variable_2);
          plot(extra,z,'--k')
       end
       
       
       xlabel(SP(spid).xlab)
       ylabel('z [km]')
       ylim(SP(spid).ylims);
       set(gca,'ydir','rev','box','on')
     end    
     
     pause(0.1); 
  
 end