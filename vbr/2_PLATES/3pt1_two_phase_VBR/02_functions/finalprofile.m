function [f1 mov]= finalprofile(z,xV1,xV2,xV3,t,it,refs,xn1,xn2,xn3,dnm,mflag,oplot,oldfig) 
  dc = refs.dc;
  z = z*dc;
  ylimits = [(min(z)-.5) (max(z)+.5)];
  xl1 = [-.01 (max(max(xV1))+.01)];
  
  yname = 'z [km]';
  
  if mflag == 1
   mov=plotnplay(z,xV1,t*refs.to,xn1,yname,xl1,ylimits,1,-1);
  else
   mov = 0; 
  end
  
  if oplot==0
   f1 = figure('color',[1 1 1]); 
  else
   f1 = oldfig; 
   figure(oldfig)
  end
  
  for iVar=1:3
      
      if iVar == 1
          if (it<=numel(xV1(1,:)))
              xV = xV1(:,it);
          else
              xV = xV1(:,end);
          end
          xn = xn1;
      elseif iVar == 2
          if (it<=numel(xV2(1,:)))
              xV = xV2(:,it);
          else
              xV = xV2(:,end);
          end
          xn = xn2;
      elseif iVar == 3
          if (it<=numel(xV3(1,:)))
              xV = xV3(:,it);
          else
              xV = xV3(:,end);
          end
          xn = xn3;
      end
      

      
      subplot(1,3,iVar)
      if oplot==1
          hold all
      end
      plot(xV,z,'displayname',dnm)
      xlabel(xn);ylabel(yname)
      title(['t =' num2str(t(it)*refs.to/1000) ' [kyrs]'])
      set(gca,'Ydir','rev')
      
      if oplot==1
          hold off      
      end
  end
end
