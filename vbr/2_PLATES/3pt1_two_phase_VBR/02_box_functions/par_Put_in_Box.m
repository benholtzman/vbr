function  [Box] = par_Put_in_Box(Box0,Box)



% stagger routine (output on nodes, not element center)
  stag=@(f) (f(2:end)+f(1:end-1))/2; 

  nvar11=numel(Box0(:,1));
  nvar22=numel(Box0(1,:));
  for ivar1 = 1:nvar11
      for ivar2=1:nvar22
          
          Info=Box0(ivar1,ivar2).Info;
          Vars=Box0(ivar1,ivar2).Vars;
          settings=Box0(ivar1,ivar2).settings;
          
          % get infos!
          nt = numel(Info.t);
          zHigh = stag(Info.z_km);
          nz = numel(zHigh);
          meth = settings.Box.DownSampleMeth;
          DownFactor=settings.Box.DownSampleFactor;
                              
          % store the Box Movie info
          Box(ivar1,ivar2).run_info.timesteps=Info.t;
          Box(ivar1,ivar2).run_info.timesteps_myrs=Info.tMyrs;
          Box(ivar1,ivar2).run_info.ssresid=Info.ssresid;
          Box(ivar1,ivar2).run_info.zLAB=Info.zLABeta;
          Box(ivar1,ivar2).run_info.zSOL=Info.zSOL;
          Box(ivar1,ivar2).run_info.zMO=Info.zMO;
          Box(ivar1,ivar2).run_info.end_of_run=Info.final_message;
          Box(ivar1,ivar2).run_info.Extruded_m=Info.Extruded;
          
          % get downsampled depth bins or points
          if strcmp(meth,'interp')
              zLow = linspace(zHigh(1),zHigh(end),nz/DownFactor)';
          end
          Box(ivar1,ivar2).run_info.Z_km=zLow;
          Box(ivar1,ivar2).run_info.settings = settings;
          
          % loop over time steps, store in Movie Frames
          for it = 1:nt
              
              %  loop over variables in the Vars structure
              Fields = fieldnames(Vars);
              for iFie = 1:numel(Fields);
                  LowRes = downsample(stag(Vars.(Fields{iFie})(:,it)),zHigh,zLow,meth);
                  Box(ivar1,ivar2).Frames(it).(Fields{iFie})=LowRes;
              end
          end
          
      end
  end
end

function Ylow = downsample(Yhigh,Xhigh,Xlow,meth)
 if strcmp(meth,'interp')
     Ylow = interp1(Xhigh,Yhigh,Xlow);
 elseif strcmp(meth,'averaging')
     disp('not implemented!')
 end
end

