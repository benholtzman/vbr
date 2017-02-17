%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script builds the initial settings to call the Driver with. Commonly
% changed parameters are set here, less frequently changed parameters are
% in initialize_settings. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% set sweep parameters
  Vbgtest =20; % [cm/yr]
  Musotest = 20; % log10(muso) [Pa s]
  muftest = 1; % [Pa s]  
%   sourcetest = [0 25 50 100 200 300 400 500]; % source fraction (either PPM H2O or percent pyroxenite)
  sourcetest =    [0 25 50 100 200 300 20 40];
  sourcetest(2,:)=[2 2  2  2   2   2   3  3]; 
% time and mesh settings used for all runs
  settings.nt = 200000; % number of timesteps to run
  settings.outk =2500;%settings.nt/4; % output frequency
  settings.dz = 0.02; % node spacing in reference compaction lengths
  settings.zmin =80; % min depth for model domain [km]  
  settings.zmax = 160; % max z depth [km]    
  settings.Gdot_H = 50;% melt production integration height [km]
  
% source productivity settings
%   settings.CompID=1; % source composition (1 batch peridotie, 2 fractional peridotite, 3 dry perid-pyrox)  
  settings.MeltingOnset = 150; % melting onset depth [km]
  
% save info  
  savedir=['data_realdFdP_Mu_' num2str(Musotest) '_Vup_' num2str(Vbgtest)]
  if exist(savedir,'dir')~=7
     mkdir(savedir); 
  end
  savebase = [savedir '/run'];
  
% set up run ID info
  nruns = numel(Vbgtest)*numel(Musotest)*numel(muftest)*numel(sourcetest(1,:));
  Runs.phiref=zeros(1,nruns);
  Runs.dc=zeros(1,nruns);
  Runs.Vup = zeros(1,nruns);
  Runs.Muso=zeros(1,nruns);
  Runs.SourceFrac=zeros(1,nruns);
  Runs.Muf=zeros(1,nruns);
  Runs.runid=zeros(1,nruns);

% parameter sweep
 irun =1;
for iVbg = 1:numel(Vbgtest)
for imuso = 1:numel(Musotest)
for imuf = 1:numel(muftest)
for iSrc = 1:numel(sourcetest(1,:))
  display(' ')
  display('______________________________________________________________')
  display(' ')
  display(['building settings for run ' num2str(irun) ' of ' num2str(nruns)])

% melt generation and domain settings  
  settings.Vbg = Vbgtest(iVbg); % upwelling velocity [cm/yr]  
  settings.muso = 10^Musotest(imuso); % reference solid shear viscosity [Pa s]
  settings.muf = muftest(imuf); % fluid shear viscosity [Pa s] 
  settings.SourceFrac=sourcetest(1,iSrc); % source water or pyroxenite fraction
  settings.CompID=sourcetest(2,iSrc);
  if exist('find_pot.m','file')==0
      addpath('../../dFdP')
  end
  settings.Tpot=find_pot(settings.CompID,settings.SourceFrac,settings.MeltingOnset,0); % potential temperature
% save options
  savename = [savebase num2str(irun)]; 
  
% initialize the arguments passed to solver  
  [scales refs settings]=initialize_settings(settings);
  display(settings)
  display(refs)
  display(scales)

% solve it!  
  display('starting run')
  Vars=Driver(scales,settings);
  display('run finished')
  if Vars.final == 1
      display('  should be at steady state')
  elseif Vars.final == 0
      display('  might not be at steady state')
      display(['  last residual:' num2str(Vars.dphi)])
  elseif Vars.final == 999
      display('  NOOOOOO! NOT THE NANs!!!')
      display(scales)
  end
  
% save it
  display('saving data') 
  save(savename,'Vars','scales','settings','refs')
  
% plotting
  display('plotting')
%  dnm = ['\mu_s^o=10^{' ...
%       num2str(log10(settings.muso)) '} Pa s'];
  if settings.CompID <= 2
   dnm = [num2str(settings.SourceFrac) ' PPM H_2O, T_p=' num2str(settings.Tpot)];
  elseif settings.CompID == 3
   dnm = ['\phi_{opx}=' num2str(settings.SourceFrac) ', T_p=' num2str(settings.Tpot)];
  end
  Gdot = Vars.Gdot(:,end)/refs.to/(365*3600*24);
  if irun == 1
  [f1 mov]= finalprofile(Vars.z,Gdot,Vars.phi,...
      Vars.P,Vars.t,numel(Vars.t),refs,'\Gamma [1/s]','\phi','P',...
      dnm,0,0,0);
  else
  [f1 mov]= finalprofile(Vars.z,Gdot,Vars.phi,...
      Vars.P,Vars.t,numel(Vars.t),refs,'\Gamma [1/s]','\phi','P',...
      dnm,0,1,f1);
  end
  pause(0.5)
  saveas(f1,[savebase '.eps'],'epsc')
      
  Runs.phiref(irun)=scales.phiref;   
  Runs.Vup(irun) = Vbgtest(iVbg); 
  Runs.dc(irun) = settings.dc; 
  Runs.Muso(irun)=Musotest(imuso);   
  Runs.Muf(irun)=muftest(imuf);
  Runs.SourceFrac(irun)=sourcetest(iSrc);
  Runs.runid(irun) = irun;
  irun = irun+1; 
end
end
end
end
save([savebase '_info'],'Runs')

%  close(f1)
% display('post processing complete')
%  quit force
