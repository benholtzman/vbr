function data = loadYT2016Q()
  dataDir='../../../../vbrWork/expt_data/3_attenuation/Yamauchi2016/';
  data=struct();

  if exist([dataDir,'YT16_41_fQinv_allT.csv'],'file')
    d=csvread([dataDir,'YT16_41_fQinv_allT.csv']);
    d=d(2:end,:);
    data.Qinv=struct();
    data.Qinv.sample=d(:,1);
    data.Qinv.T_C=d(:,2);
    data.Qinv.f=d(:,3);
    data.Qinv.Qinv=d(:,4);
    data.Qinv.sample_list=unique(data.Qinv.sample);
  end


  if exist([dataDir,'YT16_41_fE_allT.csv'],'file')
    d=csvread([dataDir,'YT16_41_fE_allT.csv']);
    d=d(2:end,:);
    data.E=struct();
    data.E.sample=d(:,1);
    data.E.T_C=d(:,2);
    data.E.f=d(:,3);
    data.E.E=d(:,4);
    data.E.sample_list=unique(data.E.sample);
  end
end
