function data = loadYT2016visc()
  dataDir='../../../../vbrWork/expt_data/3_attenuation/Yamauchi2016/';
  data=struct();
  if exist([dataDir,'table3.mat'],'file')
    disp('loading')
    load([dataDir,'table3.mat'])
    data.table3_H=table3_H;
  end

  if exist([dataDir,'viscosity_table2subset.csv'],'file')
    d=csvread([dataDir,'viscosity_table2subset.csv']);
    d=d(2:end,:);
    data.visc=struct();
    data.visc.sample=d(:,1);
    data.visc.dg_um=d(:,2);
    data.visc.T_C=d(:,3);
    data.visc.T_C_pm=d(:,4);
    data.visc.eta=d(:,5)*1e12;
    data.visc.sample_list=unique(data.visc.sample);
  end

end
