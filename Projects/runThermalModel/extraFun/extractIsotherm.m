function zIso=extractIsotherm(Box0,isoVal)

  tsteps=Box0.run_info.tMyrs;
  zIso=zeros(size(tsteps));
  z=Box0.run_info.Z_km;
  for it=1:numel(tsteps)
    T=Box0.Frames(it).T;
    zIso(it)=min(z(T>=isoVal));
  end
end
