function [G,dGdT,dGdT_ave]= YT16_E(T_C)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % YT2016 Figure 6 caption polynomial fit of G vs T in C
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  a0=2.5943;
  a1=2.6911*1e-3;
  a2=-2.9636*1e-4;
  a3 =1.4932*1e-5;
  a4 =-2.93511e-7;
  a5 =1.8997*1e-9;
  a = [a5,a4,a3,a2,a1,a0];


  G=zeros(size(T_C));
  dGdT=zeros(size(T_C));
  for iTC=1:numel(T_C)
    Gpoly = polyval(a,T_C(iTC));
    dGdTpoly = polyval(a(1:end-1),T_C(iTC));
    G(iTC)=sum(Gpoly(:));
    dGdT(iTC)=sum(dGdTpoly(:));
  end
  dGdT_ave=mean(dGdT);

end
