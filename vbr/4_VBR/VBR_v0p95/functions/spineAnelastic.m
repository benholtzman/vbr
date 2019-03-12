function [VBR,telapsed]=spineAnelastic(VBR)
  % calculates user-selected anelastic methods
  
  methods_list=VBR.in.anelastic.methods_list; % list of methods to use

  %  Extended Burgers method
  if sum(strncmp('eBurgers',methods_list,8)) > 0
     if isfield(VBR.in.anelastic,'eBurgers')==0
         VBR.in.anelastic.eBurgers=Params_Anelastic('eBurgers');
     end
     if strcmp(VBR.in.anelastic.eBurgers.method,'PointWise')
         telapsed.eBurgers=tic;
         VBR=Q_eBurgers_f(VBR);
         telapsed.eBurgers=toc(telapsed.eBurgers);
     elseif strcmp(VBR.in.anelastic.eBurgers.method,'FastBurger')
         telapsed.eBurgers=tic;
         [VBR]=Q_eFastBurgers(VBR) ;
         telapsed.eBurgers=toc(telapsed.eBurgers);
     end
  end

  %  Andrade Pesudo-Period method
  if sum(strncmp('AndradePsP',methods_list,10)) > 0
     telapsed.AndradePsP=tic;
     if isfield(VBR.in.anelastic,'AndradePsP')==0
        VBR.in.anelastic.AndradePsP=Params_Anelastic('AndradePsP');
     end
     [VBR]=Q_Andrade_PseudoP_f(VBR) ;
     telapsed.AndradePsP=toc(telapsed.AndradePsP);
  end

  % Takei Maxwell Scaling
  if sum(strncmp('YT_maxwell',methods_list,10)) > 0
    telapsed.YT_maxwell=tic;
    if isfield(VBR.in.anelastic,'YT_maxwell')==0
       VBR.in.anelastic.YT_maxwell=Params_Anelastic('YT_maxwell');
    end
    [VBR]=Q_YT_maxwell(VBR) ;
    telapsed.YT_maxwell=toc(telapsed.YT_maxwell);
  end

  % Yamauchi & Takei 2016 - solidus scaling
  if sum(strncmp('YT2016_solidus',methods_list,10)) > 0
    telapsed.YT2016_solidus=tic;
    [VBR]=Q_YT2016_solidus(VBR) ;
    telapsed.YT2016_solidus=toc(telapsed.YT2016_solidus);
  end
end
