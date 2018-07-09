function Fit_Params=Calc_bestFit_ch_2(Box,Fit_Params,Obs)

%% set W array
W = ones(length(Obs.depth),1);

for iweighted = 1 :  length(Fit_Params.weighted)
    depthmin =  Fit_Params.depthrange(iweighted,1);
    depthmax =  Fit_Params.depthrange(iweighted,2);
    ind_depth_list = find(Obs.depth >= depthmin & Obs.depth < depthmax );
    W(ind_depth_list) = Fit_Params.weighted(iweighted);
end


%% =================
% VELOCITY !
% =================


[n_Var1, n_Var2] = size(Box);
misfit=zeros(n_Var1,n_Var2);
for i_Var1 = 1:n_Var1
    for j_Var2 = 1: n_Var2

        Frames=Box(i_Var1,j_Var2).Frames;
        Z_km = Box(i_Var1,j_Var2).run_info.Z_km ;

        iFrame = Fit_Params.Frame_Selection(i_Var1,j_Var2);

        VBR = Frames(iFrame).VBR.out.anelastic.(Fit_Params.VBR_anelastic_method);
        meanV = VBR.Vave;

 %      interpolate predicted values to observation depths
        vel  = interp1(Z_km,meanV,Obs.depth);

 %      calculate cumulative misfit
        misfit(i_Var1,j_Var2) = sum( W .* (Obs.Vs - vel).^2 );
        % BH added the  '[]/Obs.Vs' because it should be there, right?
        % and explains why those numbers are so huge ! (March, 2018)
        %misfit(i_Var1,j_Var2) = sum( W .* (Obs.Vs - vel).^2 / Obs.Vs );
        % and got this error:
        % Assignment has more non-singleton rhs dimensions than non-singleton subscripts

    end
end

[i_Var1, j_Var2] = find((misfit-min(min(misfit)))==0);
if numel(j_Var2) > 1;
    i_Var1 = i_Var1(1);
    j_Var2=j_Var2(1);
end

Fit_Params.best_Var1=i_Var1;
Fit_Params.best_Var2=j_Var2;
Fit_Params.misfit=misfit;
Fit_Params.wtz=W;

end
