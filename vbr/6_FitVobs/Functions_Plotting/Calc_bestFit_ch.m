function [i_Var1,j_Var2,misfit,W]=Calc_bestFit_ch(Box,Vel_ExptSet,tsnap,...
    depthrange,weighted,obsVs,obsdepth)

% dur = [ num2str(round(Box(1).Movie.info.timesteps_myrs(tsnap))) ' Myrs' ] ;          

%% set W array 
W = ones(length(obsdepth),1);

for iweighted = 1 :  length(weighted)
    depthmin =  depthrange(iweighted,1);
    depthmax =  depthrange(iweighted,2);
    ind_depth_list = find(obsdepth >= depthmin & obsdepth < depthmax );
    W(ind_depth_list) = weighted(iweighted);
end



%% =================
% VELOCITY !
% =================


[n_Var1, n_Var2] = size(Box);

for i_Var1 = 1:n_Var1
    
    for j_Var2 = 1: n_Var2
        clear Movie
        Movie=Box(i_Var1,j_Var2).Movie;      
        Z_km = Movie.info.Z_km ; 
        
        switch Vel_ExptSet 
            case 1  % UNRELAXED (changed Feb 26, 2015, BH) 
                V_VBR = Movie.Frames(tsnap(i_Var1,j_Var2)).VBR.Vsu ./1e3 ;                
				meanV = V_VBR;
            case 2    
                V_VBR = Movie.Frames(tsnap(i_Var1,j_Var2)).VBR.Andrade.Va_comp ./1e3 ; 
                meanV = V_VBR;
            case 3
                V_VBR = Movie.Frames(tsnap(i_Var1,j_Var2)).VBR.AndradePsP.Va ./1e3;                 
				meanV = mean(V_VBR);
            case 4
                V_VBR = Movie.Frames(tsnap(i_Var1,j_Var2)).VBR.AndradePsP.Va_comp ./1e3;
                meanV = V_VBR;
            case 5
                V_VBR = Movie.Frames(tsnap(i_Var1,j_Var2)).VBR.eBurgers.V ./1e3;
                meanV = V_VBR;
        end
    
        vel  = interp1(Z_km,meanV',obsdepth);  % obsdepth is the location of the observations, Z_km is the model z.    
        misfit(i_Var1,j_Var2) = sum( W .* (obsVs - vel).^2 );
        
        
    end
end

[i_Var1, j_Var2] = find((misfit-min(min(misfit)))==0); 
if numel(j_Var2) > 1; % why?
    i_Var1 = i_Var1(1); 
    j_Var2=j_Var2(1); 
end

end

