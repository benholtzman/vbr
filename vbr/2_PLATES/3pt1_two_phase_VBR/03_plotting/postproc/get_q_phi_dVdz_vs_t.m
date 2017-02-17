function OutVars = get_q_phi_dVdz_vs_t(Box,iv1,iv2,MELT)
 
nt = numel(Box(iv1,iv2).run_info.timesteps_myrs);
phiFullColumn=zeros(1,nt); phiDBL=zeros(1,nt); phiMax=zeros(1,nt); 
heat=zeros(1,nt);
zSOLt=zeros(1,nt); 
for it = 1:nt
    phi=Box(iv1,iv2).Frames(it).phi; 
    if max(phi)>1e-8;        
        [phiFullColumn(it),phiMax(it),phiDBL(it),dc]=calc_Phi_vals(Box,iv1,iv2,it,'calc_melt');
        zSOLt(it)=Box(iv1,iv2).run_info.zSOL(it)/1e3;
    else
        [phiFullColumn(it),phiMax(it),phiDBL(it),dc]=calc_Phi_vals(Box,iv1,iv2,it,'dont_calc_melt');
    end
    [heat(it)]=calc_heatflow(Box,iv1,iv2,it);
end

OutVars.PHI_Mean=phiFullColumn;
OutVars.PHI_DBL=phiDBL;
OutVars.PHI_Max=phiMax;
OutVars.tMyr = Box(iv1,iv2).run_info.timesteps_myrs;
OutVars.zLAB = Box(iv1,iv2).run_info.zLAB(1:nt)/1e3;
OutVars.zSOL = zSOLt;
OutVars.Maxq = heat;
OutVars.nz=dc.nz; OutVars.dzkm=dc.dzkm;
end

function [phiAVE,phiMax,phiDBL,dc]=calc_Phi_vals(Box,iv1,iv2,it,MELT)

    zLAB = Box(iv1,iv2).run_info.zLAB(it)/1e3;
    sets=Box(iv1,iv2).run_info.settings;
    zkm  = Box(iv1,iv2).run_info.Z_km;
    dzkm=zkm(2)-zkm(1);
    dc.dzkm=dzkm; dc.nz=numel(zkm); 
    
    
    if strcmp(MELT,'calc_melt')
        
%   local settings
    dc.ref_average_width=15; % bottom depth averaging depths
    
%   extract frame   
    FR=Box(iv1,iv2).Frames(it);
    phiz=FR.phi; T=FR.T; Tsol=FR.Tsol; eta=FR.eta;
        
%   get constants
    a=sets.grain0;
    C=sets.C;
    muf=sets.mufo;
    
%   calc compaction length
    xi = eta.*(phiz.^sets.nxi);
    k = a*a*phiz.^sets.np/C;
    dc.dc = sqrt(xi.*k/muf)/1e3;
    
%   find melting onset
    [MOz,MOid]=max(zkm.*(T>Tsol));
    [junk,LABid]=min(abs(zkm-zLAB));
    
    dc.ref_Naves=round(dc.ref_average_width./dzkm);
    
    MO1 = MOid-dc.ref_Naves;
    MO1 = (MO1-1) * (MO1 >= 1) + 1;
    if MO1 == 1;         
        dc.ref_km=mean(dc.dc);
        phiAVE=mean(phiz(phiz>1e-6));
        phiMax=max(phiz(phiz>1e-6));
        phiDBL=mean(phiz(phiz>1e-6));
    else
        dc.ref_km=mean(dc.dc(MO1:MOid));  
        phiAVE=mean(phiz(T>Tsol));
        phiMax=max(phiz(T>Tsol));
        phiDBL=mean(phiz(LABid+1:LABid+1+round(10/dzkm)));
    end
    
     
%   calculate vertical means
    
    else
        phiAVE=0;phiMax=0;phiDBL=0;
    end
end


function [maxq,heatflow]=calc_heatflow(Box,iv1,iv2,it)

%   local settings
    dc.ref_average_width=15; % bottom depth averaging depths
    
%   extract some stuff    
    FR=Box(iv1,iv2).Frames(it);
    T=FR.T; k=FR.Kc;
       
    zkm  = Box(iv1,iv2).run_info.Z_km;
    dz=1000*(zkm(2)-zkm(1));

%   calculate it
    dTdz=(T(2:end)-T(1:end-1))./dz;
    kT=(k(2:end)+k(1:end-1))/2; % on edges
    heatflow=kT.*dTdz;
    maxq=max(abs(heatflow));
end