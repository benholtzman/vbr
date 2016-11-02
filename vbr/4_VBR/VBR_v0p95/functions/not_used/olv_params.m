% Olivine Fo90 flow law parameters
% from Hirth and Kohlstedt, 2003 (and others after) ======================= 
% c = Coble Diffusion creep (GB)
olv_Fo90.HK03.diff.Ac = 1.5e9 ; % preexponential for coble diffusion creep
olv_Fo90.HK03.diff.Qc = 375e3 ;% activation energy for coble diffusion creep
olv_Fo90.HK03.diff.Vc = 10e-6 ; % activation volume for coble diff creep
olv_Fo90.HK03.diff.pc = 3 ; % grain size exponent
olv_Fo90.HK03.diff.alf_c = 25 ; % melt factor

% d = Dislocation creep 
olv_Fo90.HK03.disl.Ad = 1.1e5 ; % preexponential  
olv_Fo90.HK03.disl.Qd = 530e3 ;% activation energy 
olv_Fo90.HK03.disl.Vd = 15e-6 ; % activation volume 
olv_Fo90.HK03.disl.nd = 3.5 ; % stress exponent
olv_Fo90.HK03.disl.alf_d = 30 ; % melt factor

% g = GBS disl accomodated
% if T<1250+273
olv_Fo90.HK03.gbs.lt1250.Ag = 6500  ; % preexponential for GBS-disl creep
olv_Fo90.HK03.gbs.lt1250.Qg = 400e3  ; % activation energy for GBS-disl creep
olv_Fo90.HK03.gbs.lt1250.Vg = 15e-6  ; % activation volume
% elseif T>=1250+273
olv_Fo90.HK03.gbs.gt1250.Ag = 4.7e10  ; % preexponential for GBS-disl creep
olv_Fo90.HK03.gbs.gt1250.Qg = 600e3  ; % activation energy for GBS-disl creep
olv_Fo90.HK03.gbs.gt1250.Vg = 15e-6  ; % activation volume
% end
olv_Fo90.HK03.gbs.pg = 2 ; % grain size exponent
olv_Fo90.HK03.gbs.ng = 2 ; % stress exponent
olv_Fo90.HK03.gbs.alf_g = 35 ; % melt factor
 
% from Hansen et al, 2012 (LH12)=======================
% c = Coble Diffusion creep (GB)
olv_Fo90.LH12.diff.Ac = 10^7.6 ; % preexponential for coble diffusion creep
olv_Fo90.LH12.diff.Qc = 375e3 ;% activation energy for coble diffusion creep
olv_Fo90.LH12.diff.Vc = 10e-6 ; % activation volume for coble diff creep
olv_Fo90.LH12.diff.pc = 3 ; % grain size exponent
olv_Fo90.LH12.diff.alf_c = 25 ; % melt factor

% d = Dislocation creep 
olv_Fo90.LH12.disl.Ad = 1.1e5 ; % preexponential 
olv_Fo90.LH12.disl.Qd = 530e3 ;% activation energy 
olv_Fo90.LH12.disl.Vd = 15e-6 ; % activation volume 
olv_Fo90.LH12.disl.nd = 3.5 ; % stress exponent
olv_Fo90.LH12.disl.alf_d = 30 ; % melt factor

% g = GBS disl accomodated
olv_Fo90.LH12.gbs.Ag = 10^4.8  ; % preexponential for GBS-disl creep
olv_Fo90.LH12.gbs.Qg = 445e3  ; % activation energy for GBS-disl creep
olv_Fo90.LH12.gbs.Vg = 15e-6  ; % activation volume
olv_Fo90.LH12.gbs.pg = 0.73 ; % grain size exponent
olv_Fo90.LH12.gbs.ng = 2.9 ; % stress exponent
olv_Fo90.LH12.gbs.alf_g = 35 ; % melt factor

