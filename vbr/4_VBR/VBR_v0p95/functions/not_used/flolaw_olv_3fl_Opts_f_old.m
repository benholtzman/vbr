function[sr_c,sr_d,sr_g,sr_tot] = flolaw_olv_3fl_Opts_f(T,P,sig_MPa,d,phi,fH2O,Opt_FloLaw, phi_c, speedfac)
% olivine diffusion creep flow law

olv_params 

% T in Kelvins
% P in Pa (input in MPa)
% sig in MPa
% d in microns
P = P.*1e6; % convert to Pa

% combined GBS Disl FLOW LAW + "diffusion" creep
% diffusion creep + gbs (H&K 2003)
%clear R
R = 8.314 ;

if Opt_FloLaw==0  ;
% Coble Diffusion creep (GB)
Ac = olv_Fo90.HK03.diff.Ac ; % preexponential for coble diffusion creep
Qc = olv_Fo90.HK03.diff.Qc ;% activation energy for coble diffusion creep
Vc = olv_Fo90.HK03.diff.Vc ; % activation volume for coble diff creep
pc = olv_Fo90.HK03.diff.pc ; % grain size exponent
alf_c = olv_Fo90.HK03.diff.alf_c ; 

% Dislocation creep 
Ad = olv_Fo90.HK03.disl.Ad ; % preexponential 
Qd = olv_Fo90.HK03.disl.Qd ;% activation energy 
Vd = olv_Fo90.HK03.disl.Vd ; % activation volume 
nd = olv_Fo90.HK03.disl.nd ; % stress exponent
alf_d = olv_Fo90.LH12.disl.alf_d ; 

% Dislocation-accommodated grain boundary sliding (GBS)
Ag = olv_Fo90.HK03.gbs.lt1250.Ag ; % preexponential for coble diffusion creep
Qg = olv_Fo90.HK03.gbs.lt1250.Qg ;% activation energy for coble diffusion creep
Vg = olv_Fo90.HK03.gbs.lt1250.Vg ; % activation volume for coble diff creep
alf_g = olv_Fo90.HK03.gbs.alf_g ; % melt fraction sensitivity
ng = olv_Fo90.HK03.gbs.ng  ; % stress exponent
pg = olv_Fo90.HK03.gbs.pg  ; % stress exponent


elseif Opt_FloLaw==1|Opt_FloLaw==2 ;
% Coble Diffusion creep (GB)
Ac = olv_Fo90.LH12.diff.Ac ; % preexponential for coble diffusion creep
Qc = olv_Fo90.LH12.diff.Qc ;% activation energy for coble diffusion creep
Vc = olv_Fo90.LH12.diff.Vc ; % activation volume for coble diff creep
pc = olv_Fo90.LH12.diff.pc ; % grain size exponent
alf_c = olv_Fo90.LH12.diff.alf_c ; 

% Dislocation creep 
Ad = olv_Fo90.LH12.disl.Ad ; % preexponential 
Qd = olv_Fo90.LH12.disl.Qd ;% activation energy 
Vd = olv_Fo90.LH12.disl.Vd ; % activation volume 
nd = olv_Fo90.LH12.disl.nd ; % stress exponent
alf_d = olv_Fo90.LH12.disl.alf_d ; 

% Dislocation-accommodated grain boundary sliding (GBS)
Ag = olv_Fo90.LH12.gbs.Ag ; % preexponential for coble diffusion creep
Qg = olv_Fo90.LH12.gbs.Qg ;% activation energy for coble diffusion creep
Vg = olv_Fo90.LH12.gbs.Vg ; % activation volume for coble diff creep
pg = olv_Fo90.LH12.gbs.pg ; % grain size exponent
ng = olv_Fo90.LH12.gbs.ng ; % stress exponent
alf_g = olv_Fo90.LH12.gbs.alf_g ; % melt fraction sensitivity

end

%  ==============================================

    P_on = 0;
    if P_on==1 ;
        % with pressure effect
        sr_c = (Ac).*sig_MPa.*(d.^(-pc)).*exp(-(Qc+P.*Vc)./(R.*T)) ;
        sr_d = (Ad).*(sig_MPa.^nd).*exp(-(Qd+P.*Vd)./(R.*T)) ;
        sr_g = (Ag).*(sig_MPa.^ng).*(d.^(-pg)).*exp(-(Qg+P.*Vg)./(R.*T)) ;
    elseif P_on==0 ;
    % without pressure effect  H&K 95 ? 
        sr_c = Ac.*sig_MPa.*(d.^(-pc)).*exp(-Qc./(R.*T)) ;
        sr_d = (Ad).*(sig_MPa.^nd).*exp(-(Qd)./(R.*T)).*exp(alf_d.*phi) ;
        sr_g = (Ag).*(sig_MPa.^ng).*(d.^(-pg)).*exp(-(Qg)./(R.*T)) ;
    end

% =============================
% EFFECTS of MELT
phi_c = phi_c ;

% COBLE CREEP
alpha = alf_c ; 
speedfac_diff = speedfac ; 
sr_c = sr_c./speedfac ; % to account for lack of truly melt free samples and the drop at the onset of melting.  
[sr_c_prime] = SR_melt_enhancement(phi,alpha,speedfac_diff,phi_c) ;
sr_c = sr_c_prime.*sr_c ; 

% GBS (disl. accommodated)
alpha = alf_g ; 
speedfac_gbs = speedfac/2 ; % assuming GBS is less sensitive to melt than diff, but still a bit. 
sr_g = sr_g./speedfac ; % to account for lack of truly melt free samples in SC... assuming effect is smaller on dislocation-acc.GBS than on coble creep..
[sr_g_prime] = SR_melt_enhancement(phi,alpha,speedfac_gbs,phi_c);
sr_g = sr_g_prime.*sr_g ; 

alpha = alf_d ; 
speedfac_disl = speedfac/speedfac ; % dislocation creep has no rapid reduction
sr_d = sr_d./speedfac ; %./2.0 ; % to account for lack of truly melt free samples in SC... assuming effect is smaller on dislocation-acc.GBS than on coble creep..
[sr_d_prime] = SR_melt_enhancement(phi,alpha,speedfac_disl,phi_c) ;
sr_d = sr_d_prime.*sr_d ; 

%  ===============================================
% if Opt_FloLaw==0|Opt_FloLaw==1 ;
%     sr_tot = sr_c + sr_d ; 
% elseif Opt_FloLaw==2 ; 
%     sr_tot = sr_c + sr_g ;
% end

sr_tot = sr_c + sr_d + sr_g ;


