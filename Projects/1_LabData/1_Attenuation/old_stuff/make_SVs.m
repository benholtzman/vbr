%% MAKE VECTORS =======================================================
% THESE WILL GET PROGRESSIVELY SELF CONSISTENT
% This script must be totally general... replace things that will change
% after implementing it. 
n = 1 ;
n_freq = length(f) ; 
%% ====================================================
% INTENSIVE state variables ==========================
%Temperature
VBR.ISV.T_K = T_K.*ones(n,1) ; % in GPa

% density
VBR.ISV.rho = rho.*ones(n,1) ; % in GPa

% Pressure
VBR.ISV.P_GPa = P_GPa.*ones(n,1) ; % in GPa

% stress 
%sig_0 = 1e6 ; % MPa (for flow law calc)
VBR.ISV.sig_Pa = sig_0.*ones(n,1) ; 
VBR.ISV.sig_MPa = (VBR.ISV.sig_Pa)./1e6 ; 
% frequency (scalar or vector-- becomes the first index in matrix that is fed through VBR.. 
% which is constructed in the subroutines for anelasticity, called by VBR_spine.m
VBR.ISV.f = f ; %.*ones(n_freq,1) ; 

template = ones(size(VBR.ISV.P_GPa)) ;
VBR.Gu_0 = Gu_0.*template.*1e9 ; 

%% ====================================================
% COMPOSITIONAL state variables ======================

% water content 
%Ch2o_0 = 0 ; % 
VBR.CSV.Ch2o = Ch2o_0.*ones(n,1) ; 

%% ====================================================
% STRUCTURAL state variables =========================

% MELT fraction
%phi_0 = 0.00 ; % 
VBR.SSV.phi = phi_0.*ones(n,1) ; 

% grain size
%dg_0 = 5e3 ; % 5 mm in microns
%dg_0 = 3.1 ; % 3.1 microns, for comparison to Jackson+Faul2010
VBR.SSV.dg_um = dg_um.*ones(n,1) ; 







