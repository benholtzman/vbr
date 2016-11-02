function [SR_phi_enh] = sr_melt_enhancement(phi,alpha,x_phi_c,phi_c)
% calculates a reduction factor for viscosity (not a function...) 
% error function for step * exp(-alpha*phi)
% with a variable parameter for the amplitude of the step. 

%boostmag = x_phi_c ;
%dropmag = (1-1/x_phi_c) ;
a = log(x_phi_c) ; 
ratefac = 1./(phi_c) ; % this works out pretty well, actually, for phi_c = 0.005 
%step = (boostmag.*erf(phi.*boostratefac)) ;
step = a.*erf(phi.*ratefac) ;
slope = alpha.*phi;

ln_SR_phi_enh = slope + step ;

SR_phi_enh = exp( ln_SR_phi_enh ) ; 

end
