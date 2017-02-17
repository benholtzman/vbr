%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phis = stag_phi(Vark,settings)
% stagger phi s.t. phi = phimin is maintained
   phis = stag(Vark.phi); Tsol = stag(Vark.Tsol); T = stag(Vark.T); 
%    phis(T<Tsol)=settings.phimin;
end