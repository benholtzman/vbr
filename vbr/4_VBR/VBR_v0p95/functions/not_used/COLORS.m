% COLORS 
global nlines
nlines

blackscale = [ 0.75 0.75 0.75 ; 0.5 0.5 0.5 ; 0 0 0 ] ; 
redscale = [ 1 0 0 ; 0.75 0 0 ; 0.5 0 0 ] ; 
greenscale = [0.2 1 0.2 ; 0.2 0.75 0.2 ; 0.2 0.5 0.2 ] ;

% COLORS
rvec = linspace(0.75,0,nlines) ;
gvec = linspace(0.75,0,nlines) ; 
bvec = linspace(0.75,0,nlines) ; 
blackscale = zeros(nlines,3) ; 
blackscale(:,1) = rvec ; 
blackscale(:,2) = gvec ; 
blackscale(:,3) = bvec ; 
% ==============
rvec = linspace(0.3,1,nlines) ;
gvec = linspace(0,0,nlines) ; 
bvec = linspace(0,0,nlines) ; 
redscale = zeros(nlines,3) ; 
redscale(:,1) = rvec ; 
redscale(:,2) = gvec ; 
redscale(:,3) = bvec ; 
% ==============
rvec = linspace(0,0,nlines) ;
gvec = linspace(0,0,nlines) ; 
bvec = linspace(1,0.5,nlines) ; 
bluescale = zeros(nlines,3) ; 
bluescale(:,1) = rvec ; 
bluescale(:,2) = gvec ; 
bluescale(:,3) = bvec ; 
% ==============
rvec = linspace(0,0,nlines) ;
gvec = linspace(1,0.5,nlines) ; 
bvec = linspace(0,0,nlines) ; 
greenscale = zeros(nlines,3) ; 
greenscale(:,1) = rvec ; 
greenscale(:,2) = gvec ; 
greenscale(:,3) = bvec ; 
