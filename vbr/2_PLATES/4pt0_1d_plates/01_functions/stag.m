%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fs = stag(F)
 Fs = (F(1:end-1)+F(2:end))/2; % stagger by average
end