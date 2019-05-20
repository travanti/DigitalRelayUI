#takes a complex number and reports it as a polar formed number
#also proivdes a vector with the magnitude and the angle and displays them
#in the command window
function [abs,angle] = reprtComplex(num)
  abs = abs(num)
  angle = rad2deg(angle(num))
  disp(['Magnitude: ' num2str(abs) ' /angle ' num2str(angle)]); disp('');
endfunction
