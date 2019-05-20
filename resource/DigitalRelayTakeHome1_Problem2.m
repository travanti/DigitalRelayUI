#code adaption for Take Home Test#1 will be created herein

#written by Nathaniel P. Travanti on 5-18-19
#use of examples available for matlab and octave standard
#documentation used extensively
#an example function for circle creation was used as provided by MatWorks Support
#Staff Original post located 
#at: https://www.mathworks.com/matlabcentral/answers/98665-how-do-i-plot-a-circle-with-a-given-radius-and-center


#needed packages
#below is commands to inport packages into octave 
pkg load signal
pkg load image

dos('del resultTakeHome#1.txt')
diary('resultTakeHome#1.txt')
clc
clear all
close all 
M = csvread("fault01_PROB2.csv");
#needs to not be the first line of file can be declared anywhere
#matlab requires at end of file

function out = rcfilt(fc,in,t)
  %Function for preforming the input RC filter
  %using good practice every function of the filter will be declared functionally
  %input variable vector of data points before filtering i[]
  %output LPF result signal o[] 
  %time vector from which to extract sample rate t[]
  %cut off frequency for ~0.707 = 1/sqrt(2) point of signal 
  RC = 1/(2* pi * fc);
  pval = 0; #helper var to hold previous value
  delt = t(2,1) - t(1,1);
  o = [];
  for k = 1:length(in)
    pval = (delt * in(k,1) + RC * pval)/(delt + RC); #recursive call with variable
    o = [o, pval];
    out = o;
  endfor
endfunction
#end

#function to handle the sin filter portion
function out = sinFilt(in,newSamp)
  #Reminder only one cycle of signal is needed for filter range
  k = (1:newSamp);
  sinCoeff = (2/newSamp) * sin(2 * pi * k/newSamp);
  out = filter(sinCoeff,1,in);
  
endfunction 
#end

#function to handle the cos filter portion
function out = cosFilt(in,newSamp)
  k = (1:newSamp);
  cosCoeff = (2/newSamp) * cos(2 * pi * k/newSamp);
  out = filter(cosCoeff,1,in); #for why one is used see feedback loop on filter
  
endfunction
#end

#used from post on mathworks, option #1 presented
#https://www.mathworks.com/matlabcentral/answers/98665-how-do-i-plot-a-circle-with-a-given-radius-and-center
function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
endfunction
#end

#helper function to keep ZcCalculation from getting sign errors introduced
function out = zcCalc(Vpol,V,I,Zr)
out = 0.5 * (Zr + ((V) - Vpol)./(I));
  
endfunction
#end

#helper function to keep ZrCalculation from getting sign errors introduced 
function out = zrCalc(Vpol,V,I,Zr)
out = 0.5 * abs(Zr - ((V) - Vpol)./(I));

endfunction
#end 
angAdj = -1;
CrOrSelf = -1;
#disp('Adjust to Va? (1 = yes), (2 = no)')

###############################setting input####################################
#setting input for if angle adjustment is being applied or not
while angAdj > 1 || angAdj < 0
angAdj = input('Reference phasor to Va? (1 = yes), (0 = no)')
endwhile
if angAdj == 1
  disp('User selected yes.')
else
  disp('User selected no.')
end
disp('');

#check to see if cross or self polarized relay element is activated
while CrOrSelf > 1 || CrOrSelf < 0
CrOrSelf = input('Cross polarized or self polarized relay (1 = Cross), (0 = Self)')
endwhile
if CrOrSelf == 1
  disp('User selected cross.')
else
  disp('User selected self.')
end
disp('');

#read in data to matrixes
time = M(:,1); #used for timing and native sample rate
VaOg = M(:,2);
VbOg = M(:,3);
VcOg = M(:,4);
IaOg = M(:,5); #prefilter, primary side.
IbOg = M(:,6);
IcOg = M(:,7);

CTR = 800/5 ; #current transformer ratio primary to secondary
VTR = (500e3 / sqrt(3)) / 67; #voltage transformer ratio primary to secondary
ZTR = VTR/CTR;

IaSec = IaOg / CTR; #convert to secondary current value
IbSec = IbOg / CTR; 
IcSec = IcOg / CTR;
VaSec = VaOg / VTR; #convert to secondary voltage value
VbSec = VbOg / VTR;
VcSec = VcOg / VTR;

#plot unfiltered results current
figure
plot(time,IaSec,time,IbSec,time,IcSec)
xlim([0 0.1])
grid on
title("Currents Ia, Ib, Ic no filter Secondary Value")
xlabel("time (seconds)")
ylabel("Current (amps)")
legend("Current Ia Secondary","Current Ib Secondary","Current Ic Secondary")

#plot unfiltered results voltage
figure
plot(time,VaSec,time,VbSec,time,VcSec)
xlim([0 0.1])
grid on
title("Volts Va, Vb, Vc no filter Secondary Value")
xlabel("time (seconds)")
ylabel("Volts (volts)")
legend("Volts Va Secondary","Volts Vb Secondary","Volts Vc Secondary")

fc = 600; #selected cut off frequency
IaLPF = rcfilt(fc,IaSec,time);
IbLPF = rcfilt(fc,IbSec,time);
IcLPF = rcfilt(fc,IcSec,time);
VaLPF = rcfilt(fc,VaSec,time);
VbLPF = rcfilt(fc,VbSec,time);
VcLPF = rcfilt(fc,VcSec,time);

figure
plot(time,IaLPF,time,IbLPF,time,IcLPF)
xlim([0 0.1])
grid on
title("Current Ia, Ib, Ic after LPF filtering")
xlabel("time (seconds)")
ylabel("Current (amps)")
legend("Current Ia LPF", "Current Ib LPF", "Current Ic LPF")

figure
plot(time,VaLPF,time,VbLPF,time,VcLPF)
xlim([0 0.1])
grid on
title("Voltage Va, Vb, Vc after LPF filtering")
xlabel("time (seconds)")
ylabel("Voltage (Volts)")
legend("Volts Va LPF", "Volts Vb LPF", "Volts Vc LPF")

#implementation of interpolation (resample)
#this makes our signal a digital signal from its virtual "analog signal"
#which is a very high sample rate digital... 16 samples per cycle at 60 cycl a second
#digital

newSamp = 16; #new number of samples per cycle of wave
fs = newSamp * 60; #convert to frequency in terms of time
newDelT = 1/fs;
timeDsmp = (0:newDelT:time(length(time)));
IaFr = interp1(time,IaLPF,timeDsmp); %resampling from original time, of source sig, to new time sample
IbFr = interp1(time,IbLPF,timeDsmp);
IcFr = interp1(time,IcLPF,timeDsmp);
VaFr = interp1(time,VaLPF,timeDsmp);
VbFr = interp1(time,VbLPF,timeDsmp);
VcFr = interp1(time,VcLPF,timeDsmp);

##figure
##plot(timeDsmp,IaFr,timeDsmp,IbFr,timeDsmp,IcFr)
##xlim([0 0.1])
##grid on;
##title("Current Ia, Ib, Ic after FIR filtering")
##xlabel("time (seconds)")
##ylabel("Current (amps)")
##legend("Current Ia postFIR", "Current Ib postFIR", "Current Ic postFIR")
##
##figure
##plot(timeDsmp,VaFr,timeDsmp,VbFr,timeDsmp,VcFr)
##xlim([0 0.1])
##grid on;
##title("Voltage Va, Vb, Vc after FIR filtering")
##xlabel("time (seconds)")
##ylabel("Voltage (Volts)")
##legend("Volts Va postFIR", "Volts Vb postFIR", "Volts Vc postFIR")

k = (1:newSamp); #to define our filter as a 1 cycle filter we force coefficent 
sinCoeff = (2/newSamp) * sin(2 * pi * k/newSamp);
cosCoeff = (2/newSamp) * cos(2 * pi * k/newSamp);
Iax = filter(cosCoeff,1,IaFr); #x component of the phasor
Iay = filter(sinCoeff,1,IaFr); #y component of the phasor

Ibx = filter(cosCoeff,1,IbFr); #etc.
Iby = filter(sinCoeff,1,IbFr); 

Icx = filter(cosCoeff,1,IcFr); #ic
Icy = filter(sinCoeff,1,IcFr);

Vax = filter(cosCoeff,1,VaFr); #Va
Vay = filter(sinCoeff,1,VaFr); 

Vbx = filter(cosCoeff,1,VbFr); #Vb
Vby = filter(sinCoeff,1,VbFr);

Vcx = filter(cosCoeff,1,VcFr); #Vc
Vcy = filter(sinCoeff,1,VcFr); 

#calculate out phasor components from x and ys given by filter
j= sqrt(-1); #imaginary
rmsConst = 1/ sqrt(2);

#declare line data
Zl1 = Zl2 = 0.073 + j * 0.8;
Zl0 = 0.1 + j*2.6;
Er = 0.994 * 500e3 *(cos(deg2rad(0)) + j * sin(deg2rad(0)));
Es = 1.002 * 500e3 * (cos(deg2rad(12.5)) + j * sin(deg2rad(12.5)));
Z1s = Z2s = 1 + j * 10;
Z0s = 3 + j * 20;
Z1r = Z2r = j * 20;
Z0r = j * 10;
k0 = (1/3) * ((Zl0 - Zl1)/(Zl1));
L2 = 100;

Ia = rmsConst * (Iax + j * Iay) .* exp(-j * 2 * pi * 60 * timeDsmp); #* is similar to convolving the signal in this case
Ib = rmsConst * (Ibx + j * Iby) .* exp(-j * 2 * pi * 60 * timeDsmp); 
Ic = rmsConst * (Icx + j * Icy) .* exp(-j * 2 * pi * 60 * timeDsmp);
Va = rmsConst * (Vax + j * Vay) .* exp(-j * 2 * pi * 60 * timeDsmp); #Va voltage
Vb = rmsConst * (Vbx + j * Vby) .* exp(-j * 2 * pi * 60 * timeDsmp); 
Vc = rmsConst * (Vcx + j * Vcy) .* exp(-j * 2 * pi * 60 * timeDsmp);
if angAdj > 0
angAdjFact = angle(Va); #in radians will be left
angAdjFact = (cos(angAdjFact) + j* sin(angAdjFact));

iter = length(Va) 
  while iter > 0
    Va(1,iter) = Va(1,iter) / angAdjFact(1,iter);
    Vb(1,iter) = Vb(1,iter) / angAdjFact(1,iter);
    Vc(1,iter) = Vc(1,iter) / angAdjFact(1,iter);
    Ia(1,iter) = Ia(1,iter) / angAdjFact(1,iter);
    Ib(1,iter) = Ib(1,iter) / angAdjFact(1,iter);
    Ic(1,iter) = Ic(1,iter) / angAdjFact(1,iter);
    iter = iter - 1;
   end
endif

C_radDeg = 180/pi;
###plot the original secondary vs the final result phasor 
###Voltage Va magnitude vs its Original
##figure;plot(time,VaSec),grid;hold on;plot(timeDsmp,abs(Va),'r')
##xlim([0 0.1])
##title('Original instantaneous Phase-A voltage vs Magnitude of Phasor')
##legend('Volts Va secondary','Volts Va phasor Mag')
##ylabel('Voltage in volts')
##xlabel('Time in seconds')
##hold off;
##
###Current Ia magnitude vs its Original
##figure;plot(time,IaSec),grid;hold on;plot(timeDsmp,abs(Ia),'r')
##xlim([0 0.1])
##title('Original instantaneous Phase-A current vs Magnitude of Phasor')
##legend('Current Ia secondary','Current Ia phasor Mag')
##ylabel('Current (amps)')
##xlabel('Time (secs)')
##hold off;

#magnitudes of all results (subplotted)
figure;
subplot(3,2,1)
plot(time,VaSec);hold on;plot(timeDsmp,abs(Va),'r')
xlabel('Time in seconds')
ylabel('Va Secondary, Va Mag')
hold off;

subplot(3,2,2)
plot(time,VbSec);hold on;plot(timeDsmp,abs(Vb),'r')
xlabel('Time in seconds')
ylabel('Vb Secondary, Vb Mag')
hold off;

subplot(3,2,3)
plot(time,VcSec);hold on;plot(timeDsmp,abs(Vc),'r')
xlabel('Time in seconds')
ylabel('Vc Secondary, Vc Mag')
hold off;

subplot(3,2,4)
plot(time,IaSec);hold on;plot(timeDsmp,abs(Ia),'r')
xlabel('Time in seconds')
ylabel('Ia Secondary, Ia Mag')
hold off;

subplot(3,2,5)
plot(time,IbSec);hold on;plot(timeDsmp,abs(Ib),'r')
xlabel('Time in seconds')
ylabel('Ib Secondary, Ib Mag')
hold off;

subplot(3,2,6)
plot(time,IcSec);hold on;plot(timeDsmp,abs(Ic),'r')
xlabel('Time in seconds')
ylabel('Ic Secondary, Ic Mag')
hold off;


#angles of all results (subplotted)
figure;
subplot(3,2,1)
plot(timeDsmp,angle(Va)*C_radDeg,'r') #no primary value to compare to
xlabel('Time in seconds')
ylabel('Va angle in deg')

subplot(3,2,2)
plot(timeDsmp,angle(Vb)*C_radDeg,'r')
xlabel('Time in seconds')
ylabel('Vb angle in deg')

subplot(3,2,3)
plot(timeDsmp,angle(Vc)*C_radDeg,'r')
xlabel('Time in seconds')
ylabel('Vc angle in deg')

subplot(3,2,4)
plot(timeDsmp,angle(Ia)*C_radDeg,'r')
xlabel('Time in seconds')
ylabel('Ia angle in deg')

subplot(3,2,5)
plot(timeDsmp,angle(Ib)*C_radDeg,'r')
xlabel('Time in seconds')
ylabel('Ib angle in deg')

subplot(3,2,6)
plot(timeDsmp,angle(Ic)*C_radDeg,'r')
xlabel('Time in seconds')
ylabel('Ic angle in deg')

#calculations for distance relay elements

#we'll assume the relay is going to grab the next to last value for calculations
#this is for analysis simplicity and in actuality would have its on subroutine
#probably chosing the sample after some predetermined amount of stability is shown
#in X number of samples... (ask Perez abut this...?)

MHOAB = MHOBC = MHOCA = MHOAG = MHOBG = MHOCG = zeros(1,length(Va)); #All initally not tripped.
##ElementIa = ''; #declare vaiables to hold info for individual compass plots below
##ElementIb = '';
##ElementIc = '';
##ElementVa = '';
##ElementVb = '';
##ElementVc = '';

Ires = Ia + Ib + Ic;
Zr = (0.8 * L2 * (Zl1))/ZTR;
disp('')
disp(['Zr setting is defined as: ' num2str(abs(Zr)) ' ohms /angle ' num2str(rad2deg(angle(Zr))) ' deg'])

#report out the final found phase current and the phase voltages
disp('The penultimate value of Ia is: ')
disp([num2str(abs(Ia(1,length(Ia) - 1))) ' A /angle ' num2str(rad2deg(angle(Ia(1,length(Ia) - 1)))) ' deg'])
disp('')
disp('The penultimate Value of Ib is: ')
disp([num2str(abs(Ib(1,length(Ib) - 1))) ' A /angle ' num2str(rad2deg(angle(Ib(1,length(Ib) - 1)))) ' deg'])
disp('')
disp('The penultimate Value of Ic is: ')
disp([num2str(abs(Ic(1,length(Ic) - 1))) ' A /angle ' num2str(rad2deg(angle(Ic(1,length(Ic) - 1)))) ' deg'])
disp('')
disp('The penultimate Value of Va is: ')
disp([num2str(abs(Va(1,length(Va) - 1))) ' V /angle ' num2str(rad2deg(angle(Va(1,length(Va) - 1)))) ' deg'])
disp('')
disp('The penultimate Value of Vb is: ')
disp([num2str(abs(Vb(1,length(Vb) - 1))) ' V /angle ' num2str(rad2deg(angle(Vb(1,length(Vb) - 1)))) ' deg'])
disp('')
disp('The penultimate Value of Vc is: ')
disp([num2str(abs(Vc(1,length(Vc) - 1))) ' V /angle ' num2str(rad2deg(angle(Vc(1,length(Vc) - 1)))) ' deg'])
disp('')


#calculate the relay values based off the setting chosen by the user at the
#beginning of the program

if CrOrSelf > 0 
#Trip determination of the Cross polarized relay
neg90Cart = (1*(cos(deg2rad(-90)) + j * sin(deg2rad(-90))));
pos90Cart = (1*(cos(deg2rad(90)) + j * sin(deg2rad(90))));

#Iag
Zag = Va ./ (Ia + k0 * Ires);
s1 = Zr * (Ia + k0 * Ires) - Va;
s2 = (Vb - Vc) * pos90Cart;
ZcIag = zcCalc(s2,(Va(1,193)),(Ia(1,193) + k0 * Ires(1,193)),Zr);
ZrIag = zrCalc(s2,(Va(1,193)),(Ia(1,193) + k0 * Ires(1,193)),Zr);

count = 1;
while count <= length(Va)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOAG(1,count) = 1;
  endif
count = count + 1;
endwhile
#end loop
#Iag impedance plane plot
ZrRange = linspace(0,real(Zr)) ;
ZrPlot = ZrRange * ((imag(Zr)/real(Zr)));
ZagRange = linspace(0,real(Zag(1,193))); #change per element (CPE)
ZagPlot = ZagRange * ((imag(Zag(1,193))/real(Zag(1,193)))); #CPE
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
pCross = circle(real(ZcIag(1,193)),imag(ZcIag(1,193)),ZrIag(1,193));
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis #CPE(?)
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis #CPE(?)
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZagRange, ZagPlot) #CPE
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Expanded Mho by Cross polarization","Imaginary Axis","Real Axis","Zag apparent","Zr Setting") #CPE
title('Iag Element Cross Polarized Mho Results') #CPE
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal") #CPE(?)

#Ibg
Zbg = Vb ./ (Ib + k0 * Ires);
s1 = Zr * (Ib + k0 * Ires) - Vb;
s2 = (Vc - Va) * pos90Cart;
ZcIbg = zcCalc(s2,(Vb(1,193)),(Ib(1,193) + k0 * Ires(1,193)),Zr);
ZrIbg = zrCalc(s2,(Vb(1,193)),(Ib(1,193) + k0 * Ires(1,193)),Zr);

count = 1
while count <= length(Vb)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOBG(1,count) = 1;
  endif
  count = count + 1;
endwhile
#end loop

#Ibg impedance plane plot
ZbgRange = linspace(0,real(Zbg(1,193)));
ZbgPlot = ZbgRange * ((imag(Zbg(1,193))/real(Zbg(1,193))));
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
pCross = circle(real(ZcIbg(1,193)),imag(ZcIbg(1,193)),ZrIbg(1,193));
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZbgRange, ZbgPlot)
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Expanded Mho by Cross polarized","Imaginary Axis","Real Axis","Zbg apparent","Zr Setting")
title('Ibg Element Cross Polarized Mho Results')
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal")

#Icg
Zcg = Vc ./ (Ic + k0 * Ires);
s1 = Zr * (Ic + k0 * Ires) - Vc;
s2 = (Va - Vb) * pos90Cart;
ZcIcg = zcCalc(s2,(Vc),(Ic + k0 * Ires),Zr);
ZrIcg = zrCalc(s2,(Vc),(Ic + k0 * Ires),Zr);

count = 1;
while count <= length(Vc)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOCG(1,count) = 1;
  endif
  count = count + 1;
endwhile
#Icg impedance plane plot
ZcgRange = linspace(0,real(Zcg(1,193))); #change per element (CPE)
ZcgPlot = ZcgRange * ((imag(Zcg(1,193))/real(Zcg(1,193)))); #CPE
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
pCross = circle(real(ZcIcg(1,193)),imag(ZcIcg(1,193)),ZrIcg(1,193));
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis #CPE(?)
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis #CPE(?)
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZcgRange, ZcgPlot) #CPE
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Expanded Mho by Cross Polarized","Imaginary Axis","Real Axis","Zcg apparent","Zr Setting") #CPE
title('Icg Element Cross Polarized Mho Results') #CPE
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal") #CPE(?)


#Vab
Zab = (Va - Vb) ./ (Ia - Ib);
s1 = Zr * (Ia - Ib)  - (Va - Vb);
s2 = (Vc) * neg90Cart;
ZcVab = zcCalc(s2,(Va(1,193) - Vb(1,193)),(Ia(1,193) - Ib(1,193)),Zr);
ZrVab = zrCalc(s2,(Va(1,193) - Vb(1,193)),(Ia(1,193) - Ib(1,193)),Zr);
count = 1;
while count <= length(Va)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOAB(1,count) = 1;
  endif
  count = count + 1;
endwhile
#Vab impedance plane plot
#Zr plot reference impedance
ZrRange = linspace(0,real(Zr)) ;
ZrPlot = ZrRange * ((imag(Zr)/real(Zr)));
#Zapparent plot apparent plot for the element in question
ZabRange = linspace(0,real(Zab(1,193))); #change per element (CPE)
ZabPlot = ZabRange * ((imag(Zab(1,193)))/real(Zab(1,193))); #CPE
#begin plotting
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
pCross = circle(real(ZcVab(1,193)),imag(ZcVab(1,193)),ZrVab(1,193));
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis #CPE(?)
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis #CPE(?)
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZabRange, ZabPlot) #CPE
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Expanded Mho by Cross Polarized","Imaginary Axis ","Real Axis","Zab apparent","Zr Setting") #CPE
title('Vab Element Cross Polarized Mho Results') #CPE
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal") #CPE(?)

#Vbc
Zbc = (Vb - Vc) ./ (Ib - Ic);
s1 = Zr * (Ib - Ic) - (Vb - Vc);
s2 = (Va) * neg90Cart;
ZcVbc = zcCalc(s2,(Vb(1,193) - Vc(1,193)),(Ib(1,193) - Ic(1,193)),Zr);
ZrVbc = zrCalc(s2,(Vb(1,193) - Vc(1,193)),(Ib(1,193) - Ic(1,193)),Zr);
count = 1;
while count <= length(Vb)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOBC(1,count) = 1;
  endif
  count = count + 1;
endwhile

#Vbc impedance plane plot
ZrRange = linspace(0,real(Zr)) ;
ZrPlot = ZrRange * ((imag(Zr)/real(Zr)));
ZbcRange = linspace(0,real(Zbc(1,193))); #change per element (CPE)
ZbcPlot = ZbcRange * ((imag(Zbc(1,193)))/real(Zbc(1,193))); #CPE
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
PCross = circle(real(ZcVab(1,193)),imag(ZcVab(1,193)),ZrVbc(1,193));
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis #CPE(?)
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis #CPE(?)
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZbcRange, ZbcPlot) #CPE
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Expanded Mho by Cross Polarized","Imaginary Axis ","Real Axis","Zbc apparent","Zr Setting") #CPE
title('Vbc Element Cross Polarized Mho Results') #CPE
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal") #CPE(?)

#Vca
Zca = (Vc - Va) ./ (Ic - Ia);
s1 = Zr * (Ic - Ia) - (Vc - Va);
s2 = (Vb) * neg90Cart;
ZcVca = zcCalc(s2,(Vc(1,193) - Va(1,193)),(Ic(1,193) - Ia(1,193)),Zr);
ZrVca = zrCalc(s2,(Vc(1,192) - Va(1,193)),(Ic(1,193) - Ia(1,193)),Zr);
count = 1;
while count <= length(Vc)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOCA(1,count) = 1;
  endif
  count = count + 1;
endwhile
#plot impedance planes for Vca
ZrRange = linspace(0,real(Zr)) ;
ZrPlot = ZrRange * ((imag(Zr)/real(Zr)));
ZcaRange = linspace(0,real(Zca(1,193)));
ZcaPlot = ZcaRange * ((imag(Zca(1,193))/real(Zca(1,193))));
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
pCross = circle(real(ZcVca(1,193)),imag(ZcVca(1,193)),ZrVca(1,193));
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZcaRange, ZcaPlot)
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Expanded Mho by Corss Polarized","Imaginary Axis ","Real Axis","Zca apparent","Zr Setting")
title('Vca Element Cross Polarized Mho Results')
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal")

else
#Trip determination of the Self polarized relay

#Iag
Zag = Va ./ (Ia + k0 * Ires);
s1 = Zr * (Ia + k0 * Ires) - Va;
s2 = Va;
count = 1;
while count <= length(Va)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOAG(1,count) = 1;
  endif
count = count + 1;
endwhile
#end loop
#Iag impedance plane plot
ZrRange = linspace(0,real(Zr)) ;
ZrPlot = ZrRange * ((imag(Zr)/real(Zr)));
ZagRange = linspace(0,real(Zag(1,193))); #change per element (CPE)
ZagPlot = ZagRange * ((imag(Zag(1,193))/real(Zag(1,193)))); #CPE
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis #CPE(?)
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis #CPE(?)
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZagRange, ZagPlot) #CPE
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Imaginary Axis ","Real Axis","Zag apparent","Zr Setting") #CPE
title('Iag Element Self Polarized Mho Results') #CPE
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal") #CPE(?)

#Ibg
Zbg = Vb ./ (Ib + k0 * Ires);
s1 = Zr * (Ib + k0 * Ires) - Vb;
s2 = Vb;
count = 1
while count <= length(Vb)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOBG(1,count) = 1;
  endif
  count = count + 1;
endwhile
#end loop

#Ibg impedance plane plot
ZrRange = linspace(0,real(Zr)) ;
ZrPlot = ZrRange * ((imag(Zr)/real(Zr)));
ZbgRange = linspace(0,real(Zbg(1,193)));
ZbgPlot = ZbgRange * ((imag(Zbg(1,193))/real(Zbg(1,193))));
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZbgRange, ZbgPlot)
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Imaginary Axis ","Real Axis","Zbg apparent","Zr Setting")
title('Ibg Element Self Polarized Mho Results')
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal")

#Icg
Zcg = Vc ./ (Ic + k0 * Ires);
s1 = Zr * (Ic + k0 * Ires) - Vc;
s2 = Vc;
count = 1;
while count <= length(Vc)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOCG(1,count) = 1;
  endif
  count = count + 1;
endwhile
#Icg impedance plane plot
ZrRange = linspace(0,real(Zr)) ;
ZrPlot = ZrRange * ((imag(Zr)/real(Zr)));
ZcgRange = linspace(0,real(Zcg(1,193))); #change per element (CPE)
ZcgPlot = ZcgRange * ((imag(Zcg(1,193))/real(Zcg(1,193)))); #CPE
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis #CPE(?)
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis #CPE(?)
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZcgRange, ZcgPlot) #CPE
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Imaginary Axis ","Real Axis","Zcg apparent","Zr Setting") #CPE
title('Icg Element Self Polarized Mho Results') #CPE
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal") #CPE(?)


#Vab
Zab = (Va - Vb) ./ (Ia - Ib);
s1 = Zr * (Ia - Ib)  - (Va - Vb);
s2 = (Va - Vb);
count = 1;
while count <= length(Va)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOAB(1,count) = 1;
  endif
  count = count + 1;
endwhile
#Vab impedance plane plot
ZrRange = linspace(0,real(Zr)) ;
ZrPlot = ZrRange * ((imag(Zr)/real(Zr)));
ZabRange = linspace(0,real(Zab(1,193))); #change per element (CPE)
ZabPlot = ZabRange * ((imag(Zab(1,193))/real(Zab(1,193)))); #CPE
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis #CPE(?)
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis #CPE(?)
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZabRange, ZabPlot) #CPE
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Imaginary Axis ","Real Axis","Zab apparent","Zr Setting") #CPE
title('Vab Element Self Polarized Mho Results') #CPE
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal") #CPE(?)

#Vbc
Zbc = (Vb - Vc) ./ (Ib - Ic);
s1 = Zr * (Ib - Ic) - (Vb - Vc);
s2 = (Vb - Vc);
count = 1;
while count <= length(Vb)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOBC(1,count) = 1;
  endif
  count = count + 1;
endwhile

#Vbc impedance plane plot
ZrRange = linspace(0,real(Zr)) ;
ZrPlot = ZrRange * ((imag(Zr)/real(Zr)));
ZbcRange = linspace(0,real(Zbc(1,193))); #change per element (CPE)
ZbcPlot = ZbcRange * ((imag(Zbc(1,193))/real(Zbc(1,193)))); #CPE
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis #CPE(?)
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis #CPE(?)
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZbcRange, ZbcPlot) #CPE
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Imaginary Axis ","Real Axis","Zbc apparent","Zr Setting") #CPE
title('Vbc Element Self Polarized Mho Results') #CPE
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal") #CPE(?)

#Vca
Zca = (Vc - Va) ./ (Ic - Ia);
s1 = Zr * (Ic - Ia) - (Vc - Va);
s2 = (Vc - Va);
count = 1;
while count <= length(Vc)
  if(real(s1(1,count) * conj(s2(1,count))) > 0)
    MHOCA(1,count) = 1;
  endif
  count = count + 1;
endwhile
#plot impedance planes for Vca
ZrRange = linspace(0,real(Zr)) ;
ZrPlot = ZrRange * ((imag(Zr)/real(Zr)));
ZcaRange = linspace(0,real(Zca(1,193)));
ZcaPlot = ZcaRange * ((imag(Zca(1,193))/real(Zca(1,193))));
figure
p = circle(real(Zr)/2,imag(Zr)/2,abs(Zr)/2);
xL = xlim;
yL = ylim;
YLine = line([0 0], yL);  %y-axis
set(YLine,'Color','black') #set the origin lines to black
XLine = line(xL, [0 0]);  %x-axis
set(XLine,'Color','black') #set the origin lines to black

hold on
plot(ZcaRange, ZcaPlot)
plot(ZbgRange, ZbgPlot)
plot(ZrRange,ZrPlot)
hold off
legend("Mho Circle","Imaginary Axis ","Real Axis","Zca apparent","Zr Setting")
title('Vca Element Self Polarized Mho Results')
xlabel('R Real ohms')
ylabel('X Imaginary ohms')
axis("equal")

endif

#plot relay results
#I a to ground fault
figure
subplot(3,2,1)
plot(timeDsmp, MHOAG)
title('MHOAG at times')
xlabel('time seconds')
ylabel('relay MHOAG status')
axis([0.02, 0.1,0,1]) #reminder form xstrt, xend, ystrt, yend
grid on

#I b to ground fault
subplot(3,2,2)
plot(timeDsmp, MHOBG)
title('MHOBG at times')
xlabel('time seconds')
ylabel('relay MHOBG status')
axis([0.02,0.1,0,1])
grid on

#I c to ground fault
subplot(3,2,3)
plot(timeDsmp, MHOCG)
title('MHOCG at times')
xlabel('time seconds')
ylabel('relay MHOCG status')
axis([0.02,0.1,0,1])
grid on


#V phase a to phase b fault
subplot(3,2,4)
plot(timeDsmp, MHOAB)
title('MHOAB at times')
xlabel('time seconds')
ylabel('relay MHOAB status')
axis([0.02,0.1,0,1])
grid on

#V phase b to phase c fault
subplot(3,2,5)
plot(timeDsmp, MHOBC)
title('MHOBC at times')
xlabel('time seconds')
ylabel('relay MHOBC status')
axis([0.02,0.1,0,1])
grid on

#V phase c to phase a fault
subplot(3,2,6)
plot(timeDsmp, MHOCA)
title('MHOCA at times')
xlabel('time seconds')
ylabel('relay MHOCA status')
axis([0.02,0.1,0,1])
grid on

diary off
#fin