#assumes main relay loop has run and read in variables to be accessed by this script
#adds m plots to the report
#line data is hard coded and will need to be set y user each time.

j = sqrt(-1);
i = sqrt(-1);
kk=180/pi;
a=1*exp(j*120/kk);
CTR = 800/5 ; #current transformer ratio primary to secondary
VTR = (500e3 / sqrt(3)) / 67; #voltage transformer ratio primary to secondary
ZTR = VTR/CTR;


Zs1 = Zs2 = 1 + j *10;
Zs0 = 2 + j * 30;
Zr1 = Zr2 = j *20;
Zr0 = j * 10;
Es = 1.002 * 500e3 *( cos(deg2rad(12.5)) + j * sin(deg2rad(12.5)));
Er = 0.994 * 500e3 * (cos((deg2rad(0))) + j * sin(deg2rad(0)));
line1 = 50; # transmissin line side one length
line2 = 100; #transmission line side 2 length
d2F2 = 10;
d1F1 = 20;

#line data
zl1 = zl2 = 0.073+j*0.8; #per mile
zl0 = 0.1 + j * 2.6; #per mile
cl  = 0.013e-6; #per mile
ZL11 = ZL21 = zl1 * line1;
ZL12 = ZL22 = zl1 * line2;
ZL01 = zl0 * line1; 
ZL02 = zl0 * line2;

k01 = (1/3) *(((ZL01) - (ZL11)) / (ZL11));
k02 = (1/3) *(((ZL02) - (ZL12)) / (ZL12));
Vab = Va - Vb;
Vbc = Vb - Vc;
Vca = Vc - Va;
Iab = Ia - Ib;
Ibc = Ib - Ic;
Ica = Ic - Ia;
Ires = Ia + Ib + Ic;
#assume Kc is 1 and that the compensation angle tau is 0.


#AG element 
mAG1 =abs(imag(Va./(Ia + k01 .* Ires))./imag(ZL11));
figure
plot(timeDsmp, mAG1)
xlabel("time in sec")
ylabel("value of m")
title("AG element m estimate plot Line1")
disp(['estimated m is ' num2str(mAG1(1,length(mAG1)-1))])
disp(['estimated distance (miles) is ' num2str(mAG1(1,length(mAG1)-1) * line1)])


mAG2 =abs(imag(Va./(Ia + k02 .* Ires))./imag(ZL12));
figure
plot(timeDsmp, mAG2)
xlabel("time in sec")
ylabel("value of m")
title("AG element m estimate plot Line2")
disp(['estimated m is ' num2str(mAG2(1,length(mAG2)-1))])
disp(['estimated distance (miles) is ' num2str(mAG2(1,length(mAG2)-1) * line2)])


#BG element
mBG1 =abs(imag(Vb./(Ib + k0 .* Ires))./imag(ZL11));
figure
plot(timeDsmp, mBG1)
xlabel("time in sec")
ylabel("value of m")
title("BG element m estimate plot Line1")
disp(['estimated m is ' num2str(mBG1(1,length(mBG1)-1))])
disp(['estimated distance (miles) is ' num2str(mBG1(1,length(mBG1)-1) * line1)])

mBG2 =abs(imag(Vb./(Ib + k0 .* Ires))./imag(ZL12));
figure
plot(timeDsmp, mBG2)
xlabel("time in sec")
ylabel("value of m")
title("BG element m estimate plot Line2")
disp(['estimated m is ' num2str(mBG2(1,length(mBG2)-1))])
disp(['estimated distance (miles) is ' num2str(mBG2(1,length(mBG2)-1) * line2)])


#CG element
mCG1 =abs(imag(Vc./(Ic + k0 .* Ires))./imag(ZL11));
figure
plot(timeDsmp, mCG1)
xlabel("time in sec")
ylabel("value of m")
title("CG element m estimate plot Line1")
disp(['estimated m is ' num2str(mCG1(1,length(mCG1)-1))])
disp(['estimated distance (miles) is ' num2str(mCG1(1,length(mCG1)-1) * line1)])

mCG2 =abs(imag(Vc./(Ic + k0 .* Ires))./imag(ZL12));
figure
plot(timeDsmp, mCG2)
xlabel("time in sec")
ylabel("value of m")
title("CG element m estimate plot Line2")
disp(['estimated m is ' num2str(mCG2(1,length(mCG2)-1))])
disp(['estimated distance (miles) is ' num2str(mCG2(1,length(mCG2)-1) * line2)])


#AB element


#BC element

