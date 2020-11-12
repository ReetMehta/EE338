%Butterworth Analog LPF parameters
Wc = 1.07;              %cut-off frequency
N = 8;                  %order 

%poles of Butterworth polynomial of degree 8 in the open CLHP 
p1 = Wc*cos(pi/2 + pi/16) + 1i*Wc*sin(pi/2 + pi/16);
p2 = Wc*cos(pi/2 + pi/16) - 1i*Wc*sin(pi/2 + pi/16);
p3 = Wc*cos(pi/2 + pi/16+pi/8) + 1i*Wc*sin(pi/2 + pi/16+pi/8);
p4 = Wc*cos(pi/2 + pi/16+pi/8) - 1i*Wc*sin(pi/2 + pi/16+pi/8);
p5 = Wc*cos(pi/2 + pi/16+2*pi/8) + 1i*Wc*sin(pi/2 + pi/16+2*pi/8);
p6 = Wc*cos(pi/2 + pi/16+2*pi/8) - 1i*Wc*sin(pi/2 + pi/16+2*pi/8);
p7 = Wc*cos(pi/2 + pi/16+3*pi/8) + 1i*Wc*sin(pi/2 + pi/16+3*pi/8);
p8 = Wc*cos(pi/2 + pi/16+3*pi/8) - 1i*Wc*sin(pi/2 + pi/16+3*pi/8);

%Band Edge speifications
fs1 = 34.9;
fp1 = 38.9;
fp2 = 58.9;
fs2 = 62.9;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 330;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi) ;
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

                                                        
%Evaluating Frequency Response of Final Filter
syms s z o;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bpf(s) = analog_lpf((s*s + W0*W0)/(B*s));        %bandstop transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

analog_lpf_f(o) = analog_lpf(1i*o);                 %finding H_analog,lpf as a funtion of angular frequency 
func(o)= abs(analog_lpf_f);
figure, fplot(func)                                 %plotting _analog,lpf against angular frequency
xline(1,':');
xline(-1,':');
xline(1.361,':');
xline(-1.361,':');
yline(0.85,':');
yline(0.15,':');
title('Analog LPF')
ylabel('|H(w)|');
xlabel('Angular Frequency(w)');

analog_bpf_f(o) = analog_bpf(1i*o);                  %finding H_analog,bpf as a funtion of angular frequency 
func(o)= abs(analog_bpf_f);
figure, fplot(func)                                 %plotting _analog,bpf against angular frequency
xline(0.345,':');
xline(0.3883,':');
xline(0.6276,':');
xline(0.6827,':');
yline(0.85,':');
yline(0.15,':');
xlim([-1,1]);
title('Analog BPF')
ylabel('|H(w)|');
xlabel('Angular Frequency(w)');


%coeffs of discrete bsf
[nz, dz] = numden(discrete_bpf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 330e3);
figure, plot(f,abs(H));
xline(34900,':');
xline(38900,':');
xline(58900,':');
xline(62900,':');
yline(0.85,':');
yline(0.15,':');
ylabel('|H(f)|');
xlabel('Frequency(f)');
grid