%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 4;

% Open CLHP Poles of the Chebyshev Polynomial of order 4
B = asinh(1/epsilon)/N ;
p1 = -sin(pi/(2*N))*sinh(B)+1i*cos(pi/(2*N))*cosh(B);
p2 = -sin(pi/(2*N))*sinh(B)-1i*cos(pi/(2*N))*cosh(B);
p3 = -sin(3*pi/(2*N))*sinh(B)+1i*cos(3*pi/(2*N))*cosh(B);
p4 = -sin(3*pi/(2*N))*sinh(B)-1i*cos(3*pi/(2*N))*cosh(B);        

%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
den = conv(n1,n2); %multiply n1 and n2, which are the two quadratic factors in the denominator
num = [den(5)*sqrt(1/(1+D1))];        % even order, DC Gain set as 1/(1+D1)^0.5

%Band Edge specification
fp1 = 31.1;
fs1 = 35.1;
fs2 = 55.1;
fp2 = 59.1;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 260;
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lsf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bsf(s) = analog_lsf((B*s)/(s*s +W0*W0));     %bandpass transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BPF
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

analog_lsf_f(o) = analog_lsf(1i*o);                 %finding H_analog,lpf as a function of angular frequency
func(o)= abs(analog_lsf_f);
figure, fplot(func)                                 %plotting H_analog,lpf against angular frequency
xline(1,':');
xline(-1,':');
xline(1.349,':');
xline(-1.349,':');
yline(0.85,':');
yline(0.15,':');
title('Analog LPF')
ylabel('|H(w)|');
xlabel('Angular Frequency(w)');hold;


analog_bsf_f(o) = analog_bsf(1i*o);                 %finding H_analog,bsf as a funtion of angular frequency
func(o)= abs(analog_bsf_f);
figure, fplot(func)                                 %plotting H_analog,lpf against angular frequency
xline(0.3945,':');
xline(0.4515,':');
xline(0.5847,':');
xline(0.7853,':');
xline(0.8667,':')
yline(0.85,':');
yline(0.15,':');
xlim([-3,3]);
title('Analog BSF')
ylabel('|H(w)|');
xlabel('Angular Frequency(w)');

%coeffs of discrete BPF
[nz, dz] = numden(discrete_bsf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                       %frequency response in dB

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024,260e3);
figure, plot(f,abs(H),'r');
yline(0.85,':');
yline(0.15,':');
xline(31100,':');
xline(35100,':');
xline(55100,':');
xline(59100,':');
title('Discrete BSF')
ylabel('|H(f)|');
xlabel('Frequency(f)');
grid