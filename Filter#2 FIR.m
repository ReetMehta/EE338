f_samp = 260e3;

%Band Edge speifications
fp1 = 31.1e3;
fs1 = 35.1e3;
fs2 = 55.1e3;
fp2 = 59.1e3;


Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;
Ws1 = fs1*2*pi/f_samp;
Ws2  = fs2*2*pi/f_samp;
%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-8) / (2.285*(Ws1-Wc1)));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min+14;
disp(n);

%Ideal bandstop impulse response of length "n"

bs_ideal =  ideal_lp(pi,n) -ideal_lp((Ws2+Wc2)/2,n) + ideal_lp((Ws1+Wc1)/2,n);
x=[-26:26];
%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response
%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
figure, plot(f,abs(H));
yline(0.85,':');
yline(0.15,':');
xline(fp1,':');
xline(fs1,':');
xline(fs2,':');
xline(fp2,':');
yline(1.15,':');

title('Discrete BSF - FIR filter')
ylabel('|H(f)|');
xlabel('Frequency(f)');

figure, stem(x,bs_ideal,'.')
xlabel('n');
ylabel('y[n]');
xlim([-27,27])
title('FIR filter time domain sequence');
grid