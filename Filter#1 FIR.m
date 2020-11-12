f_samp = 330e3;

%Band Edge speifications
fs1 = 34.9e3;
fp1 = 38.9e3;
fp2 = 58.9e3;
fs2 = 62.9e3;

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

N_min = ceil((A-8) / (2.285*(Wc1-Ws1)));           %empirical formula for N_min

%Window length for Kaiser Window
n=N_min+20;
disp(n);
%Ideal bandpass impulse response of length "n"
bp_ideal = - ideal_lp((Wc1+Ws1)/2,n) + ideal_lp((Wc2+Ws2)/2,n);

x=[-34:34];

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandPass,1,1024, f_samp);
figure, plot(f,abs(H));
yline(0.85,':');
yline(1.15,':');
yline(0.15,':');
xline(fs1,':');
xline(fp1,':');
xline(fp2,':');
xline(fs2,':');
title('Discrete BPF - FIR filter')
ylabel('|H(f)|');
xlabel('Frequency(f)');

figure, stem(x,bp_ideal,'.')
xlabel('n');
ylabel('y[n]');
xlim([-35,35])
title('FIR filter time domain sequence');
grid