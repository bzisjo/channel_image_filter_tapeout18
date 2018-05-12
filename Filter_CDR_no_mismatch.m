close all
clear all
clc

%% Some parameters

Fs_dig = 40e6;       % This is the clock frequency of the digital module
dT_dig = 1/Fs_dig;
downrate = 200;      % This is used for the MATLAB filters
Fs_matlab = downrate * Fs_dig;
dT_matlab = 1/Fs_matlab;
RF_ampl = 2;         
LOI_ampl = 1;        % LO ampl for In-phase arm
LOQ_ampl = 1;        % LO ampl for Quadrature arm
RF_freq = 2.406e9;   % RF frequency. Refer to BLE specs
IF = 2.5e6;          % Intermediate frequency

SNR_dB = 4;         % Pretty obvious what this is

LO_freq = RF_freq - IF; % Low-side injection: LO < RF

%% Complex BPF with Custom Specs
samplingFrequency = 40e6;       % 40MHz
centerFrequency = 1e6  ;        % 1MHz
transitionWidth = 750e3;        % 750kHz
passbandRipple = 1;             % 1dB Ripple
stopbandAttenuation = 17;       % 17dB stopband attenuation

designSpec = fdesign.lowpass('Fp,Fst,Ap,Ast',...
    centerFrequency-transitionWidth/2, ...
    centerFrequency+transitionWidth/2, ...
    passbandRipple,stopbandAttenuation, ...
    samplingFrequency);
LPF = design(designSpec,'equiripple','SystemObject',true);

Hbp = LPF.Numerator;
N = length(Hbp)-1;
Fc = 0.125;
j = complex(0,1);
Hbp = Hbp.*exp(j*Fc*pi*(0:N));
Hbpreal = real(Hbp);
Hbpimag = imag(Hbp);

%% Quantize the filter to 8-bits
scalingfactorreal = min(abs(-128/min(real(Hbp))), abs(127/max(real(Hbp))));
scalingfactorimag = min(abs(-128/min(imag(Hbp))), abs(127/max(imag(Hbp))));
maxscalingfactor = min(scalingfactorreal, scalingfactorimag);
Hbpreal = round(real(Hbp) * maxscalingfactor);
Hbpimag = round(imag(Hbp) * maxscalingfactor);

%% Export filter coeffs
% These are the image rejection/channel selection filter coefficients
csvwrite("rcoeffs.csv", Hbpreal)
csvwrite("icoeffs.csv", Hbpimag)

%% Generate Data and modulate as RF MSK
num_bits = 220;
T=1e-6;             % This is the symbol rate
freq_shift = 250e3; % This is the frequency deviation or tone spacing for fsk
% Generate some random data
data = randi([0 1],num_bits-6,1);
data = [1 0 1 0 1 0 data']';
data_exp = data;

%Convert data to +/- 1
data = data*2-1;

% Generate first cycle of RF
t = 0:1/Fs_matlab:T-1/Fs_matlab;
phase = 2*pi*(RF_freq+data(1)*freq_shift)*t;
phase_last = phase(end);

% Generate remaining modulated RF cycles
t = 1/Fs_matlab:1/Fs_matlab:T;
for x=2:1:length(data)
    % The phase of each cycle starts where the last one ended
    phase = [phase 2*pi*(RF_freq+data(x)*freq_shift)*t + phase_last];
    phase_last = phase(end);
end

% Scale amplitude to get desired power
FSK = RF_ampl * sin(phase);
t = 0:dT_matlab:(length(FSK)-1)*dT_matlab;

FSK_pure = FSK;

%% Generate side channel tone
sidechannel_freq = 2e6;
side = RF_ampl * sin(2*pi*(RF_freq + sidechannel_freq)*t);
%% Add side channel to main signal
FSK = FSK + side;

%% Generate image channel tone
RF_freq_image = LO_freq-IF;
image = RF_ampl * sin(2*pi*(RF_freq_image)*t);
%% Add image channel to main signal
FSK = FSK + image;

%% Mix and match zone
% FSK_pure, side, image can be added at will to form waveform-of-interest

%% Add RF Noise
FSK_noiseless = FSK;
FSK = awgn(FSK, SNR_dB, 'measured');

% To verify snr
FSKsnr = snr(FSK_noiseless, FSK-FSK_noiseless)
plot(t(1:2000), [FSK(1:2000)' FSK_noiseless(1:2000)'])
legend("Noisy Signal", "Original Signal")
xlabel("Time (s)")
ylabel("Amplitude")
%% Create LO

t = 0:dT_matlab:(length(FSK)-1)*dT_matlab;
LOI = LOI_ampl * sin(2*pi*LO_freq*t) + LOI_ampl;            % In phase: carrier frequency
LOQ = LOQ_ampl * sin(2*pi*LO_freq*t -90*pi/180) + LOQ_ampl; % Quadrature: carrier frequency - 90degrees

%% Downconvert with LO
IFI = LOI .* FSK;
IFQ = LOQ .* FSK;

% Remove high frequency mixing content
[b,a] = butter(5, Fs_dig/(Fs_matlab/2));
IFI_filt = filter(b,a,IFI);     % These should be in analog, dealt with by the mixers
IFQ_filt = filter(b,a,IFQ);

% Downsample to 40 MHz (Fs_dig) before further processing
Iout = IFI_filt(1:downrate:end);
Qout = IFQ_filt(1:downrate:end);
t = t(1:downrate:end);

%% Add IF Noise
% I_noiseless = Iout;
% Q_noiseless = Qout;
% Iout = awgn(Iout, SNR_dB, 'measured');
% Qout = awgn(Qout, SNR_dB, 'measured');

% To find power dB of input signal
% pow2db(rms(I_noiseless)^2)
% To verify snr
% Isnr = snr(I_noiseless, Iout-I_noiseless)
% Qsnr = snr(Q_noiseless, Qout-Q_noiseless)

%% Quantize to 5-bit resolution
Iout = Iout/max(abs(Iout));
Qout = Qout/max(abs(Qout));

Iout = floor(Iout*15);
Qout = floor(Qout*15);

Iout = max(Iout,-16);
Iout = min(Iout,15);
Qout = max(Qout,-16);
Qout = min(Qout,15);

%% 2-Sided FFT for checking purposes
X=fft(Iout);
X_shift=fftshift(X);
N=length(Iout);
w=fftshift((0:N-1)/N*2*pi);
w(1:N/2)=w(1:N/2)-2*pi;
f=Fs_dig/(2*pi)*w;
plot(f, abs(X_shift))
%% CSV Write
% This is the downconverted signal that needs to be filtered, and
% demodulated
csvwrite("I_out.csv",Iout);
csvwrite("Q_out.csv",Qout);
csvwrite("data_exp.csv",data_exp);

%% If you implemented I/Q compensation, it would go here

%% Filter signal with the complex bandpass filters

% Original
% Rcoeff = real(Hbp);
% Icoeff = imag(Hbp);

% Quantized
Rcoeff = Hbpreal;
Icoeff = Hbpimag;

x1 = conv(Iout, Rcoeff, 'same');
x2 = conv(Qout, Icoeff, 'same');

x3 = conv(Iout, Icoeff, 'same');
x4 = conv(Qout, Rcoeff, 'same');

Iout_prefilter = Iout;
Iout = x1 + x2;
Qout_prefilter = Qout;
Qout = x3 - x4;

% Plot the I and Q channels after filtering
plot(Iout);
figure;
plot(Qout);

%% 2-Sided FFT for checking purposes
X=fft(Iout);
X_shift=fftshift(X);
N=length(Iout);
w=fftshift((0:N-1)/N*2*pi);
w(1:N/2)=w(1:N/2)-2*pi;
f=Fs_dig/(2*pi)*w;
plot(f, abs(X_shift))

%% If we want to export post-filter signals, do so here

%% Downscale back to 5-bit resolution
scalingfactor = 2^11;           % This is set to a hard value so that when there's noise, the remaining signal is less
I = round(Iout/2^11);
Q = round(Qout/2^11);

%% 2-Sided FFT for checking purposes
X=fft(I);
X_shift=fftshift(X);
N=length(I);
w=fftshift((0:N-1)/N*2*pi);
w(1:N/2)=w(1:N/2)-2*pi;
f=Fs_dig/(2*pi)*w;
plot(f, abs(X_shift))
%% Matched Filter Demodulation

% Using these previous variables:
% T is symbol rate
% IF is IF frequency
% dT_dig is digital sampling period
tt = 0:dT_dig:T-dT_dig;

% m = round(T/dT_dig);
tau = 7*ones(1,num_bits);

% -1
template1Q_v = round(15*sin(2*pi*((IF-250e3)*tt)));
template1I_v = round(15*cos(2*pi*((IF-250e3)*tt)));

% +1
template2Q_v = round(15*sin(2*pi*((IF+250e3)*tt)));
template2I_v = round(15*cos(2*pi*((IF+250e3)*tt)));

% Matched Filter 1
I1_v = conv(I,template1I_v);
Q1_v = conv(I,template1Q_v);
I1_v = I1_v(0.5*T/dT_dig:end-0.5*T/dT_dig);
Q1_v = Q1_v(0.5*T/dT_dig:end-0.5*T/dT_dig);
mf1_v = sqrt(I1_v.^2 + Q1_v.^2);

% Matched Filter 2
I2_v = conv(I,template2I_v);
Q2_v = conv(I,template2Q_v);
I2_v = I2_v(0.5*T/dT_dig:end-0.5*T/dT_dig);
Q2_v = Q2_v(0.5*T/dT_dig:end-0.5*T/dT_dig);
mf2_v = sqrt(I2_v.^2 + Q2_v.^2);
compared_v = mf1_v < mf2_v;

data_v = zeros(1,num_bits);

for k = 2:num_bits-4
    idx = 40*k + tau(k);
    data_v(k) = compared_v(idx);
end

%% Plot, for debug
plot(t(1:1000), [mf1_v(1:1000)' mf2_v(1:1000)'])
legend("Matched Filter 1", "Matched Filter 2")
xlabel("Time (s)")
%% Compare with initial bitstream
data_orig = (data(2:end)+1)/2;
data_orig = [data_orig' 0]';
plot([data_v' data_orig])
legend('Demodulated', 'Original')
figure;plot(data_v'-data_orig);

%% Compare with output generated by chisel
data_chisel = csvread("data_out_demod.csv");
data_chisel = data_chisel(2:end);
data_exp = csvread("data_exp.csv");
data_exp = data_exp(1:end-1);
plot(data_chisel-data_exp)
%% CSV Write to export the templates

csvwrite("template1I.csv",template1I_v);
csvwrite("template1Q.csv",template1Q_v);
csvwrite("template2I.csv",template2I_v);
csvwrite("template2Q.csv",template2Q_v);