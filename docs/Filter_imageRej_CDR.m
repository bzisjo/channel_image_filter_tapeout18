close all
clear all
clc

format shortg
clock1 = clock;
%% Some parameters
Fs_dig = 40e6;       % 40MHz: clock frequency of the digital module
dT_dig = 1/Fs_dig;
downrate = 200;      % This is used for the MATLAB filters
Fs_matlab = downrate * Fs_dig;
dT_matlab = 1/Fs_matlab;
RF_ampl = 2;         
LOI_ampl = 1;        % LO ampl for In-phase arm
LOQ_ampl = 1;        % LO ampl for Quadrature arm
RF_freq = 2.406e9;   % 2.406 GHz: RF frequency. Refer to BLE specs, channel 1
IF = 2.5e6;          % Intermediate frequency
SNR_dB = 6;          % Analog SNR spec requirement 10 dB; Here gives some margins

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

%% Quantize the filter to 5-bits
scalingfactorreal = min(abs(-16/min(real(Hbp))), abs(15/max(real(Hbp))));
scalingfactorimag = min(abs(-16/min(imag(Hbp))), abs(15/max(imag(Hbp))));
maxscalingfactor = min(scalingfactorreal, scalingfactorimag);
Hbpreal = round(real(Hbp) * maxscalingfactor);
Hbpimag = round(imag(Hbp) * maxscalingfactor);

%% Export filter coeffs
% These are the image rejection/channel selection filter coefficients
csvwrite("rcoeffs.csv", Hbpreal)
csvwrite("icoeffs.csv", Hbpimag)

%% Generate Data and modulate as RF MSK
num_bits = 500;
T=1e-6;             % This is the symbol rate
freq_shift = 250e3; % This is the frequency deviation or tone spacing for fsk

% Generate some random data
data = randi([0 1],num_bits-8,1);
data = [1 0 1 0 1 0 1 0 data']';    % Add preamble
data_exp = data;    % Data expected at the end of CDR block
data = data*2-1;    % Convert data to +/- 1

num_sampling = downrate*T/dT_dig;  % #(sampling times for each bit in Matlab)
phase = zeros(1,num_sampling*num_bits);
% Generate first cycle of RF
t = 0:1/Fs_matlab:T-1/Fs_matlab;
phase(1:num_sampling) = 2*pi*(RF_freq+data(1)*freq_shift)*t;
phase_last = phase(num_sampling);

% Generate remaining modulated RF cycles
t = 1/Fs_matlab:1/Fs_matlab:T;
for x=2:num_bits
    % The phase of each cycle starts where the last one ended
    phase((x-1)*num_sampling+1:(x*num_sampling)) = 2*pi*(RF_freq+data(x)*freq_shift)*t + phase_last;
    phase_last = phase(x*num_sampling);
end


% Scale amplitude to get desed power
FSK = RF_ampl * sin(phase);
FSK_pure = FSK;


%% Generate data and modulate as RF side channel (dependent on prev block)
sidechannel_freq = 2e6;
data_side = randi([0 1],num_bits,1);
data_side = data_side*2-1;

phase = zeros(1,num_sampling*num_bits);
t = 0:1/Fs_matlab:T-1/Fs_matlab;
phase(1:num_sampling) = 2*pi*((RF_freq+sidechannel_freq)+data_side(1)*freq_shift)*t;
phase_last = phase(num_sampling);

t = 1/Fs_matlab:1/Fs_matlab:T;
for x=2:num_bits
    phase((x-1)*num_sampling+1:(x*num_sampling)) = 2*pi*((RF_freq+sidechannel_freq)+data_side(x)*freq_shift)*t + phase_last;
    phase_last = phase(x*num_sampling);
end

side = RF_ampl * sin(phase);
side = side(1:length(FSK_pure));
FSK = FSK + side;


%% Generate data and modulate as RF image (dependent on prev block)
num_bits_im = ceil(num_bits*1.5);
RF_freq_image = LO_freq-IF;
data_image = randi([0 1],num_bits_im,1);
data_image = data_image*2-1;

phase = zeros(1,num_sampling*num_bits_im);
t = 0:1/Fs_matlab:T-1/Fs_matlab;
phase(1:num_sampling) = 2*pi*((RF_freq_image)+data_image(1)*freq_shift)*t;
phase_last = phase(num_sampling);

t = 1/Fs_matlab:1/Fs_matlab:T;
for x=2:num_bits_im
    phase((x-1)*num_sampling+1:(x*num_sampling)) = 2*pi*((RF_freq_image)+data_image(x)*freq_shift)*t + phase_last;
    phase_last = phase(x*num_sampling);
end

image = RF_ampl * sin(phase);
image = image(1:length(FSK_pure));
FSK = FSK + image;

%% Mix and match zone
% FSK_pure, side, image can be added at will to form waveform-of-interest

%% Add White Gaussian Noise
FSK_noiseless = FSK;
FSK = awgn(FSK, SNR_dB, 'measured');

% To verify snr, for debug
FSKsnr = snr(FSK_noiseless, FSK-FSK_noiseless);
figure; plot(t(1:2000), [FSK(1:2000)' FSK_noiseless(1:2000)']);
legend("Noisy Signal", "Original Signal")
xlabel("Time (s)")
ylabel("Amplitude")


%% Create LO
t = 0:dT_matlab:(length(FSK)-1)*dT_matlab;
LOI = LOI_ampl * sin(2*pi*LO_freq*t) + LOI_ampl;        % In phase: carrier frequency
LOQ = LOQ_ampl * sin(2*pi*LO_freq*t - pi/2) + LOQ_ampl; % Quadrature: carrier frequency - pi/2

% Downconvert with LO
IFI = LOI .* FSK;
IFQ = LOQ .* FSK;

% Remove high frequency mixing content
[b,a] = butter(5, Fs_dig/(Fs_matlab/2));
IFI_filt = filter(b,a,IFI);     % These should be dealt with by the mixers
IFQ_filt = filter(b,a,IFQ);

% Downsample to 40 MHz (Fs_dig) before further processing
Iout = IFI_filt(1:downrate:end);
Qout = IFQ_filt(1:downrate:end);
t = t(1:downrate:end);


%% Add IF Noise
I_noiseless = Iout;
Q_noiseless = Qout;
Iout = awgn(Iout, SNR_dB, 'measured');
Qout = awgn(Qout, SNR_dB, 'measured');

% To find power dB of input signal
pow2db(rms(I_noiseless)^2);
% To verify snr
Isnr = snr(I_noiseless, Iout-I_noiseless);
Qsnr = snr(Q_noiseless, Qout-Q_noiseless);

%% Quantize to 5-bit resolution (ADC)
Iout = Iout/max(abs(Iout)); % Normalize
Qout = Qout/max(abs(Qout));

Iout = floor(Iout*15);      % 5-bit quantization: [-16,15]
Qout = floor(Qout*15);

% Iout = max(Iout,-16);
% Iout = min(Iout,15);
% Qout = max(Qout,-16);
% Qout = min(Qout,15);      % Redundant; ensure the data range is [-16,15]

%% CSV Write
% Downconverted signal that needs to be filtered and demodulated
csvwrite("I_out_withImage.csv",Iout);
csvwrite("Q_out_withImage.csv",Qout);
csvwrite("data_exp.csv",data_exp);

%% 2-Sided FFT for checking purposes
Xi=fft(Iout);
Xi_shift_before=fftshift(Xi);
Xq=fft(Qout);
Xq_shift=fftshift(Xq);
N=length(Iout);
w=fftshift((0:N-1)/N*2*pi);
w(1:N/2)=w(1:N/2)-2*pi;
f=Fs_dig/(2*pi)*w;
figure; plot(f, [abs(Xi_shift_before)' abs(Xq_shift)']);
legend("Iout", "Qout")

%figure; plot(Iout);
%% If you implemented I/Q compensation, it would go here

%% Filter signal with the complex bandpass filters: Image Rejection
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

%figure; plot(Iout);

%% 2-Sided FFT for checking purposes
Xi=fft(Iout);
Xi_shift_after=fftshift(Xi);
Xq=fft(Qout);
Xq_shift=fftshift(Xq);
N=length(Iout);
w=fftshift((0:N-1)/N*2*pi);
w(1:N/2)=w(1:N/2)-2*pi;
f=Fs_dig/(2*pi)*w;
figure; plot(f, [abs(Xi_shift_after)' abs(Xq_shift)']);
legend("Iout", "Qout")

scalingfactor = 2^8;
figure; plot(f, [abs(Xi_shift_before)' abs(Xi_shift_after/scalingfactor)']);
legend("Iout before", "Iout after")
%% If we want to export post-filter(image rejection) signals, do so here

%% Downscale back to 5-bit resolution
% scalingfactor = 2^11;         % This is set to a hard value so that when there's noise, the remaining signal is less
scalingfactor = 2^8;            % 5-bit filters
I = round(Iout/scalingfactor);
Q = round(Qout/scalingfactor);

%% Data Recovery: Matched Filter Demodulation
% T is symbol rate, 1Mbps
% IF is IF frequency
% dT_dig is digital sampling period, 1/40MHz
tt = 0:dT_dig:T-dT_dig;

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

%figure; plot(compared_v(1:2000));
%% Plot, for debug
% plot(t(1:2000), [mf1_v(1:2000)' mf2_v(1:2000)']);
% legend("Matched Filter 1", "Matched Filter 2")
% xlabel("Time (s)")

%% CSV Write to export the templates
csvwrite("template1I.csv",template1I_v);
csvwrite("template1Q.csv",template1Q_v);
csvwrite("template2I.csv",template2I_v);
csvwrite("template2Q.csv",template2Q_v);

%% Clock Recovery: Multibit Shift Register
% Use the proper index to simulate multibit shift register

compared_v = compared_v*2 - 1;  % mapping to +/- 1

% If here would like to add white noise for evaluation
% SNR_db = 10;
% compared_v = awgn(compared_v, SNR_db);

data_v = zeros(1, num_bits);

ratio = T/dT_dig;   % the convert ratio of clock recovery
                    % Here made one decision every 40 bits
                    
ptr_pst = 0;       % Current ptr position, the start of sum
ptr_step = 2;       % Step size of ptr position tuning
idx = 1;

while idx < num_bits
    data_sum = sum(compared_v((idx-1)*ratio+1+ptr_pst : idx*ratio+ptr_pst));
    data_first_half = sum(compared_v((idx-1)*ratio+1+ptr_pst : (idx-0.5)*ratio+ptr_pst));
    data_second_half = data_sum - data_first_half;
    ptr_flag = 0;   % -1 moves backwards; 0 stays; 1 moves forwards;
    
    % Determine the output based on data_sum
    % Determine whether ptr moves forward or backward
    if data_sum == ratio || data_sum == -ratio
        data_v(idx) = data_sum > 0;
    elseif data_sum > 0
        data_v(idx) = 1;
        if data_first_half > data_second_half  
            ptr_flag = -1;
        else
            ptr_flag = 1;
        end    
    else
        data_v(idx) = 0;
        if data_first_half < data_second_half  
            ptr_flag = -1;
        else
            ptr_flag = 1;
        end
        
    end
    
    % ptr position tuning based on ptr_flag and ptr_step
    if ptr_flag == 1        % ptr moves forward
        ptr_pst = ptr_pst + ptr_step;
    elseif ptr_flag == -1   % ptr moves backward
        ptr_pst = ptr_pst - ptr_step;
    % else ptr_flag == 0
    end
    
    idx = idx + 1;
end


error_rate_model = mean(abs(data_v'-data_exp));


%% Compare with initial bitstream
figure; plot(data_v'-data_exp);
ylim([-1 1])
title("Matlab model: Demodulated - Original");
txt = ['BER: ' num2str(error_rate_model)];
text(0.75*num_bits,0.8,txt,'FontSize', 14);

%% for debug
format shortg
exec_time = clock-clock1;
fprintf("Execution time: %d hours %d minutes %f seconds\n", [exec_time(4)  exec_time(5)  exec_time(6)]);

%% Compare with output generated by chisel
% data_chisel = csvread("data_out_demod_mfcdr.csv");
% data_chisel = data_chisel(2:end);
% data_exp = csvread("data_exp.csv");
% data_exp = data_exp(1:end);
% error_rate_chisel = mean(abs(data_chisel-data_exp));
% figure; plot(data_chisel-data_exp);
% ylim([-1 1])
% title("chisel: Demodulated - Original");
% txt = ['BER: ' num2str(error_rate_chisel)];
% text(0.75*length(data_exp),0.8,txt,'FontSize', 14);



