%% 0. Program Info
% Program Author: Mangesh Sonawane
% Program Details: A basic implementation of the Phase Vocoder as described
% in the MAFTDSP assignment 1 instructions.

%% 1 - Variable Declaration 
% Variables to be defined by the user

% Define file to be manipulated, TEST_FILE, Window Size in sample, N, 
% Analysis Hop Size in samples, HA, Time-Strech Factor, Q.

TEST_FILE = 'Mozart.wav';               
N = 1024;
HA = N/4;                               
Q = 1.5;
if Q <= 0 
    error('Error. Time stretch factor Q must be a positve number')
end
if N < 512 || N > 2048
    error('Error. Window size N must be between 512 and 2048')
end
if (N/HA) < 4
    error('Error. Hop Analysis size HA must be at most %25 of window size')
end
if rem(N,1) || rem(HA,1)
    error('Error. Window Size and Hop Size must be integers')
end

%% 2 - Determine Synthesis Hop Size

% Calculate Synthesis Hop Size in samples, HS, and round to integer.
HS = round(Q*HA);                       
% Reset Time-Strech Factor
Q = HS/HA;

%% 3 - Generate Hann Window

% Generate periodic Hann Window, win.
n = transpose(0:N-1);
win = 0.5*(1-cos(2*pi*n/N));    

%% 4 - Create Normalised Angular Frequency Storage Matrix

% The normalised frequency associated with each of the N FFT bins, bins 
% from k = 0 to k = N-1;
% The bin values for 0 < k <= N/2 (or k <= (N-1)/2 if N is odd) represent
% positive frequencies
posFreq = (1 : (N - mod(N,2)) / 2)';
% The bin values for k > N/2 (or k > (N-1)/2 if N is odd) represent 
% negative frequencies.
negFreq = (-(N + mod(N,2))/2 + 1  : -1)';
% omega contains the normalised frequency associated with each bin.
omega = [0; 2*pi*posFreq/N ; 2*pi*negFreq/N];

%% 5 - Read Input File

% Read wav file and assign the samples to 'x' and the sample rate to 'Fs'
try 
   [x,Fs] = audioread(TEST_FILE);
catch ME
    if (strcmp(ME.identifier,'MATLAB:audiovideo:audioread:fileNotFound'))
        error('Error. Wav file not found.')
    end
end
% Checks for stereo files and converts to mono
if size(x,2) == 2
    x = 0.5*sum(x,2);
end

%% 6 - Zero Pad Input

% Zero pad x at start with N zeros, so the first window's phase is 0
x = [zeros(N,1) ; x];
% Assign the number of samples in x to 'L'
L = length(x);
if mod(L-N,HA) ~= 0
    % Zero pad x if needed so a whole number of analysis frames fit into x
    Diff = HA - mod(L-N,HA);
    x = [x ; zeros(Diff, 1)];
    L = length(x);
end
% Find number the number of frames required to cover the zero padded signal
NF = (L-N)/HA + 1;

%% 7 - Initialise Matrices

% X will contain the frames of the input signal, each N samples long
X = zeros(N,NF);
% instFreq will contain the instantaneous frequencies for each frame
instFreq = zeros(N,NF);
% thetaMat will contain the modified phase values for each frame
thetaMat = zeros(N,NF);
% YF will contain the frequency domain frames of the output signal
YF = zeros(N,NF);
% Y will contain the frames of the output signal
Y = zeros(N,NF);

%% 8 - Create Output Vector

% y will contain the output, and will be constructed from the Y matrix
y = zeros(N + (NF-1)*HS ,1);

%% 9 - Load X

% Load each analysis frame of the input into the columns of X and window
for c = 0:NF-1
    X(:,c+1) = win.*x((c*HA)+1:(c*HA)+N,1);
end

%% 10 - Compute DFT of X

% XF is filled with the fourier transform for each of the frames in X
XF = zeros(N,NF);
for c = 0:NF-1
    XF(:,c+1) = fft(X(:,c+1));
end
% Store Magnitude and Phase of XF in XFM and XFP
XFM = abs(XF);
XFP = angle(XF);

%% 11 - Compute Phase Modification

for c = 2:NF-1
    % The Instantaneous Frequencies are calculated
    instFreq(:,c) = XFP(:,c) - XFP(:,c-1) - omega*HA;
    % Phase wrapping to find equivalent phases in the range of -pi to pi is 
    % done by dividing by 2pi and rounding to the nearest integer, this 
    % integer is multiplied by 2pi and subtracted from the phase.
    instFreq(:,c) = instFreq(:,c) - round(instFreq(:,c)/(2*pi))*2*pi;
    instFreq(:,c) = instFreq(:,c)./HA + omega;
    % The phase angle adjustment is calculated from the Instantaneous
    % Frequencies and the Synthesis Hop Size
    thetaMat(:,c) = thetaMat(:,c-1) + HS.*instFreq(:,c);
    % The phase at the k=0 (and k=N/2 if N is even) frequency bins are
    % reset to their initial values, they should both be purely real.
    % If negative they are read as having a phase angle of pi, which
    % results in undesired moficiation to their phases.
    thetaMat(1,c) = XFP(1,c);
    if mod(N,2) == 0
        thetaMat(N/2 + 1,c) = XFP(N/2 + 1,c);
    end
    % The modified phases are combined with the initial fourier magnitudes
    YF(:,c) = XFM(:,c).*exp(1i*thetaMat(:,c));
    % The modified frames have the IFFT applied and then are windowed
    Y(:,c) = win.*ifft(YF(:,c));
    % The output frames are then overlap added into the output vector
    y(((c-2)*HS)+1:((c-2)*HS)+N,1) = ...
        y(((c-2)*HS)+1:((c-2)*HS)+N,1) + Y(:,c);
end
% Seperate the rounding error from the output so it is purely real
y = real(y);
% rms of the 2 input/output is compared to determine required gain
gain = rms(x)/rms(y);
y = y*gain;
soundsc(y,Fs)

%% 12 - Plot Data

% Plot 1 - The Input Waveform
time_axis_in = (0:1/Fs:(L-1)/44100);
fig1 = figure(1);
set(fig1, 'Position', get(0,'Screensize'));
ax(1) = subplot(2,2,1);
plot(time_axis_in,x)
if Q > 1
    xlim([0,ceil(length(y)/Fs)])
end
title('Input Signal') ; xlabel('Time (s)')
hold on

% Plot 2 - The STFT Frame Magnitudes
ax(2) = subplot(2,2,2);
imagesc(XFM(1:ceil(N/2),:), [0 4])
 %ax(2) = gca;
 set(gca,'YDir','normal')
% xticks(0:round(NF/8):NF);
% yticks(0:round(N/16):N/2);
% yticklabels(Fs*yticks/N)
title('STFT Frame Magnitudes vs Time and Frequency')
xlabel('Frame Number') ; ylabel('Frequency(Hz)')

% Plot 3 - The Output Waveform
ax(3) = subplot(2,2,3);
time_axis_out = (0:1/Fs:(length(y)-1)/44100);
plot(time_axis_out,y)
if Q < 1
    xlim([0,ceil(length(x)/Fs)])
end
title('Output Signal')
xlabel('Time (s)')

% Plot 4 - The Instantaneous Frequencies
ax(4) = subplot(2,2,4);
imagesc(instFreq)
ax(4) = gca;
set(gca,'YDir','normal')
xticks(0:round(NF/8):NF);
yticks(0:round(N/4):N);
title('Instantaneous Frequency vs Frequency Bin and Time')
xlabel('Frame Number') ; ylabel('Frequency Bin (k)')
ax(5) = colorbar('Ticks',-pi:pi/2:pi,...
    'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'});
set(ax, 'FontSize', 20)
hold off
