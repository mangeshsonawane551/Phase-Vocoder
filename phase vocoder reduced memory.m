%% 0. Program Info

% Program Details: An implementation of the Phase Vocoder 
% which doesn't utilise any large matrices and
% reduce memory usage.

%% 1 - Variable Declaration 

% Define file to be manipulated, TEST_FILE, Window Size in sample, N, 
% Analysis Hop Size in samples, HA, Time-Strech Factor, Q.

TEST_FILE = 'Mozart.wav';               
N = 2^10;                               
HA = N/4;                               
Q = 2;
if Q <= 0 
    error('Error. Input must be a positve number')
end

%% 2 - Determine Synthesis Hop Size

% Calculate Synthesis Hop Size in samples, HS, and round to integer.
HS = round(Q*HA);                       
% Reset Time-Strech Factor
Q = HS/HA;
% Determine number of unique values required to reconstruct signal using
% conjugate symmetry
uniques = floor(N/2) + 1;

%% 3 - Generate Hann Window

% Generate periodic Hann Window, win.
n = transpose(0:N-1);
win = 0.5*(1-cos(2*pi*n/N));    

%% 4 - Create Normalised Angular Frequency Storage Matrix

% The bin values for 0 < k <= N/2 (or k <= (N-1)/2 if N is odd) represent
% positive frequencies, negative frequencies are not required due to
% conjugate symmetry
posFreq = (1 : (N - mod(N,2)) / 2)';
% omega contains the normalised frequency associated with each bin up to 
% the Nyquist bin.
omega = [0; 2*pi*posFreq/N];

%% 5 - Read Input File

% Read wav file and assign the samples to 'x' and the sample rate to 'Fs'
[x,Fs] = audioread('Mozart.wav');
% Checks for stereo files and converts to mono
if size(x,2) == 2
    x = 0.5*sum(x,2);
end

%% 6 - Zero Pad Input

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

%% 7 - Create Output Vector

% y will contain the output, and will be constructed from the Y matrix
y = zeros(N + (NF-1)*HS ,1);

%% 8 - Carry Out Frame by Frame Phase Modification

% Set initial value of zero for previous frame phase and modified phase
prevXFP = zeros(uniques,1);
prevThetaMat = zeros(uniques,1);

% Load each analysis frame of the input in turn and make phase modification
for c = 1:NF-1
    % window input
    X = win.*x((c*HA)+1:(c*HA)+N);
    % compute fft of frame and seperate phase/magnitude
    XF = fft(X);
    XFM = abs(XF);
    XFP = angle(XF);
    % throw away the negative frequencies of the phases
    XFP = XFP(1:uniques);
    % calculate the instantaneous frequencies for the frame
    instFreq = XFP - prevXFP - omega*HA;
    % phase wrap to find equivalent phases in the range of -pi to pi. Done 
    % by dividing by 2pi and rounding to the nearest integer, this integer 
    % is multiplied by 2pi and subtracted from the phase.
    instFreq = instFreq - round(instFreq/(2*pi))*2*pi;
    instFreq = instFreq./HA + omega;
    % save the phase angles to prevXFP for use in the next iteration
    prevXFP = XFP;
    % calculate the required phase modification
    thetaMat = prevThetaMat + HS.*instFreq;
    % The phase at the k=0 (and k=N/2 if N is even) frequency bins are
    % reset to their initial values, they should both be purely real.
    % If negative they are read as having a phase angle of pi, which
    % results in undesired moficiation to their phases.
    thetaMat(1) = XFP(1);
    if mod(N,2) == 0
        thetaMat(N/2 + 1) = XFP(N/2 + 1);
    end
    % save the phase modification to prevThetaMat for the next iteration
    prevThetaMat = thetaMat;
    % construct output using the input magnitudes and the modified phases
    YF = XFM.*exp(1i*[thetaMat;-flipud(thetaMat(2:end-1+mod(N,2)))]);
    Y = win.*ifft(YF);
    % overlap add the output frame into the output vector
    y((c-1)*HS+1:(c-1)*HS+N) = y((c-1)*HS+1:(c-1)*HS+N) + Y;
end

y = real(y);
% rms of the 2 input/output is compared to determine required gain
gain = rms(x)/rms(y);
y = y*gain;
soundsc(y,Fs)
