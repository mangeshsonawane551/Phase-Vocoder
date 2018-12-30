%% 0. Program Info
% Program Author: Mangesh Sonawane
% Program Details: A GUI for the Phase Vocoder 


%% 1 - Create Layout
% Create the GUI figure, user controls and output plots.
close all
% Create main figure which will contain the GUI
handles.main = figure('Name', 'Phase Vocoder','NumberTitle', 'off',...
    'Position', [184   176   919   508]);
% Set the default value for Q and save a list of the files in the pwd
setappdata(gcf,'Q',1.0);
pwd_files = dir('*.wav');
% Create play input button
handles.play_in = uicontrol('style','pushbutton','string','Play Input File',...
    'units','normalized','position',[0.25 0.54 0.2 0.16],...
    'Enable','off','callback', {@play,'input'});
% Create list of wav files in pwd also with option to select a file
handles.file_menu = uicontrol('style','listbox','units','normalized',...
    'string',{pwd_files.name,'select file'},...
    'position',[0.05 0.545 0.2 0.15], 'callback',@load);
% Create label for the stretch factor, 'Q'
handles.Q_label_h = uicontrol('style','text','units','normalized',...
    'string','Time Stretch Factor','position',[0.4 0.85 0.2 0.1],...
    'fontsize',14);
% Create readout for the stretch factor, 'Q'
handles.Q_readout = uicontrol('style','text','units','normalized',...
    'string','1','position',[0.475 0.87 0.05 0.04],...
    'fontsize',14,'fontweight','bold');
% Create slider for setting the stretch factor, 'Q'
handles.strch_fac = uicontrol('style','slider','units','normalized','max',8,...
    'min',0.1,'value',1,'units','normalized',...
    'position',[0.4 0.72 0.2 0.1],'callback',@Q_Upd);
% Create button to process the audio 
handles.create_out =  uicontrol('style','pushbutton','string','Process Audio',...
    'units','normalized','position',[0.55 0.54 0.2 0.16],...
    'Enable','off','callback', @process);
% Create button to play the processed audio 
handles.play_out = uicontrol('style','pushbutton','string','Play Output File',...
    'units','normalized','position',[0.75 0.54 0.2 0.16],...
    'Enable','off','callback', {@play,'output'});
% Create axes which display the input/output signal
handles.ax(1) = subplot('position',[0.05 0.25 0.4 0.2],...
    'nextplot','replacechildren');
title('Input Signal:');xlabel('time(s)');

handles.ax(2) = subplot('position',[0.55 0.25 0.4 0.2],...
    'nextplot','replacechildren');
title('Output Signal:');xlabel('time(s)');
set(handles.ax,'fontsize', 15, 'ylim',[-1,1])
% Save the handles to guidata so they can be accessed by callback functions
guidata(handles.main,handles)
%% 2 - Create Functions
% Create the callback functions for the GUI's buttons and slider.

%% 'load' reads the chosen audio file or opens a dialog for file selection
function load(self_handle, ~ )
    % Get handles for other buttons/plots
    handles = guidata(self_handle);
    ax1 = handles.ax(1); 
    play_in = handles.play_in; 
    create_out = handles.create_out;
    play_out = handles.play_out;
    % Create string containing the selected option
    file_list = get(self_handle,'string');
    file_index = get(self_handle,'value');
    file_name = char(file_list(file_index));
    % Define default empty path
    path = '';
    % If 'select file' is chosen open dialog for user to choose file
    if strcmp(strtrim(file_name),'select file') == 1  
        [file_name,path] = uigetfile...
            ('../*.wav','Select the wav audio file');
    end
    % Read audio file
    [x,Fs] = audioread([path,file_name]);
    % Convert to mono if input is stereo
    if size(x,2) == 2
        x = 0.5*sum(x,2);
    end
    % Save input to app data
    setappdata(gcf,'input',x)
    setappdata(gcf,'Fs',Fs)
    setappdata(gcf,'file_name',file_name)
    % Plot input signal
    time_axis_in = (0:1/Fs:(length(x) - 1)/Fs)';
    plot(ax1,time_axis_in,x); 
    title(ax1,['Input Signal: ',file_name],'interpreter', 'none');
    % Set pushbuttons on/off
    set(play_in,'Enable','on')
    set(create_out,'Enable','on')
    set(play_out,'Enable','off')
end   

%% 'play' plays the input/output signal depending on which button called it
function play(self_handle,~,input)
    % Get all handles
    handles = guidata(self_handle);
    % Get all handles which can be enabled (buttons/sliders etc)
    uicontrols = findobj(handles.main,'Enable','on');
    % Load the signal to be played and the sample rate
    signal = getappdata(gcf,input);
    Fs = getappdata(gcf,'Fs');
    % Play the signal
    soundsc(signal,Fs);
    % Disable other buttons whilst playing
    set(uicontrols,'Enable','off');
    pause(length(signal)/Fs);
    % Try statement preventing an error being thrown if user closes the GUI
    % while the audio is playing
    try
        set(uicontrols,'Enable','on');
    catch ME
        if (strcmp(ME.identifier,'MATLAB:class:InvalidHandle'))
        end
    end
end

%% 'Q_Upd' updates the guidata for stretch factor as well as the readout
function Q_Upd(self_handle,~)
    % Get required handles
    handles = guidata(self_handle);
    Q_readout = handles.Q_readout; 
    % Get chosen stretch factor, 'Q', then round to 2 dp.
    Q = get(self_handle,'Value');
    Q = round(Q,2);
    % Set chosen Q to the guidata and update readout
    setappdata(gcf,'Q',Q);
    set(Q_readout,'string',Q);
end


%% 'process' carries out the time stretching
function process(self_handle,~)
%% 1 - Initialise
%Get the required handles
handles = guidata(self_handle);
play_out = handles.play_out; 
ax2 = handles.ax(2); 
% Get stretch factor and initialise constants
Q = getappdata(gcf,'Q');
N = 2^10;
% Set HA based on stretch factor chosen
if Q >= 0.1 && Q <= 1.7
    HA = N/4;
elseif Q > 1.7 && Q <= 2.4
    HA = N/8;
elseif Q > 2.5 && Q <= 5
    HA = N/16;
else
    HA = N/32;
end
HS = round(Q*HA);
% Determine number of unique values required to reconstruct signal using
% conjugate symmetry
uniques = floor(N/2) + 1;

%% 2 - Generate Hann Window
n = transpose(0:N-1);
win = 0.5*(1-cos(2*pi*n/N));    

%% 3 - Create Normalised Angular Frequency Storage Matrix

% The bin values for 0 < k <= N/2 (or k <= (N-1)/2 if N is odd) represent
% positive frequencies, negative frequencies are not required due to
% conjugate symmetry
posFreq = (1 : (N - mod(N,2)) / 2)';
% omega contains the normalised frequency associated with each bin up to 
% the Nyquist bin.
omega = [0; 2*pi*posFreq/N];

%% 4 - Get Input

x = getappdata(gcf,'input');
Fs = getappdata(gcf,'Fs');

%% 5 - Zero Pad Input

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

%% 6 - Create Output Vector

% y will contain the output, and will be constructed from the Y matrix
y = zeros(N + (NF-1)*HS ,1);

%% 7 - Compute Phase Modification

% Set initial value of zero for previous frame phase and modified phase
prevXFP = zeros(uniques,1);
prevThetaMat = zeros(uniques,1);

% Load each analysis frame of the input in turn and make phase modification
for c = 1:NF-1
    % window input
    X = win.*x((c*HA)+1:(c*HA)+N);
    % compute fft of frame and seperate phase/magnitude
    XF = fft(X);
    % throw away the negative frequencies of the phases
    XF = XF(1:uniques);
    XFM = abs(XF);
    XFP = angle(XF);
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
    % save the phase modification to prevThetaMat for the next iteration
    prevThetaMat = thetaMat;
    % construct output using the input magnitudes and the modified phases
    YF = XFM.*exp(1i*thetaMat);
    % IFFT and window the newly phase modified output, using the symmetric
    % flag so MATLAB treats YF as conjugate symmetric
    Y = win.*ifft(YF,N,'symmetric');
    % overlap add the output frame into the output vector
    y((c-1)*HS+1:(c-1)*HS+N) = y((c-1)*HS+1:(c-1)*HS+N) + Y;
end

% Seperate the rounding error from the output so it is purely real
y = real(y);
% rms of the 2 input/output is compared to determine required gain
gain = rms(x)/rms(y);
y = y*gain;
% Save output signal to app data
setappdata(gcf,'output',y)
% Plot output signal
time_axis_in = (0:1/Fs:(length(y) - 1)/Fs)';
plot(ax2,time_axis_in,y); 
title(ax2,['Output Signal: ',getappdata(gcf,'file_name')],'interpreter', 'none');
% Enable play output button
set(play_out,'Enable','on')
end

