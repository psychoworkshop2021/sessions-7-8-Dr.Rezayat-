%% In His Exalted Name
% Psychophysic Workshope 2021  (Held in: Psychology & education school
% University of Tehran)
%
%
%written by : Helia Thaghavi
%
%
%heliabsb@gmail.com
%
%last edition : 2021/11/8
%

path_code = cd;
addpath(genpath([path_code '\sup']))

%% Priority
whichScreen = max(Screen('Screens'));
maxPriorityLevel = MaxPriority(whichScreen);
Priority(maxPriorityLevel);
if IsOSX
    fprintf('Your operating system is OSX and it''s Priority levels range from 0-9\n'); 
elseif IsWindows
    fprintf('Your operating system is Windows and it''s Priority levels range from 0-2\n');

end

close all;
clear all;
clc;

%% DAQ toolbox example
%
% %Find the DAQ device number
% devices = PsychHID('devices');
% USBdeviceNum = 0;
% DAQFound = 0;
% for i = 1:length(devices)
%     if strcmp(devices(i).product,'USB-1024LS')
%         USBdeviceNum = i;
%     end
% end
% %Let the user know if we've found it
% if (USBdeviceNum < 1)
% fprintf('ERROR: Could not locate USB device.\n');
% else
%     fprintf('FOUND DAQ DEVICE AT DEVICE %d',USBdeviceNum);
%     DAQFound = 1;
% end
% %Initialize USB box to zero
% if (DAQFound), Daq0ConfigPort(USBdeviceNum,4,0), end;
% %find key code for trigger key, which is a 5
% triggerCode =KbName('5%');
% keyIsDown = 0;
% %Make sure no keys are disabled
% DisableKeysForKbCheck([]);
% %wait for trigger
% while 1
% [ keyIsDown, pressedSecs, keyCode ] = KbCheck(-1);
% if keyIsDown
%     if find(keyCode)==triggerCode
%         break;
%     end
% end
% end
% %Send a pulse to the biopac
% if (DAQFound)
%     DaqDOut(USBdeviceNum,4,255);
%     DaqDOut(USBdeviceNum,4,0);
% end
% %Record trigger time for future reference
% triggerTime = pressedSecs;
% fprintf('Trigger detected\n');

%% IO Port Initialization
try
    fclose(instrfind)
catch
end

porta = serial('COM1', 'BaudRate', 57600, 'DataBits', 8);

fopen(porta);

% when you  want to send event
 eve = 1;    % a number in [0 :255]
 fwrite(porta,eve);

fclose ('all');

%% Shuffeling Trials  
nContrast = 7;
repeats   = 15;
% initialize trial numbers from 1 to nTrials
contrastIndex = repmat(1 : nContrast, repeats, 1); 
        % use repmat to cycle thru contrast #s
Cond_trilas  = Shuffle(contrastIndex(:)); 

Cond_trilas  = randsample(contrastIndex(:),nContrast*repeats); 

rng
%
randi(10,10,1)
randperm(10,5)
%
f=1:100;
Shuffle(f)

% Other Example codes
% RandSample()
% ChooseKFromN(100,200)
% RandSel()
% URandSel()
% CoinFlip()

%% Control timing of task
tic1= tic;
WaitSecs(3)
toc(tic1)

start_time = GetSecs();

present_time = GetSecs()-start_time;

%% Display Setup Module
% Define display parameters
whichScreen = max(Screen('screens'));
p.ScreenDistance = 30;  % in inches
p.ScreenHeight = 15;    % in inches
p.ScreenGamma = 2;  % from monitor calibration
p.maxLuminance = 100; % from monitor calibration
p.ScreenBackground = 0.5; 

% Open the display window, set up lookup table, and hide the 
% mouse cursor
if exist('onCleanup', 'class'), oC_Obj = onCleanup(@()sca); end  
        % close any pre-existing PTB Screen window
% Prepare setup of imaging pipeline for onscreen window. 
PsychImaging('PrepareConfiguration'); % First step in starting
                                      % pipeline
PsychImaging('AddTask', 'General',  'FloatingPoint32BitIfPossible');  
        % set up a 32-bit floatingpoint framebuffer
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');   
        % normalize the color range ([0, 1] corresponds 
        % to [min, max])
PsychImaging('AddTask', 'General', 'EnablePseudoGrayOutput'); 
        % enable high gray level resolution output with 
        % bitstealing
PsychImaging('AddTask','FinalFormatting',  'DisplayColorCorrection','SimpleGamma');  
        % setup Gamma correction method using simple power 
        % function for all color channels 
[windowPtr p.ScreenRect] = PsychImaging('OpenWindow', whichScreen, p.ScreenBackground );  
        % Finishes the setup phase for imaging pipeline
        % creates an onscreen window, performs all remaining  
        % configuration steps
PsychColorCorrection('SetEncodingGamma', windowPtr, 1/ p.ScreenGamma);  
% set Gamma for all color channels
HideCursor;  % Hide the mouse cursor 

% Get frame rate and set screen font
p.ScreenFrameRate = FrameRate(windowPtr); 
        % get current frame rate
Screen('TextFont', windowPtr, 'Times'); 
        % set the font for the screen to Times
Screen('TextSize', windowPtr, 24); % set the font size 
                                   % for the screen to 24
                                   
%%%%%%%%%%%%% Experimental Module %%%%%%%%%%%%

% Specify general experiment parameters

nTrials =   15  ; % number of trials
p.randSeed = ClockRandSeed; % use clock to set random number generator

% Specify the stimulus

p.stimSize = [8 6];      % horizontal and vertical stimulus size 
                         % in visual angle
p.stimDuration = 0.250;  % stimulus duration in seconds
p.ISI = 0.5;             % duration between response and next trial onset
p.contrast = 0.2;        % grating contrast
p.tf = 4;                % drifting temporal frequency in Hz
p.sf = 1;                % spatial frequency in cycles/degree

% grating contrast
% drifting temporal frequency in Hz
% spatial frequency in cycles/degree
% Compute stimulus parameters
ppd = pi/180 * p.ScreenDistance / p.ScreenHeight * ...
p.ScreenRect(4); % pixels per degree 
nFrames = round(p.stimDuration * p.ScreenFrameRate);
% stimulus frames
m = 2 * round(p.stimSize * ppd / 2);% horizontal and vertical
                                    % stimulus size in pixels 
sf = p.sf / ppd; % cycles per pixel
phasePerFrame = 360 * p.tf / p.ScreenFrameRate; % phase drift per frame
fixRect = CenterRect([0 0 1 1] * 8, p.ScreenRect); % 8 x 8 fixation
params = [0 sf p.contrast 0];
% parameters for DrawTexture: initial phase, spatial 
% frequency, contrast, 0
% procedural sinewavegrating on GPU sets up a texture for
% the sine wave
tex = CreateProceduralSineGrating(windowPtr, m(1), m(2), ...
[1 1 1 0] * 0.5, [], 0.5);


% CreateProceduralSineGrating is a psychtoolbox
% function that creates a procedural texture that
% allows you to draw sine grating stimulus patches in a very fast
% and efficient manner on modern graphics hardware.

% Initialize a table to set up experimental conditions

p.recLabel = {'trialIndex' 'motionDirection' 'respCorrect' ...
 'respTime'};   
rec = nan(nTrials, length(p.recLabel)); % matrix rec is nTrials x 4 of NaN
rec(:, 1) = 1 : nTrials;                % label the trial type numbers from 1 to nTrials
rec(:, 2) = -1;                         % -1 for left motion direction
rec(1 : nTrials/2, 2) = 1 ;             % half of the trials set to +1 for
                                        % right motion Direction 
rec(:, 2) = Shuffle(rec(:, 2));         % randomize motion direction over trials


% Prioritize display to optimize display timing 

Priority(MaxPriority(windowPtr));


% Start experiment with instructions

str = sprintf('Left/Right arrow keys for direction.\n\n Press SPACE to start.');
    
DrawFormattedText(windowPtr, str, 'center', 'center', 1);
% Draw instruction text string centered in window 
Screen('Flip', windowPtr);
% flip the text image into active buffer 
WaitTill('space'); % wait till space bar is pressed 
Screen('FillOval', windowPtr, 0, fixRect);
% create fixation box as black (0) 
Secs = Screen('Flip', windowPtr);
% flip the fixation image into active buffer 
p.start = datestr(now); % record start time

% Run nTrials trials 

for i = 1 : nTrials
    params(1) = 360 * rand;     % set initial phase randomly
    Screen('DrawTexture', windowPtr, tex, [], [], 0, ...
        [], [], [], [], [], params);
    % call to draw or compute the texture pointed to by tex
    % with texture parameters of the initial phase, the
    % spatial frequency, the contrast, and fillers required 
    % for 4 required auxiliary parameters
    t0 = Screen('Flip', windowPtr, Secs + p.ISI); % initiate first frame after p.ISI secs
    for j = 2 : nFrames % For each of the next frames one by one
        params(1) = params(1) - phasePerFrame * rec(i, 2);
        % change phase
        Screen('DrawTexture', windowPtr, tex, [], [], 0, ...
            [], [], [], [], [], params);
        % call to draw/compute the next frame
        Screen('Flip', windowPtr); % show frame
        % each new computation occurs fast enough to show % all nFrames at the framerate
    end
    Screen('FillOval', windowPtr, 0, fixRect); % black fixation for response interval
    Screen('Flip', windowPtr);
    [key Secs] = WaitTill({'left' 'right' 'esc'}); % wait till response
    if iscellstr(key), key = key{1}; end
    % take the 1st key in case of multiple key presses
    if strcmp(key, 'esc'), break; end
    % stop the trial sequence if keypress = <esc>
    respCorrect = strcmp(key, 'right') == (rec(i, 2) == 1); % compute if correct or incorrect
    rec(i, 3 : 4) = [respCorrect Secs-t0]; % record correctness and RT in rec
    if rec(i, 3), Beeper; end % beep if correct
end
p.finish = datestr(now); % record finish time

% Save Results
save DriftingSinewave_rst.mat rec p; % save the results
clear Screen; 

%% Plot results
clear all
load DriftingSinewave_rst.mat
cci = rec(:,2);
bhv = rec(:,3);
performance = [];
conditions = unique(cci);
for ci = 1: length(conditions)
    ind_h = ismember(rec(:,2),conditions(ci));
    
   performance(ci) = sum(ind_h&bhv)/sum(ind_h)*100 ;
   RT_cr (ci) = nanmean(rec(ind_h&bhv,4));
   RT_wr (ci) = nanmean(rec(ind_h&(~bhv),4));

end
figure
hold on
bar(conditions,performance)
xlabel('Contrast (a.u.)')
ylabel('directions (%)')

figure
hold on
bar(conditions,RT_cr,'g')
bar(conditions,RT_wr,'b')
xlabel('directions (a.u.)')
ylabel('RT (sec.)')

%% play sound

%Mark script start time
scriptStart = GetSecs();
%Initialize the sound driver
InitializePsychSound;
%Read in the sound data from file
wavfilename=[cd '/sup/Applause.wav']
[soundData, freq] = audioread(wavfilename);
%Prepare sound data (make it two rows for stereo playback)
soundData = [soundData, soundData];
numChannels=2;
%Open the audio driver
pahandle = PsychPortAudio('Open',[], [], 0, freq,numChannels);
%Fill the buffer
PsychPortAudio('FillBuffer',pahandle,soundData');
%Play the sound
playTime = PsychPortAudio('Start',pahandle);
fprintf('\nSound started playing %.2f seconds after start of script\n',playTime,scriptStart);
%Close the audio driver
PsychPortAudio('Stop', pahandle, 1,0);
PsychPortAudio('Close',pahandle);

%% recording sound

%Initialize sound driver
InitializePsychSound;
duration = 5;
%Open audio channel for recording using mode 2
freq = 44100;
pahandle = PsychPortAudio('Open', [], 2, 0, freq, 2);
%Set up buffer for recording
PsychPortAudio('GetAudioData', pahandle, duration);
%Start recording
PsychPortAudio('Start', pahandle, 0, 0, 1);
%Go until keypress
fprintf('Recording...\n');
WaitSecs(duration);
fprintf('Done recording.\n');
%Stop Recording
PsychPortAudio('Stop', pahandle);
%Get the audio data we recorded
audioData = PsychPortAudio('GetAudioData',pahandle);

%Close the audio channel
PsychPortAudio('Close',pahandle);
