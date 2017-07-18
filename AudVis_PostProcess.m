function AudVis_PostProcess
%% File selection

% Select files
[task_file, tpath] = uigetfile('*.mat', 'Task performance data');
[stim_file, spath] = uigetfile('*.mat', 'Stimulus data');

% Load files
task = load([tpath task_file]);
stim = load([spath stim_file]);

% Data in task file should be of format: 
    % Variable name:  task_table
    % Columns: Trial #, Visual Mode, Coherence, Target Speaker, Choice,
    %          Success, Stimulus Start Time, Stimulus Stop Time, Response 
    %          Time Stamp, RT
    % Also contains meta_data:
        % subject, date, lowFreq, hiFreq, toneDur, trialDur, fs, nTrials

% Data in stim file should be of format:
    % Variable name:  stim_table
    % Columns: Trial #, Waveform, Masker, Frequency Vector, isHigh

%% First must calculate: played coh, played hi or lo, and num tones played

nTrials = meta_data.nTrials;

fs = meta_data.fs;
hiFreq = meta_data.hiFreq;
loFreq = meta_data.loFreq;

delay = 1600; %ms
delaysamp = floor(delay / 1000 * fs);

isH_played = zeros(nTrials,1);
coh_played = zeros(nTrials,1);
numTones_played = zeros(nTrials,1);
        
%get the number of tones played
freq = stim(:, 4);
rt = task(:, 10);

for j = 1:nTrials
    
    %convert reaction time to number of samples
    samples = floor(rt(j) / 1000 * fs);
    
    %get frequencies corresponding to played samples
    curFreq = freq(j);
    playedFreq = curFreq(delaysamp:samples+delaysamp, :);
    numSamples = sum(playedFreq ~= 0);
    numTones = floor(numSamples / 1220);
    numHi = floor(sum(playedFreq == hiFreq)/1220);
    numLo = floor(sum(playedFreq == loFreq)/1220);

    %determine more popular pitch of played tones
    numTones_played(j) = numTones;
    coh_played(j) = numHi/numTones;
    isH_played(j) = numHi > numLo;
end

%% Create one final data file and data file for DDM
data_folder = '/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/Data/';
save_filename = ['All_Data_' meta_data.subject '_' meta_data.date];

visualModes = task(:, 2);
cohLevels = task(:, 3);
speaker = task(:, 4);
waveforms = stim(:, 2);
maskers = stim(:, 3);
isH = stim(:, 5);
choices = task(:, 5);
success = task(:, 6);
stimStart = task(:, 7);
stimStop = task(:, 8);
responseTimeStamp = task(:,9);

%matlab data table that includes waveform
all_data = table((1:nTrials)', visualModes, cohLevels, coh_played,....
    numTones_played, speaker, waveforms, maskers, isH, isH_played, choices, success, stimStart,...
    stimStop, responseTimeStamp, rt, 'VariableNames', {'trialID','visualMode',...
    'cohLevel','coh_played','numTones_played','speaker', 'waveform','isH',...
    'isH_played','choice','success','stimStartTime','stimStopTime',...
    'responseTimeStamps','RT'});
%csv file for DDM code - only includes specific columns
ddm_data = table(coh_played, choices - 1, rt, success);

%save matlab data table
save([data_folder save_filename '.mat'], 'all_data', 'meta_data');
save([data_folder 'DDM_' save_filename '.mat'], 'ddm_data', 'meta_data');

end

