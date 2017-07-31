function AudVis_PostProcess
%% File selection

% Select files
[task_file, tpath] = uigetfile('*.mat', 'Task performance data');
[stim_file, spath] = uigetfile('*.mat', 'Stimulus data');

% Load files
load([tpath task_file]);
load([spath stim_file]);
task = task_table;
stim = stim_table;

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

delay = meta_data.delay; %s
delaysamp = floor(delay * fs);

isH_played = zeros(nTrials,1);
coh_played = zeros(nTrials,1);
numTones_played = zeros(nTrials,1);
        
%get the number of tones played
freq = stim.tone_frequencies;
rt = task.RT;

for j = 1:nTrials
    
    %convert reaction time to number of samples
    samples = floor(rt(j) / 1000 * fs);
    if ~isnan(samples)
        %get frequencies corresponding to played samples
        curFreq = freq{j};
        playedFreq = curFreq(1:samples);
        numSamples = sum(playedFreq ~= 0);
        numTones = floor(numSamples / 1220);
        numHi = floor(sum(playedFreq == hiFreq)/1220);
        numLo = floor(sum(playedFreq == loFreq)/1220);

        %determine more popular pitch of played tones
        numTones_played(j) = numTones;
        coh_played(j) = numHi/numTones;
        isH_played(j) = numHi > numLo;
    else 
        numTones_played(j) = NaN;
        coh_played(j) = NaN;
        isH_played(j) = NaN;
    end
end

%% Create one final data file and data file for DDM
data_folder = '/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/Data/';
save_filename = ['All_Data_' meta_data.subject '_' meta_data.date];

visualModes = task.visualMode;
cohLevels = task.cohLevel;
speaker = task.speaker;
waveforms = stim.waveform;
maskers = stim.masker;
isH = stim.isHigh;
choices = task.choice;
success = task.success;
stimStart = task.stimStartTime;
stimStop = task.stimStopTime;
responseTimeStamp = task.responseTimeStamps;

%matlab data table that includes waveform
all_data = table((1:nTrials)', visualModes, cohLevels, coh_played,....
    numTones_played, speaker, waveforms, maskers, isH, isH_played, choices, success, stimStart,...
    stimStop, responseTimeStamp, rt, 'VariableNames', {'trialID','visualMode',...
    'cohLevel','coh_played','numTones_played','speaker', 'waveform','masker', 'isH',...
    'isH_played','choice','success','stimStartTime','stimStopTime',...
    'responseTimeStamps','RT'});
%csv file for DDM code - only includes specific columns
ddm_data = [coh_played, choices - 1, rt, success];

%save matlab data table
save([data_folder save_filename '.mat'], 'all_data', 'meta_data');
save([data_folder 'DDM_' save_filename '.mat'], 'ddm_data');

end

