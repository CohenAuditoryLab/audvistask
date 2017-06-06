%% prepare to run AudVisTask_v1

%clean all existing data
clear;
close all;
topsDataLog.flushAllData();

%dispInd = 0 for small screen, 1 for full screen, >1 for external monitors
[task, list] = AudVisTask_v1(1); 
task.run();

%Post-Processing 
data_folder = '/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/Data/';
save_filename = list{'meta'}{'saveFilename'};

% create data table 
meta_data.subject = list{'meta'}{'subjID'};
meta_data.date = list{'meta'}{'date'};

hd = list{'Stimulus'}{'header'};
meta_data.loFreq = hd.loFreq;
meta_data.hiFreq = hd.hiFreq;
meta_data.toneDur = hd.toneDur;
meta_data.trialDur = hd.trialDur;
meta_data.fs = hd.fs;
nTrials = list{'Counter'}{'trial'};
meta_data.nTrials = nTrials;

% trial data 
visualModes = list{'control'}{'visualModes'};
cohLevels = list{'control'}{'cohLevels'};
coh_played = list{'Stimulus'}{'coh_played'};
isH = list{'Stimulus'}{'isH'};
isH_played = list{'Stimulus'}{'isH_played'};
numTones_played = list{'Stimulus'}{'numTones_played'};
waveforms = list{'Stimulus'}{'waveforms'};
stimStart = list{'Timestamps'}{'stim_start'};
stimStop = list{'Timestamps'}{'stim_stop'};
responseTimeStamp = list{'Timestamps'}{'choices'};
choices = list{'Input'}{'choices'};
success = list{'Input'}{'corrects'};
rt = list{'Input'}{'RT'};

%matlab data table that includes waveform
data_table_stim = table((1:nTrials)', visualModes,cohLevels,coh_played,....
    numTones_played,waveforms,isH,isH_played,choices,success,stimStart,...
    stimStop,responseTimeStamp,rt,'VariableNames',{'trialID','visualMode',...
    'cohLevel','coh_played','numTones_played','waveform','isH',...
    'isH_played','choice','success','stimStartTime','stimStopTime',...
    'responseTimeStamps','RT'});
%csv file that does not include waveform
data_table_nostim = table((1:nTrials)',visualModes,cohLevels,coh_played,...
    numTones_played,isH,isH_played,choices,success,stimStart,...
    stimStop,responseTimeStamp,rt,'VariableNames',{'trialID','visualMode',...
    'cohLevel','coh_played','numTones_played','isHigh',...
    'playedHigh','choice','success','stimStartTime','stimStopTime',...
    'responseTimeStamps','RT'});
%csv file for DDM code - only includes specified columns
%data_table_ddm = table(coh_played, choices - 1, rt, success);

cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/Data/');
%save matlab data table (to keep track of waveform)
save([data_folder save_filename '_table.mat'], 'data_table_stim', 'meta_data');
%save csv file with all data except waveform
writetable(data_table_nostim, strcat(save_filename,'.csv'));
%save csv file for DDM (column 1 - coherence from 0 to 1, column 2 - choice
%minus 1, column 3 - RT in ms) 
%writetable(data_table_ddm, strcat('DDM_', save_filename, '.csv'));

%clear again
clear
close all;
