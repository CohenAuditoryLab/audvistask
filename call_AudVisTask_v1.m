%% prepare to run AudVisTask_v1

clear;
close all;
topsDataLog.flushAllData();

[task, list] = AudVisTask_v1(0); 
task.run();

%Post-Processing 
data_folder = '/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/';
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

data_table_stim = table((1:nTrials)', visualModes,cohLevels,coh_played,....
    numTones_played,waveforms,isH,isH_played,choices,success,stimStart,...
    stimStop,responseTimeStamp,rt,'VariableNames',{'trialID','visualMode',...
    'cohLevel','coh_played','numTones_played','waveform','isH',...
    'isH_played','choice','success','stimStartTime','stimStopTime',...
    'responseTimeStamps','RT'});
data_table_nostim = table((1:nTrials)',visualModes,cohLevels,coh_played,...
    numTones_played,isH,isH_played,choices,success,stimStart,...
    stimStop,responseTimeStamp,rt,'VariableNames',{'trialID','visualMode',...
    'cohLevel','coh_played','numTones_played','isHigh',...
    'playedHigh','choice','success','stimStartTime','stimStopTime',...
    'responseTimeStamps','RT'});

save([data_folder save_filename '_table.mat'], 'data_table_stim', 'meta_data');
writetable(data_table_nostim, strcat(save_filename,'.csv'));

clear
close all;
