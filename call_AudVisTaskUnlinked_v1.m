%% prepare to run AudVisTaskUnlinked_v1

%clean all existing data
clear;
close all;
topsDataLog.flushAllData();

%dispInd = 0 for small screen, 1 for full screen, >1 for external monitors
[task, list] = AudVisTaskUnlinked_v1(2); 
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
meta_data.delay = hd.delay;
nTrials = list{'Counter'}{'trial'};
meta_data.nTrials = nTrials;

% trial data 
visualModes = list{'control'}{'visualModes'};
cohLevels = list{'control'}{'cohLevels'};
speaker = list{'Stimulus'}{'speaker'};
led = list{'Stimulus'}{'LED'};
stimStart = list{'Timestamps'}{'stim_start'};
stimStop = list{'Timestamps'}{'stim_stop'};
responseTimeStamp = list{'Timestamps'}{'choices'};
choices = list{'Input'}{'choices'};
success = list{'Input'}{'corrects'};
rt = list{'Input'}{'RT'};

%matlab data table that includes waveform
task_table = table((1:nTrials)', visualModes, cohLevels, speaker, led, ....
   choices,success,stimStart,stimStop,responseTimeStamp,rt,'VariableNames',...
   {'trialID','visualMode','cohLevel','speaker', 'LED', 'choice','success',...
   'stimStartTime','stimStopTime','responseTimeStamps','RT'});

cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/Data/');
%save matlab data table (to keep track of waveform)
save([data_folder save_filename '_table.mat'], 'task_table', 'meta_data');

%clear again
clear
close all;
