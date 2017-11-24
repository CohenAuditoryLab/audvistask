%% prepare to run AudVisTaskUnlinked_Practice

%clean all existing data
clear;
close all;
topsDataLog.flushAllData();

%dispInd = 0 for small screen, 1 for full screen, >1 for external monitors
[task, list] = AudVisTaskUnlinked_Practice(2); 
task.run();

%clear again
clear
close all;
