clear all;

% CREATE AUDIOPLAYERS
player = dotsPlayableWave_TDT();
player.sampleFrequency = 24414;
player.intensity = 1;
% player2 = dotsPlayableWave_TDT2();
% player2.sampleFrequency = 24414;
% player2.intensity = 1;

[waveform, full_target, f, h, masker_waveform, full_masker] = ...
    CreateStimulus(100, 200, .5, 'High');
%add speaker cue to stimuli
[new_target, new_masker] = AddSpeakerCue(full_target, full_masker);

%decide which speaker gets target and which gets masker
% r = 2 * rand;
% if r <= 1
%     chosen = 'Speaker 1';
    player.wave = new_target;
    player.masker = new_masker;
% else
%     chosen = 'Speaker 2';
%     player.wave = new_masker;
%     player2.wave = new_target;
% end

%play stimuli in correct speakers  
player.prepareToPlay;
%player2.prepareToPlay;
drawnow;
player.play;
%player2.play;

tic;
pause(5)
disp(toc);

player.stop;
%player2.stop;

% ----------------------------------------------------- %
% %RCXFile = fullfile('C:\','work', 'AudVis_Task','Speaker_Test.rcx');
% RCXFile = fullfile('C:\','work', 'AudVis_Task','AudVis_Circuit.rcx');
% 
% RP = actxcontrol('RPco.x',[5 5 26 26]);
% 
% if RP.ConnectRX6('GB', 1)
%     disp('Connected to RX6!');
% else
%     disp('Unable to connect to RX6');
% end
% 
% % load rcx file
% RP.LoadCOF(RCXFile);
% RP.Run;
% 
% RP.WriteTagVEX('lightin1', 0, 'F32', full(2,:));
% RP.WriteTagVEX('datain1', 0, 'F32', data);
% drawnow;
% RP.WriteTagVEX('lightin2', 0, 'F32', full2(2,:));
% RP.WriteTagVEX('datain2', 0, 'F32', data2);
% drawnow;
% 
% RP.SoftTrg(3); %Ch2
% RP.SoftTrg(1); %Ch1
% 
% tic;
% pause(4)
% disp(toc);
% 
% RP.SoftTrg(4); %Ch2
% RP.SoftTrg(2); %Ch1
% RP.Halt;
% RP.ClearCOF();
