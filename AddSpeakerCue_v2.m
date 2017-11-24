function [new_target, new_masker] = AddSpeakerCue_v2(target_stimulus, masker_stimulus)
%Adds a series of tone bursts before the target_stimulus and then a 1000ms
%delay so that the subject knows which speaker to attend to. Adds an
%equivalent delay (zeros) before the masker_stimulus so that the two
%stimuli are synchronous.

samplingFreq = 24414;
%create sounds
s1 = PreStimulus_Unlinked(440);
s2 = PreStimulus_Unlinked(880);
s3 = PreStimulus_Unlinked(1760);
blank = zeros(1, 3*length(s1));

%create complete auditory cue
tcue = [s1, blank, s2, blank, s3];
mcue = zeros(1, length(tcue));

%create 1000ms delay
delay = 1000; %ms
sampDelay = floor(delay / 1000 * samplingFreq);
delayvec = zeros(1, sampDelay);

%full cue 
tcue = [tcue delayvec];
mcue = [mcue delayvec];
acue = zeros(1, length(tcue));

tcue = cat(1, tcue, acue);
mcue = cat(1, mcue, acue);

new_target = [tcue, target_stimulus];
new_masker = [mcue, masker_stimulus];

sound(new_target, samplingFreq);

figure();
time = linspace(0, 4000, 2205*40);
nt = new_target(1, :);
na = new_target(2, :);
plot(time, nt(1:2205*40));
hold on;
plot(time, na(1:2205*40), 'LineWidth', 3, 'Color', 'r');

figure();
nm = new_masker(1,:);
ma = new_masker(2,:);
plot(time, nm(1:2205*40));
hold on;
plot(time, ma(1:2205*40), 'LineWidth', 3, 'Color', 'r');

end

