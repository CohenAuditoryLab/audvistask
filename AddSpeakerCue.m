function [new_target, new_masker] = AddSpeakerCue(target_stimulus, masker_stimulus)
%Adds a series of LED bursts before the target_stimulus and then a 500ms
%delay so that the subject knows which speaker to attend to. Adds an
%equivalent delay (zeros) before the masker_stimulus so that the two
%stimuli are synchronous.

%number of samples for one 50ms visual burst
duration = 50; %ms
samplingFreq = 44100; %Hz
samples = duration / 1000 * samplingFreq;
vburst = ones(1, samples);
aburst = zeros(1, samples);

%create complete visual cue
vcue = [aburst, vburst, aburst, vburst, aburst, vburst];
acue = repmat(aburst, 1, 6);

%create 500ms delay
delay = 500; %ms
sampDelay = delay / 1000 * samplingFreq;
delayvec = zeros(1, sampDelay);

%full cue 
vcue = [vcue delayvec];
acue = [acue delayvec];

%append onto beginning of stimuli
targ = cat(1, acue, vcue);
mask = cat(1, acue, acue);

new_target = [targ, target_stimulus];
new_masker = [mask, masker_stimulus];

% time = linspace(0, 4000, 2205*40);
% nt = new_target(1, :);
% na = new_target(2, :);
% plot(time, 100.*nt(1:2205*40));
% hold on;
% plot(time, na(1:2205*40), 'LineWidth', 3);
% 
% figure();
% nm = new_masker(1,:);
% ma = new_masker(2,:);
% plot(time, 100.*nm(1:2205*40));
% hold on;
% plot(time, ma(1:2205*40), 'LineWidth', 3);

end

