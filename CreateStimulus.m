function [target_auditory, target_stimulus, frequencies, isHigh,...
    masker_auditory, masker_stimulus] = CreateStimulus(loFreq, hiFreq, c, mode)
%Create the target stimulus (loFreq and hiFreq set, c and mode random) and
%the masker stimulus (loFreq and hiFreq same, c = 50% and mode always
%'None') before adding speaker cue.

[target_auditory, target_stimulus, frequencies, isHigh] = ...
    VisualTones2(loFreq, hiFreq, c, mode);

[masker_auditory, masker_stimulus, ~, ~] = VisualTones2(loFreq, hiFreq, ...
    0.50, 'None');

end