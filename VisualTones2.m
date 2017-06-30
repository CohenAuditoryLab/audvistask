function [auditory, stimulus, frequencies, isHigh] = VisualTones2(loFreq, hiFreq, c, mode)

%VISUALTONES uses frequency boundaries of the form [loFreq, hiFreq]
%and a coherence factor (c) to produce a sound stimulus of 2 seconds.  The
%coherence is a value between 0 and 1, where 0 means all frequencies are low
%and 1 means all frequencies are high; 0.5 represents 50% low sound bursts
%and 50% high sound bursts.  Mode is a string input that takes the form
%"None", "Low", "High", or "All" and corresponds to the sounds during which
%the LED should be on (1) versus off (0). No visual stimulus will be
%produced if another string is entered for this value.
%% Constants and time vector

%duration in ms
duration = 50;
%sampling frequency in samples/sec
samplingFreq = 44100; %Hz
%number of samples needed
samples = duration / 1000 * samplingFreq;
%time vector
t = linspace(0, duration, samples);

%% Create 2 possible auditory stimuli

%sine vector for one sound burst
%create HIGH FREQUENCY sine wave
xHi = sin((2 * pi * (hiFreq/1000)) .* t);
%create LOW FREQUENCY sine wave
xLo = sin((2 * pi * (loFreq/1000)) .* t);

%gate the sine wave
%since 5ms up, 50ms total -> need 1/10 of samples
gatedSineHi = zeros(1, samples);
gatedSineLo = zeros(1, samples);
%add 0.5 to make whole number of samples
numRise = (5 / 1000 * samplingFreq);
numRiseR = numRise + 0.5;
%ramp up
gatedSineHi (1 : numRiseR)= xHi(1 : numRiseR) ...
    .* (0.2 .* t(1 : numRiseR));
gatedSineLo (1 : numRiseR) = xLo(1 : numRiseR) ...
    .* (0.2 .* t(1 : numRiseR));
%flat
gatedSineHi(numRiseR + 1 : samples - numRiseR) = ...
    xHi(numRiseR + 1 : samples - numRiseR);
gatedSineLo(numRiseR + 1 : samples - numRiseR) = ...
    xLo(numRiseR + 1 : samples - numRiseR);
%ramp down
gatedSineHi(samples - numRiseR + 1: end) = ...
    xHi(samples - numRiseR + 1 : end) .* ...
    (-0.2 .* t(samples - numRiseR + 1 : end) + 10);
gatedSineLo(samples - numRiseR + 1: end) = ...
    xLo(samples - numRiseR + 1 : end) .* ...
    (-0.2 .* t(samples - numRiseR + 1 : end) + 10);

%% Initialize PPP

%created the Poisson Point Process
rate = 16.7; %Poisson Process rate (lambda) in Hz
dur = 4; %in seconds
ppp = rand(1, dur * samplingFreq) < (rate / samplingFreq);

%initialize auditory and visual vectors
aud = zeros(1, length(ppp));
vis = zeros(1, length(ppp));
%list of frequencies
frequencies = zeros(length(ppp),1);

%% Conditions to fill in the auditory and visual stimulus vectors

%find the ones in the ppp vector
i = find(ppp == 1);

%for each of these indices
for jj = 1:length(i)
    %decide based on coherence whether the stimulus should be high or low
    n = rand;
    if n <= c
        auditory = gatedSineHi;
        freq = hiFreq;
    else
        auditory = gatedSineLo;
        freq = loFreq;
    end
    
    %% Create 4 possible visual stimuli
    %create visual stimulus
    %if input string is None
    if strcmp(mode, 'None')
        visual = zeros(1, length(gatedSineHi));
    end
    
    %if input string is Low
    if strcmp(mode, 'Low')
        if freq == loFreq
            visual = ones(1, length(gatedSineLo));
        else
            visual = zeros(1, length(gatedSineLo));
        end
    end
    
    %if input string is High
    if strcmp(mode, 'High')
        if freq == hiFreq
            visual = ones(1, length(gatedSineHi));
        else
            visual = zeros(1, length(gatedSineHi));
        end
    end
    
    %if input string is All
    if strcmp(mode, 'All')
        visual = ones(1, length(gatedSineHi));
    end
    %% Create visual and auditory vectors 
    
    ind_end = min(length(ppp), i(jj) + length(auditory) -1);
    if length(aud(i(jj) : ind_end)) < length(auditory)
        auditory = auditory(1:length(aud(i(jj) : ind_end)));
        visual = visual(1:length(vis(i(jj) : ind_end)));
    end
    auditory = auditory + aud(i(jj) : ind_end);
    aud(i(jj) : ind_end) = 0.015 .* auditory;
    vis(i(jj) : ind_end) = visual;
    
    frequencies(i(jj) : ind_end, :) = freq;
end
%% Outputs

%concatenate
stimulus = cat(1, aud, vis);

numLo = sum(frequencies == loFreq);
numHi = sum(frequencies == hiFreq);

isHigh = numHi > numLo;

%% Graphing and playing for testing purposes
% time = linspace(0, 4000, 2205*40);
% plot(time, 100.*aud(1:2205*40));
% hold on;
% plot(time, vis(1:2205*40), 'LineWidth', 3);
% ylim([-3.5 3.5])

%sound(stimulus, 44100);
%disp(stimulus);
end

