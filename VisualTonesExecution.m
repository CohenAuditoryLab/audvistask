%assuming each trial has length 2 seconds with about 4 seconds in
%between for a total of 6 seconds
%in 20 minutes (1200 seconds), there can be 200 trials
%this works out to 50 trials per each of the 4 block options

%generate the order of the visual stimuli for participant
order = randperm(4);

%add the corresponding modes to a cell array in the random order 
counter = 1;
visual = cell(1,4);
for j = order
    if j == 1
        type = 'None';
    elseif j == 2
        type = 'Low';
    elseif j == 3
        type = 'High';
    else
        type = 'All';
    end
    visual{counter} = type;
    counter = counter + 1;
end

%for each of the 4 modes
for h = 1 : 4
    % for 50 trials per mode
    for i = 1:50
        
        %select an auditory coherence from Tsunada paper
        %coherences used were 0%, 50%, 100%
        %%%%% NO OTHER VALUES WERE SPECIFIED IN THE METHODS SECTION OF THE
        %%%%% PAPER
        %%%%% I SUBSTITUTED 25 AND 75 FOR VARIATION BUT THESE WILL LIKELY NEED
        %%%%% TO BE ADJUSTED
        
        %randomly generate a coherence 
        coherences = [0, .25, .5, .75, 1];
        index = randsample(5,1);
        c = coherences(index);
        
        %play the tone for each coherence 
        [~, ~] = VisualTones(8000, 24000, c, visual{h});
        pause(4);
        
    end
end