function [task, list] = AudVisTask_v2(dispInd)
%% 05-22-2017 created by Brianna - Auditory Visual Task
% all 5 blocks with a break after block 3
%% Setting up the screen

isClient = false;

sc = dotsTheScreen.theObject;
%dispInd = 0 for small screen, 1 for full screen, >1 for external monitors
sc.reset('displayIndex', dispInd);

%Call GetSecs to load up the Mex files for getting time, so no delays later
GetSecs;

%Subject ID and date/time info
subj_id = input('Subject ID: ','s');
cur_date = datestr(now,'yymmdd');
cur_time = datestr(now,'HHMM');
cur_task = mfilename;
save_filename = [cur_task '_' subj_id '_' cur_date '_' cur_time];
%% Setting up a list structure

%create list structure
list = topsGroupedList(cur_task);

%Save subject/date/time info as meta data 
list{'meta'}{'subjID'} = subj_id;
list{'meta'}{'date'} = cur_date;
list{'meta'}{'time'} = cur_time;
list{'meta'}{'task'} = cur_task;
list{'meta'}{'saveFilename'} = save_filename;

%% Settings for generating the sequence of conditions

% number visual modes 
block_size = 5;
% number of trials per visual mode
block_rep = 80; %1 %15 %50 %75
% possible visual values to select from
vis_vals = {'None', 'Low', 'High', 'All', 'Random'};

%% Visual conditions for each trial

taskConditions = topsConditions(cur_task);

vis_parameter = 'visualMode';
taskConditions.addParameter(vis_parameter, vis_vals);
likesVisMode = topsFoundation();
taskConditions.addAssignment('visualMode', likesVisMode, '.', 'name');

nTrials = block_rep * block_size;
visualModes = cell(nTrials, 1);
%select each visual mode once in a random order
%each mode will only be selected once per trial
taskConditions.setPickingMethod('shuffledEach',1);

keepGoing = true;
counter = 0;
while keepGoing
    taskConditions.run();
    %add the visual modes in for the length of the block
    for k = counter + 1 : counter + block_rep
         visualModes{k} = likesVisMode.name;
    end 
    %continue until task conditions are finished being selected
    keepGoing = ~taskConditions.isDone;
    %increment by the size of the block
    counter = counter + block_rep;
end 

list{'control'}{'task conditions'} = taskConditions;
list{'control'}{'visualModes'} = visualModes;

%% Generate coherence for each trial randomly

cohLevels = zeros(nTrials, 1);

%Create COUNTER object
list{'Counter'}{'trial'} = 0;

%possible coherences
%most data should be collected between about 25 and 75
coherences = [0, .25, .33, .40, .50, .60, .67, .75, 1];

for i = 1:nTrials
    %add coherence as a condition
    tmp_coh = topsConditions('coh');
    tmp_coh.addParameter('cohLevel', coherences);
    likesCohLevel = topsFoundation();
    tmp_coh.addAssignment('cohLevel', likesCohLevel, '.', 'name');

    %randomly generate a coherence
    index = randsample(9,1);
    c = coherences(index);
    cohLevels(i) = c;
end 

list{'control'}{'cohLevels'} = cohLevels;
%% Audio Settings

hd.loFreq = 5000; %hz      312.5 |  625 | 1250 | 2500 |  5000
hd.hiFreq = 20000; %hz     625   | 1250 | 2500 | 5000 | 10000
hd.toneDur = 50; %ms
hd.toneSOA = 10; %ms, actually poisson random number centered around 10 
hd.trialDur = 4000; %ms
hd.fs = 44100; %samples/sec

% INPUT PARAMETERS
responsewindow = hd.trialDur/2; %time allowed to respond = half trial duration, ms
list{'Input'}{'responseWindow'} = responsewindow;

% CREATE AUDIOPLAYER
player = dotsPlayableWave_2Channel();
player.sampleFrequency = hd.fs;
player.duration = hd.trialDur; %ms
player.intensity = 1;
%% Time Variables
iti = 2; %seconds
list{'timing'}{'intertrial'} = iti; %intertrial interval, between tones
%% Input Settings

% Set up gamepad object
gp = dotsReadableHIDGamepad();

if gp.isAvailable
    
    %use gamepad if connected
    ui = gp;
    
    %define movements, must be held down
    %map x-axis -1 to left and +1 to right
    isLeft = [gp.components.ID] == 9;
    isA = [gp.components.ID] == 3;
    isRight = [gp.components.ID] == 10;
    
    Left = gp.components(isLeft);
    A = gp.components(isA);
    Right = gp.components(isRight);
    
    gp.setComponentCalibration(Left.ID, [], [], [0 +2]);
    gp.setComponentCalibration(A.ID, [], [], [0 +3]);
    gp.setComponentCalibration(Right.ID, [], [], [0 +4]);
    
    %undefine any default events
    IDs = gp.getComponentIDs();
    for k = 1:numel(IDs)
        gp.undefineEvent(IDs(k));
    end
    
    %define values for relevant button presses
    gp.defineEvent(Left.ID, 'left', 0, 0, true);
    gp.defineEvent(A.ID, 'continue', 0, 0, true);
    gp.defineEvent(Right.ID, 'right', 0, 0, true);
    
else
    
    %if gamepad not available, use keyboard
    kb = dotsReadableHIDKeyboard();
    
    %define movements, must be held down
    %left = +2, up = +3, right = +4
    isLeft = strcmp({kb.components.name}, 'KeyboardF');
    isSpace = strcmp({kb.components.name}, 'KeyboardSpacebar');
    isRight = strcmp({kb.components.name}, 'KeyboardJ');
    
    Left = kb.components(isLeft);
    Space = kb.components(isSpace);
    Right = kb.components(isRight);
    
    kb.setComponentCalibration(Left.ID, [], [], [0 +2]);
    kb.setComponentCalibration(Space.ID, [], [], [0 +3]);
    kb.setComponentCalibration(Right.ID, [], [], [0 +4]);
    
    %undefine default keyboard events
    IDs = kb.getComponentIDs();
    for j = 1:numel(IDs)
        kb.undefineEvent(IDs(j));
    end
    
    %define keyboard events
    %fire once event if held down
    kb.defineEvent(Left.ID, 'left', 0, 0, true);
    kb.defineEvent(Space.ID, 'continue', 0, 0, true);
    kb.defineEvent(Right.ID, 'right', 0, 0, true);
    
    ui = kb;
end

%Make sure the UI is running on the same clock as everything else
%Use operating system time as absolute clock
ui.clockFunction = @GetSecs;

%Store UI in list bin to access from functions
ui.isAutoRead = 1;
list{'Input'}{'controller'} = ui;
%% Store data in the list structure

%STIMULUS INFORMATION
list{'Stimulus'}{'header'} = hd;
list{'Stimulus'}{'player'} = player;
list{'Stimulus'}{'waveforms'} = cell(nTrials,1);
list{'Stimulus'}{'freq'} = cell(nTrials,1);
list{'Stimulus'}{'isH'} = zeros(nTrials,1);
list{'Stimulus'}{'bursts'} = zeros(nTrials,136);
list{'Stimulus'}{'isH_played'} = zeros(nTrials,1);
list{'Stimulus'}{'coh_played'} = zeros(nTrials,1);
list{'Stimulus'}{'numTones_played'} = zeros(nTrials,1);

%TIMESTAMPS
list{'Timestamps'}{'stim_start'} = zeros(nTrials,1);
list{'Timestamps'}{'stim_stop'} = zeros(nTrials,1);
list{'Timestamps'}{'choices'} = zeros(nTrials,1);

%INPUT INFORMATION
list{'Input'}{'choices'} = zeros(nTrials,1);
list{'Input'}{'corrects'} = zeros(nTrials,1);
list{'Input'}{'RT'} = zeros(nTrials,1);
%% Graphics

%Define colors for targets
list{'Graphics'}{'gray'} = [0.5 0.5 0.5];
list{'Graphics'}{'red'} = [0.75 0.25 0.1];
list{'Graphics'}{'green'} = [.25 0.75 0.1];

%Text prompts
%Low pitch label (right side)
lowlabel = dotsDrawableText();
lowlabel.string = 'Low';
lowlabel.fontSize = 36;
lowlabel.typefaceName = 'Calibri';
lowlabel.isVisible = false;
lowlabel.x = 5;
lowlabel.y = 3;

%High pitch label (left side)
highlabel = dotsDrawableText();
highlabel.string = 'High';
highlabel.fontSize = 36;
highlabel.typefaceName = 'Calibri';
highlabel.isVisible = false;
highlabel.x = -5;
highlabel.y = 3;

%Block Label (top)
blocklabel = dotsDrawableText();
blocklabel.string = sprintf('Block Number 1 of 5');
blocklabel.typefaceName = 'Calibri';
blocklabel.isVisible = false;
blocklabel.x = 0;
blocklabel.y = 5.5;

readyprompt = dotsDrawableText();
readyprompt.string = 'Ready?';
readyprompt.fontSize = 42;
readyprompt.typefaceName = 'Calibri';
readyprompt.isVisible = true;

buttonprompt = dotsDrawableText();
buttonprompt.string = 'press A to get started';
buttonprompt.fontSize = 24;
buttonprompt.typefaceName = 'Calibri';
buttonprompt.y = -2;
buttonprompt.isVisible = true;

readyprompt2 = dotsDrawableText();
readyprompt2.string = 'Congratulations! Your performance is ';
readyprompt2.fontSize = 30;
readyprompt2.typefaceName = 'Calibri';
readyprompt2.isVisible = false;

buttonprompt2 = dotsDrawableText();
buttonprompt2.string = 'press A to quit';
buttonprompt2.fontSize = 24;
buttonprompt2.typefaceName = 'Calibri';
buttonprompt2.y = -2;
buttonprompt2.isVisible = false;

%Create a cursor dot to indicate user selection/provide feedback
cursor = dotsDrawableTargets();
cursor.colors = list{'Graphics'}{'gray'};
cursor.width = 1.5;
cursor.height = 1.5;
cursor.xCenter = 0;
cursor.yCenter = 0;
cursor.isVisible = false;
list{'Graphics'}{'cursor'} = cursor;

%Graphical ensemble
ensemble = dotsEnsembleUtilities.makeEnsemble('drawables', isClient);
target = ensemble.addObject(cursor);
ready = ensemble.addObject(readyprompt);
button = ensemble.addObject(buttonprompt);
ready2 = ensemble.addObject(readyprompt2);
button2 = ensemble.addObject(buttonprompt2);
lowlabel = ensemble.addObject(lowlabel);
highlabel = ensemble.addObject(highlabel);
blocklabel = ensemble.addObject(blocklabel);

list{'Graphics'}{'ensemble'} = ensemble;
list{'Graphics'}{'target'} = target;
list{'Graphics'}{'ready2'} = ready2;
list{'Graphics'}{'button2'} = button2;
list{'Graphics'}{'low'} = lowlabel;
list{'Graphics'}{'high'} = highlabel;
list{'Graphics'}{'block'} = blocklabel;

% tell the ensembles how to draw a frame of graphics
% the static drawFrame() takes a cell array of objects
ensemble.automateObjectMethod(...
    'draw', @dotsDrawable.drawFrame, {}, [], true);

% also put dotsTheScreen into its own ensemble
screen = dotsEnsembleUtilities.makeEnsemble('screen', isClient);
screen.addObject(dotsTheScreen.theObject());
list{'Graphics'}{'screen'} = screen;

% automate the task of flipping screen buffers
screen.automateObjectMethod('flip', @nextFrame);
%% Control

% a batch of function calls that apply to all the trial types below
% start- and finishFevalable get called once per trial
% addCall() accepts fevalables to be called repeatedly during a trial

trialCalls = topsCallList();
trialCalls.addCall({@read, ui}, 'read input');
list{'control'}{'trial calls'} = trialCalls;
%% State Machine

show = @(index) ensemble.setObjectProperty('isVisible', true, index); %show asset
hide = @(index) ensemble.setObjectProperty('isVisible', false, index); %hide asset

%Prepare Machine - used in antetask
prepareMachine = topsStateMachine();
prepStates = {'name', 'entry', 'input', 'exit', 'timeout', 'next';
    'Ready', {},      {},      {@waitForCheckKey list},     0,       'Hide';
    'Hide', {hide [ready button]}, {}, {}, 0, 'Finish'
    'Finish', {}, {}, {}, 0, '';};
prepareMachine.addMultipleStates(prepStates);

list{'control'}{'prepareMachine'} = prepareMachine;

% State Machine - used in maintask
mainMachine = topsStateMachine();
mainStates = {'name', 'entry', 'input', 'exit', 'timeout', 'next';
    'CheckReady', {@startTrial list block_rep}, {}, {@waitForCheckKey list}, 0, 'Stimulus';
    'Stimulus', {@playstim list}, {}, {@waitForChoiceKey list}, 0, 'Feedback';
    'Feedback', {@showFeedback list}, {}, {}, 0, 'Exit';
    'Exit',{@finishTrial list}, {}, {}, iti,''};
mainMachine.addMultipleStates(mainStates);

list{'control'}{'mainMachine'} = mainMachine;

% End Machine - used in post-task
endMachine = topsStateMachine();
endStates = {'name', 'entry', 'input', 'exit', 'timeout', 'next';
    'Ready', {@startEndTask list},      {},      {@waitForCheckKey list},     0,       'Hide';
    'Hide', {hide [ready2 button2]}, {}, {}, 0, 'Finish';
    'Finish', {}, {}, {}, 0, '';};
endMachine.addMultipleStates(endStates);

list{'control'}{'endMachine'} = endMachine;

prepareConcurrents = topsConcurrentComposite();
prepareConcurrents.addChild(ensemble);
prepareConcurrents.addChild(prepareMachine);
prepareConcurrents.addChild(screen);

% add a branch to the tree trunk to lauch a Fixed Time trial
prepareTree = topsTreeNode();
prepareTree.addChild(prepareConcurrents);

mainConcurrents = topsConcurrentComposite();
mainConcurrents.addChild(ensemble);
mainConcurrents.addChild(trialCalls);
mainConcurrents.addChild(mainMachine);
mainConcurrents.addChild(screen);

mainTree = topsTreeNode();
mainTree.iterations = nTrials;
mainTree.addChild(mainConcurrents);

endConcurrents = topsConcurrentComposite();
endConcurrents.addChild(ensemble);
endConcurrents.addChild(endMachine);
endConcurrents.addChild(screen);

% add a branch to the tree trunk to lauch a Fixed Time trial
endTree = topsTreeNode();
endTree.addChild(endConcurrents);

% Top Level Runnables
task = topsTreeNode();
task.startFevalable = {@callObjectMethod, screen, @open};
task.finishFevalable = {@callObjectMethod, screen, @close};
task.addChild(prepareTree);
task.addChild(mainTree);
task.addChild(endTree);
end

%% Accessory Functions

function startEndTask(list)
    ensemble = list{'Graphics'}{'ensemble'};

    corrects = list{'Input'}{'corrects'};
    %calculate performance
    perf = 100*sum(corrects)/length(corrects);

    %prepare text + performance
    ready2 = list{'Graphics'}{'ready2'};
    button2 = list{'Graphics'}{'button2'};
    tmp_str = ensemble.getObjectProperty('string', ready2);
    tmp_str = [tmp_str num2str(perf) ' %'];
    ensemble.setObjectProperty('string', tmp_str, ready2);

    %make visible
    ensemble.setObjectProperty('isVisible', true, ready2);
    ensemble.setObjectProperty('isVisible', true, button2);
end

function startTrial(list, block_rep)
    %clear last trial data
    ui = list{'Input'}{'controller'};
    ui.flushData();

    %increment counter to label trial
    counter = list{'Counter'}{'trial'};
    counter = counter + 1;
    list{'Counter'}{'trial'} = counter;

    %set high and low labels to be visible 
    ensemble = list{'Graphics'}{'ensemble'};
    low = list{'Graphics'}{'low'};
    high = list{'Graphics'}{'high'};
    block = list{'Graphics'}{'block'};
    
    ensemble.setObjectProperty('isVisible', false, block);
    
    b = int16(counter/(block_rep));
    block = dotsDrawableText();
    block.string = sprintf('Block %d of 5', b);
    block.typefaceName = 'Calibri';
    block.isVisible = false;
    block.x = 0;
    block.y = 5.5;
    block = ensemble.addObject(block);
    list{'Graphics'}{'block'} = block;

    ensemble.setObjectProperty('isVisible', true, low);
    ensemble.setObjectProperty('isVisible', true, high);
    ensemble.setObjectProperty('isVisible', true, block);
end

function finishTrial(list)
    %draw the target
    ensemble = list{'Graphics'}{'ensemble'};
    target = list{'Graphics'}{'target'};
    ensemble.setObjectProperty('isVisible', false, target);

    %time between trials
    pause(list{'timing'}{'intertrial'});
end

function showFeedback(list)
    %hide the fixation point and cursor
    ensemble = list{'Graphics'}{'ensemble'};
    target = list{'Graphics'}{'target'};
    counter = list{'Counter'}{'trial'};

    %compare stimulus direction to choice direction
    isCorrect = list{'Input'}{'corrects'};

    %indicate correct or incorrect by coloring in the targets
    if isnan(isCorrect(counter))
        ensemble.setObjectProperty('colors', list{'Graphics'}{'gray'}, target);
        isCorrect(counter) = 0;
    elseif isCorrect(counter)
        ensemble.setObjectProperty('colors', list{'Graphics'}{'green'}, target);
    else
        ensemble.setObjectProperty('colors', list{'Graphics'}{'red'}, target);
    end

    list{'Input'}{'corrects'} = isCorrect;
end

function string = waitForChoiceKey(list)
    %Get list items
    counter = list{'Counter'}{'trial'};
    ensemble = list{'Graphics'}{'ensemble'};
    target = list{'Graphics'}{'target'};
    ui = list{'Input'}{'controller'};
    player = list{'Stimulus'}{'player'};
    freq = list{'Stimulus'}{'freq'};
    hd = list{'Stimulus'}{'header'};
    stim_start = list{'Timestamps'}{'stim_start'};
    responsewindow = list{'Input'}{'responseWindow'};
    choices = list{'Input'}{'choices'};

    % whether it's a high-freq trial
    isH = list{'Stimulus'}{'isH'}; 
    bursts = list{'Stimulus'}{'bursts'};
    isH_played = list{'Stimulus'}{'isH_played'};
    coh_played = list{'Stimulus'}{'coh_played'};
    numTones_played = list{'Stimulus'}{'numTones_played'};

    %clear existing data 
    ui.flushData 

    %initialize variable 
    press = '';

    %wait for keypress
    %start timer 
    tic 
    while ~strcmp(press, 'left') && ~strcmp(press, 'right')
        %Break loop if responsewindow time expires and move to next trial
        if toc > responsewindow 
            choice = NaN;
            timestamp = NaN;
            break
        end 

        %Check for button press 
        press = '';
        read(ui);
        [~, ~, eventname, ~] = ui.getHappeningEvent();
        events = cell(length(coh_played), 1);

        if ~isempty(eventname) && length(eventname) == 1
            press = eventname;
            events(counter) = press;
            %stop the stimulus once a response is detected 
            %avoid error by only stopping if left or right is pressed
            if ~strcmp(events{counter}, 'continue')
                player.stop;
            end

            %get the timestamp of the stimulus stop time 
            stim_stop = list{'Timestamps'}{'stim_stop'};
            stim_stop(counter) = player.stopTime;
            list{'Timestamps'}{'stim_stop'} = stim_stop;
        end 
    end

    %Get intertrial interval
    iti = list{'timing'}{'intertrial'};

    %Get timestamp of button press time 
    if ~isempty(press)
        timestamp = ui.history;
        %to ensure timestamp from a pressed key/button
        timestamp = timestamp(timestamp(:, 2) > 1, :);
        timestamp = timestamp(end);

        %calculate reaction time 
        rt = (timestamp - stim_start(counter)) * 1000; %ms
        %record current choice 
        cur_choice = press{1};

        %get the number of tones played
        b = bursts(counter, :);
        cumT = 0;
        count = 1;
        n = 0;
        while cumT < rt && count <= sum((b ~= 0))
            n = n + 1;
            cumT = cumT + bursts(counter,count) + bursts(counter, count+1);
            count = count + 2;
        end
        numTones_played(counter) = n;

        %determine more popular pitch of played tones
        playedTones = freq{counter}(1:n);
        isLo = sum(playedTones == hd.loFreq);
        isHi = sum(playedTones == hd.hiFreq);
        coh_played(counter) = isHi/n;
        isH_played(counter) = isHi > isLo;
    else 
        rt = NaN;
        cur_choice = NaN;
    end 

    cur_f = isH(counter) + 1; %isH : 2 - high | 1 - low

    %Update choices list 
    timestamps = list{'Timestamps'}{'choices'};
    timestamps(counter) = timestamp;
    list{'Timestamps'}{'choices'} = timestamps;

    list{'Stimulus'}{'isH_played'} = isH_played;
    list{'Stimulus'}{'coh_played'} = coh_played;
    list{'Stimulus'}{'numTones_played'} = numTones_played;

    %assign choices and set target positions
    if strcmp(press, 'right')
        choice = 1;
        ensemble.setObjectProperty('xCenter', 5, target);
    elseif strcmp(press, 'left')
        choice = 2; 
        ensemble.setObjectProperty('xCenter', -5, target);
    elseif isempty(press)
        choice = NaN;
        if isH(counter)
            ensemble.setObjectProperty('xCenter', -5, target);
        else 
            ensemble.setObjectProperty('xCenter', 5, target);
        end 
    end 
    ensemble.setObjectProperty('isVisible', true, target);

    %add choice to list 
    choices(counter) = choice;
    list{'Input'}{'choices'} = choices; 

    %check whether or not choice was correct 
    if isempty(press)
        correct = NaN;
        string = 'Incorrect';
    elseif cur_f == choice
        correct = 1;
        string = 'Correct';
    else
        correct = 0;
        string = 'Incorrect';
    end
    
    %store data in list
    corrects = list{'Input'}{'corrects'};
    corrects(counter) = correct;
    list{'Input'}{'corrects'} = corrects;

    reac_times = list{'Input'}{'RT'};
    reac_times(counter) = rt;
    list{'Input'}{'RT'} = reac_times;

    fprintf('Trial %d complete. Choice: %s (%s). RT: %3.3f \n', ...
        counter, cur_choice, string, rt);
end 

function waitForCheckKey(list) 
    %Get list items 
    ui = list{'Input'}{'controller'}; 
    ui.flushData;

    %Initialize variable 
    press = '';

    %Wait for keypress to occur 
    while ~strcmp(press, 'continue')
        press = '';
        read(ui);
        [~, ~, eventname, ~] = ui.getHappeningEvent();
        if ~isempty(eventname) && length(eventname) == 1
            press = eventname;
        end 
    end 
end 

function playstim(list) 
    %Add current iteration to counter 
    counter = list{'Counter'}{'trial'};
    coh_list = list{'control'}{'cohLevels'};
    visualModes = list{'control'}{'visualModes'};

    hd = list{'Stimulus'}{'header'};
    
    %produce the stimulus
    [waveform, full_stimulus, f, h, bursts] = VisualTones(hd.loFreq, hd.hiFreq,...
        coh_list(counter), visualModes{counter});

    %player information 
    player = list{'Stimulus'}{'player'};    
    player.wave = full_stimulus; %both auditory and visual stimuli 
    %player.wave = [waveform; waveform]; %auditory only
    player.prepareToPlay;
    player.play;

    %log stimulus timestamps 
    stim_start = list{'Timestamps'}{'stim_start'};
    if ~isempty(player.playTime)
        stim_start(counter) = player.playTime;
    end 
    list{'Timestamps'}{'stim_start'} = stim_start;
    
    %store waveform
    waveforms = list{'Stimulus'}{'waveforms'};
    waveforms{counter} = waveform;
    list{'Stimulus'}{'waveforms'} = waveforms;

    %store frequency data 
    freq = list{'Stimulus'}{'freq'};
    freq{counter} = f;
    list{'Stimulus'}{'freq'} = freq;
    
    isH = list{'Stimulus'}{'isH'};
    isH(counter) = h;
    list{'Stimulus'}{'isH'} = isH;
    
    b = list{'Stimulus'}{'bursts'};
    b(counter, 1:length(bursts)) = bursts;
    list{'Stimulus'}{'bursts'} = b;
end 