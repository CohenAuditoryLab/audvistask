function [task, list] = AudVisTask_v3(dispInd)
%% 05-22-2017 created by Brianna - Auditory Visual Task
%%Using VisualTones2, Create Stimulus, AddSpeakerCue
%%Sending triggers to ActiveX on PC using Serial Port
%%Stimuli are produced and saved on PC while all task controls are executed
%%here
%%Responses and timing occurs on Mac 
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

%% Settings for generating they sequence of conditions

% number visual modes
block_size = 3; 
% number of trials per visual mode
block_rep = 130; %1 %15 %50 %75
% possible visual values to select from
vis_vals = {'Low', 'High', 'None'}; %, 'All'};  %{'None', 'Low', 'High', 'All', 'Random'};

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
coherences = [0.00, 0.10, 0.25, 0.33, 0.40, 0.45, 0.50, 0.55, 0.60, 0.67, 0.75, 0.90, 1.00];

for i = 1:nTrials
    %add coherence as a condition
    tmp_coh = topsConditions('coh');
    tmp_coh.addParameter('cohLevel', coherences);
    likesCohLevel = topsFoundation();
    tmp_coh.addAssignment('cohLevel', likesCohLevel, '.', 'name');

    %randomly generate a coherence
    index = randsample(13,1);
    c = coherences(index);
    cohLevels(i, 1) = c;
end 

list{'control'}{'cohLevels'} = cohLevels;
%% Audio Settings

hd.loFreq = 1000; %hz      312.5 |  625 | 1250 | 2500 |  5000
hd.hiFreq = 2500; %hz     625   | 1250 | 2500 | 5000 | 10000
hd.toneDur = 50; %ms
hd.toneSOA = 10; %ms, actually poisson point process number centered around 10 
hd.trialDur = 4000; %ms
hd.fs = 24414; %samples/sec
hd.delay = 1.6; %sec

% INPUT PARAMETERS
responsewindow = hd.trialDur/1000 + hd.delay; %time allowed to respond = trial duration, s
                      %4 s for stimulus, .6 second for cue and 1 s delay
list{'Input'}{'responseWindow'} = responsewindow;

% OPEN SERIAL PORT CONNECTION TO PC
delete(instrfindall);
s = serial('/dev/tty.usbserial', 'BaudRate', 9600, 'DataBits', 8, ...
    'Parity', 'None', 'StopBits', 1, 'FlowControl', 'None', 'Terminator',...
    'LF');
fopen(s);

%% Time Variables

iti = 1; %seconds
list{'timing'}{'intertrial'} = iti; %intertrial interval
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
list{'Stimulus'}{'speaker'} = cell(nTrials,1);
list{'Stimulus'}{'bursts'} = zeros(nTrials,136);

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
blocklabel.string = sprintf('Block Number 1 of 3');
blocklabel.typefaceName = 'Calibri';
blocklabel.isVisible = false;
blocklabel.y = 5.5;

readyprompt1 = dotsDrawableText();
readyprompt1.string = 'This task will consist of 3 blocks.';
readyprompt1.fontSize = 32;
readyprompt1.typefaceName = 'Calibri';
readyprompt1.y = 4;
readyprompt1.isVisible = true;

readyprompt3 = dotsDrawableText();
readyprompt3.string = 'Each block takes ~10 minutes.';
readyprompt3.fontSize = 32;
readyprompt3.typefaceName = 'Calibri';
readyprompt3.y = 2;
readyprompt3.isVisible = true;

readyprompt = dotsDrawableText();
readyprompt.string = 'Feel free to take breaks between blocks.';
readyprompt.fontSize = 30;
readyprompt.typefaceName = 'Calibri';
readyprompt.isVisible = true;

infoprompt = dotsDrawableText();
infoprompt.string = 'Press A to play each sound and L or R to respond.';
infoprompt.fontSize = 24;
infoprompt.typefaceName = 'Calibri';
infoprompt.y = -2;
infoprompt.isVisible = true;

buttonprompt = dotsDrawableText();
buttonprompt.string = 'Ready? Press A to get started! ->';
buttonprompt.fontSize = 24;
buttonprompt.typefaceName = 'Calibri';
buttonprompt.y = -4;
buttonprompt.isVisible = true;

readyprompt2 = dotsDrawableText();
readyprompt2.string = 'Congratulations! You have completed the session.';
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
ready1 = ensemble.addObject(readyprompt1);
ready3 = ensemble.addObject(readyprompt3);
info = ensemble.addObject(infoprompt);

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
    'Hide', {hide [ready button ready1 ready3 info]}, {}, {}, 0, 'Finish'
    'Finish', {}, {}, {}, 0, '';};
prepareMachine.addMultipleStates(prepStates);

list{'control'}{'prepareMachine'} = prepareMachine;

% State Machine - used in maintask
mainMachine = topsStateMachine();
mainStates = {'name', 'entry', 'input', 'exit', 'timeout', 'next';
    'CheckReady', {@startTrial list block_rep s}, {}, {@waitForCheckKey list}, 0, 'Stimulus';
    'Stimulus', {@playstim list s}, {}, {@waitForChoiceKey list s}, 0, 'Feedback';
    'Feedback', {@showFeedback list}, {}, {}, 0, 'Exit';
    'Exit',{@finishTrial list}, {}, {}, iti,''};
mainMachine.addMultipleStates(mainStates);

list{'control'}{'mainMachine'} = mainMachine;

% End Machine - used in post-task
endMachine = topsStateMachine();
endStates = {'name', 'entry', 'input', 'exit', 'timeout', 'next';
    'Ready', {@startEndTask list s},      {},      {@waitForCheckKey list},     0,       'Hide';
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

function startEndTask(list, s)
    ensemble = list{'Graphics'}{'ensemble'};
    low = list{'Graphics'}{'low'};
    high = list{'Graphics'}{'high'};
    block = list{'Graphics'}{'block'};
    
    ensemble.setObjectProperty('isVisible', false, low);
    ensemble.setObjectProperty('isVisible', false, high);
    ensemble.setObjectProperty('isVisible', false, block);

    %prepare text + performance
    ready2 = list{'Graphics'}{'ready2'};
    button2 = list{'Graphics'}{'button2'};
    tmp_str = ensemble.getObjectProperty('string', ready2);
    ensemble.setObjectProperty('string', tmp_str, ready2);

    %make visible
    ensemble.setObjectProperty('isVisible', true, ready2);
    ensemble.setObjectProperty('isVisible', true, button2);
    
    %send cue to break outer while loop
    fprintf(s, 'done');
    
    %close serial port
    fclose(s);
end

function startTrial(list, block_rep, s)
    %clear last trial data
    ui = list{'Input'}{'controller'};
    ui.flushData();

    %increment counter to label trial
    counter = list{'Counter'}{'trial'};
    counter = counter + 1;
    list{'Counter'}{'trial'} = counter;
    
    coh_list = list{'control'}{'cohLevels'};
    visualModes = list{'control'}{'visualModes'};
    hd = list{'Stimulus'}{'header'};
    
    %decide which speaker gets target and which gets masker 
    r = 2 * rand;
    if r <= 1
        chosen = 'ONE';
    else 
        chosen = 'TWO';
    end 
    %store target speaker 
    speaker = list{'Stimulus'}{'speaker'};
    speaker{counter} = chosen;
    list{'Stimulus'}{'speaker'} = speaker;
    
    %Create string to send to serial port
    data = ['START.', num2str(hd.loFreq), '.', num2str(hd.hiFreq), '.', ...
        num2str(coh_list(counter) * 100), '.', visualModes{counter}, '.', ...
        chosen, '.STOP\n'];
    
    fprintf(s, data);   

    %set high and low labels to be visible 
    ensemble = list{'Graphics'}{'ensemble'};
    low = list{'Graphics'}{'low'};
    high = list{'Graphics'}{'high'};
    block = list{'Graphics'}{'block'};
    
    ensemble.setObjectProperty('isVisible', false, block);
   
    %Wait for cue that data was received 
    cue = '';
    while length(cue) < 1
        cue = fscanf(s, '%s');
    end 
    
    %Start adding things to the screen
    b = int16(ceil((counter/block_rep)));
    block = dotsDrawableText();
    block.string = sprintf('Block %d of 3', b);
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

function string = waitForChoiceKey(list, s)
    %Get list items
    counter = list{'Counter'}{'trial'};
    ensemble = list{'Graphics'}{'ensemble'};
    target = list{'Graphics'}{'target'};
    ui = list{'Input'}{'controller'};
    hd = list{'Stimulus'}{'header'};
    stim_start = list{'Timestamps'}{'stim_start'};
    responsewindow = list{'Input'}{'responseWindow'};
    choices = list{'Input'}{'choices'};
    coh_list = list{'control'}{'cohLevels'};

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
           	fprintf(s, '%s\n', 'no');
            choice = NaN;
            timestamp = NaN;
            break
        end 

        %Check for button press 
        press = '';
        read(ui);
        [~, ~, eventname, ~] = ui.getHappeningEvent();
        events = cell(length(choices), 1);

        if ~isempty(eventname) && length(eventname) == 1
            press = eventname;
            events(counter) = press;
            %stop the stimulus once a response is detected 
            %avoid error by only stopping if left or right is pressed
            if ~strcmp(events{counter}, 'continue')
                fprintf(s, '%s\n', 'no');
            end

            %get the timestamp of the stimulus stop time 
            stim_stop = list{'Timestamps'}{'stim_stop'};
            if length(GetSecs) == 1
                stim_stop(counter) = GetSecs;
            end
            list{'Timestamps'}{'stim_stop'} = stim_stop;
        end 
    end

    %Get timestamp of button press time 
    if ~isempty(press)
        timestamp = ui.history;
        %to ensure timestamp from a pressed key/button
        timestamp = timestamp(timestamp(:, 2) > 1, :);
        timestamp = timestamp(end);

        %calculate reaction time 
        rt = (timestamp - stim_start(counter)) * 1000; %ms
        delay = hd.delay * 1000; %ms
        rt = rt - delay;
        if rt <= 0 
            rt = NaN;
        end 
        %record current choice 
        cur_choice = press{1};
    else 
        rt = NaN;
        cur_choice = NaN;
    end 

    cur_f = (coh_list(counter) > 0.50) + 1; %isH : 2 - high | 1 - low

    %Update choices list 
    timestamps = list{'Timestamps'}{'choices'};
    timestamps(counter) = timestamp;
    list{'Timestamps'}{'choices'} = timestamps;

    %assign choices and set target positions
    if strcmp(press, 'right')
        choice = 1;
        ensemble.setObjectProperty('xCenter', 5, target);
    elseif strcmp(press, 'left')
        choice = 2; 
        ensemble.setObjectProperty('xCenter', -5, target);
    elseif isempty(press)
        choice = NaN;
        if (coh_list(counter) > 0.50)
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

function playstim(list, s) 
    %Add current iteration to counter 
    counter = list{'Counter'}{'trial'};
    
    %Play stimulus
    fprintf(s, '%s\n', 'go');
    
    %log stimulus timestamps 
    stim_start = list{'Timestamps'}{'stim_start'};
    stim_start(counter) = GetSecs;
    list{'Timestamps'}{'stim_start'} = stim_start;
end 