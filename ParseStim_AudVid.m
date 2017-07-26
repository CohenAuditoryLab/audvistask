function ParseStim_AudVid
%% Set Up RX6 and ActiveX Control

% set location of RCXfile
RCXFile = fullfile('C:\','work','AudVis_Task','AudVis_Circuit.rcx');

% setup connection with TDT
RP = actxcontrol('RPco.x',[5 5 26 26]);

% attempt to connect to the RX6
if RP.ConnectRX6('GB', 1)
    disp('Connected to RX6!');
else
    disp('Unable to connect to RX6');
end

% load rcx file
if RP.LoadCOF(RCXFile);
    disp('Circuit loaded!');
else
    disp('Failed to load circuit.')
end

% begin running circuit
if RP.Run;
    disp('Circuit running.');
else
    disp('Error running circuit.');
end

%% Initialize Serial Port to Receive Data

delete(instrfindall);
s = serial('COM1', 'BaudRate', 9600, 'DataBits', 8);
fopen(s);

% set up data structures to store stimulus info
nTrials = 3;
waveforms = cell(nTrials,1);
maskers = cell(nTrials,1);
freq = cell(nTrials,1);
isH = zeros(nTrials,1);
counter = 0;

%% Begin session

sessionRun = 1;

while sessionRun == 1
    %% Wait for data from serial port
    disp('===================== Start of Session ========================');
    
    checkingForStim = true;
    
    while checkingForStim
        counter = counter + 1;
        data = fscanf(s);
        
        if length(data) >= 6
            try
                [l, h, c, v, spk] = parse(data);
            catch
                disp('Stimulus string formatted incorrectly.')
                continue
            end
            if paramFlag==1
                checkingForStim = false;
            end
        end
    end
    
    %% Construct and upload stimulus to RX6
    [target_auditory, target_stimulus, frequencies, isHigh,...
        masker_auditory, masker_stimulus] = CreateStimulus(l, h, c, v);
    
    [new_target, new_masker] = AddSpeakerCue(target_stimulus, masker_stimulus);
    
    % save stimulus info
    waveforms(counter) = target_auditory;
    maskers(counter) = masker_auditory;
    freq(counter) = frequencies;
    isH(counter) = isHigh;
    
    %% Trigger stimulus start
    
    abortFlag = findTrialStart;
    % if received the 'no' signal
    if abortFlag == 1
        % break out of trial and begin while loop over again
        disp('Aborting Trial');
        continue;
    end
    
    % LOAD stimuli to speakers/LEDs
    % if the first speaker is the target
    if strcmp(spk, 'ONE')
        RP.WriteTagVEX('datain1', 0, 'F32', new_target(1,:));
        RP.WriteTagVEX('lightin1', 0, 'F32', new_target(2,:));
        RP.WriteTagVEX('datain2', 0, 'F32', new_masker(1,:));
        RP.WriteTagVEX('lightin2', 0, 'F32', new_masker(2,:));
    end
    % if the second speaker is the target
    if strcmp(spk, 'TWO')
        RP.WriteTagVEX('datain1', 0, 'F32', new_masker(1,:));
        RP.WriteTagVEX('lightin1', 0, 'F32', new_masker(2,:));
        RP.WriteTagVEX('datain2', 0, 'F32', new_target(1,:));
        RP.WriteTagVEX('lightin2', 0, 'F32', new_target(2,:));
    end
    
    % Play the stimulus
    RP.SoftTrg(1); %Ch1
    RP.SoftTrg(3); %Ch2
    
    % If you receive the stop stimulus
    abortFlag = findTrialStart;
    if abortFlag == 1
        % break out of trial and begin while loop over again
        disp('Aborting Trial');
        RP.SoftTrg(2); %Ch1
        RP.SoftTrg(4); %Ch2
        RP.ClearCOF();
        continue;
    end
end

% Close serial port and RX6
RP.Halt;
fclose(s);

%% Save stimulus information to file

stim_table = table((1:counter)', waveforms, maskers, freq, isH, 'VariableNames',...
    {'trial', 'waveform','masker','tone_frequencies','isHigh'});
data_folder = fullfile('C:', 'work', 'AudVis_Task', 'StimData');
c = clock;
save_filename = ['Stimulus_Data_', num2str(c(1)), num2str(c(2)), ...
    num2str(c(3)), '_', num2str(c(4)), num2str(c(5))];
save([data_folder save_filename '_table.mat'], 'stim_table');

%% Helper function for parsing data

    function [low, high, coh, vis, spk] = parse(d)
        % Suppose string has form:
        % START.LLLL.HHHH.C.CC.VVV(V).SSS.STOP
        
        % Initialize values
        low = '';
        high = '';
        coh = '';
        vis = '';
        spk = '';
        paramFlag = 0;
        
        % Parse string
        tmpstr1 = strsplit(d,'START.');
        tmpstr2 = strsplit(tmpstr1{2},'.STOP');
        alldat = strsplit(tmpstr2{1},'.');
        low = str2double(alldat{1});
        high = str2double(alldat{2});
        coh = str2double([alldat{3} alldat{4}])/100;
        vis = alldat{5};
        spk = alldat{6};
        paramFlag = 1;
    end

%% Monitor start signal
    function abortFlag = findTrialStart
        checkAbortFlag = 1;
        while checkAbortFlag == 1
            abortFlag = [];
            
            dat2 = fscanf(sPort);
            if strfind(dat2, 'no')
                abortFlag = 1;
                checkAbortFlag = 0;
            end
            if strfind(dat2, 'go')
                checkAbortFlag = 0;
            end
        end
    end

end