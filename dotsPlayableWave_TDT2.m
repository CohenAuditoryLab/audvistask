classdef dotsPlayableWave_TDT2 < dotsPlayable
    % @class dotsPlayableTone
    % Play a pure sinusoudal tone
    
    properties
        % waveform vector
        wave;
        playTime;
        stopTime;
    end
    
    properties (SetAccess = protected)
        % Matlab audioplayer object -> Psychtoolbox audioplayer object
        player;
    end
    
    methods
        % Constructor takes no arguments.
        function self = dotsPlayableWave_TDT2()
            self = self@dotsPlayable();
            %InitializePsychSound
            
            %%%%%CHANGE THIS PATH
            % set location of RCXfile
            RCXFile = fullfile('C:\','Users','Amanda Yung','Documents',...
                'MATLAB','AudCategoryTask','audCat_StimGen.rcx');
            
            % setup connection with TDT
            RP = actxcontrol('RPco.x',[5 5 26 26]);
            
            if RP.ConnectRX6('GB', 1)
                disp('Connected to RX6!');
            else
                disp('Unable to connect to RX6');
            end
            
            % load rcx file
            RP.LoadCOF(RCXFile);
            RP.Run;
        end
        
        % Compute a sinusoidal wavform to play.
        function prepareToPlay(self)
            self.waveform = self.wave*self.intensity;
            RP.WriteTagVEX('datain2', 0, 'F32', self.waveform(1,:));
            RP.WriteTagVEX('lightin2', 0, 'F32', self.waveform(2,:));
            drawnow;
        end
        
        % Play the tone.
        function play(self)
            RP.Run();
            if ~isempty(self.wave)
                %trigger for on
                RP.SoftTrg(4);
                RP.SoftTrg(3);
            end
        end
        % Stop the tone
        function stop(self)
            %trigger for off 
            RP.SoftTrg(4);
            RP.Halt;
            RP.ClearCOF();
        end
    end
end