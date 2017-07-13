classdef dotsPlayableWave_TDT < dotsPlayable
    % @class dotsPlayableTone
    % Play a pure sinusoudal tone
    
    properties
        % waveform vector
        wave;
        masker;
        playTime;
        stopTime;
        RP;
    end
    
    properties (SetAccess = protected)
        % Matlab audioplayer object -> Psychtoolbox audioplayer object
        player;
    end
    
    methods
        % Constructor takes no arguments.
        function self = dotsPlayableWave_TDT()
            self = self@dotsPlayable();
            
            % set location of RCXfile
            RCXFile = fullfile('C:\','work','AudVis_Task','AudVis_Circuit.rcx');
            
            % setup connection with TDT
            self.RP = actxcontrol('RPco.x',[5 5 26 26]);
            
            if self.RP.ConnectRX6('GB', 1)
                disp('Connected to RX6!');
            else
                disp('Unable to connect to RX6');
            end
            
            % load rcx file
            self.RP.LoadCOF(RCXFile);
            self.RP.Run;
        end
        
        % Compute a sinusoidal wavform to play.
        function prepareToPlay(self)
            self.waveform = self.wave*self.intensity;
            self.masker = self.masker*self.intensity;
            %RP = actxcontrol('RPco.x',[5 5 26 26]);
            self.RP.WriteTagVEX('datain1', 0, 'F32', self.waveform(1,:));
            self.RP.WriteTagVEX('lightin1', 0, 'F32', self.waveform(2,:));
            self.RP.WriteTagVEX('datain2', 0, 'F32', self.masker(1,:));
            self.RP.WriteTagVEX('lightin2', 0, 'F32', self.masker(2,:));            
        end
        
        % Play the tone.
        function play(self)
            if ~isempty(self.wave)
                %trigger for on
                self.RP.SoftTrg(1); %Ch1
                self.RP.SoftTrg(3); %Ch2
            end
        end
        
        % Stop the tone
        function stop(self)
            %trigger for off 
            self.RP.SoftTrg(2); %Ch1
            self.RP.SoftTrg(4); %Ch2
            self.RP.Halt();
            self.RP.ClearCOF();
        end
    end
end