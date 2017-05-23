classdef dotsPlayableWave_2Channel < dotsPlayable
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
        function self = dotsPlayableWave_2Channel()
            self = self@dotsPlayable();
            InitializePsychSound
            
        end
        
        % Compute a sinusoidal wavform to play.
        function prepareToPlay(self)
            self.player = PsychPortAudio('Open', [], [], 1, self.sampleFrequency, 2);
            self.waveform = self.wave*self.intensity;
%             self.player = audioplayer(self.wave, ...
%                 self.sampleFrequency, self.bitsPerSample);
        end
        
        % Play the tone.
        function play(self)
            if ~isempty(self.wave)
                % play is async, playblocking would be sync
                PsychPortAudio('FillBuffer', self.player, self.wave);
                self.playTime = PsychPortAudio('Start', self.player, 1, GetSecs, 1, inf, 0);
            end
        end
        % Stop the tone
        function stop(self)
            [~, ~, ~, self.stopTime] = PsychPortAudio('Stop', self.player, 0, 1);
            PsychPortAudio('Close', self.player);
        end
    end
end