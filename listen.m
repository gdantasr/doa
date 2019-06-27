function listen(audio_struct, fs)
% Plays audio data stored in audio_struct

audio_data = audio_struct.soundPressure;
audio_data = [audio_data(:, 4) audio_data(:, 10)]; % L and R channels
soundsc(audio_data, fs) 

end

