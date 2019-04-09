function [ vShifted ] = track_speed( filename, tAudio, tCentralAudio )
%   Speed tracking
%   filename - Nome do arquivo CSV
%   tAudio - Vetor de tempo dos dados de áudio
%   tCentralAudio - Instante da passagem do carro em frente ao array

    % Read GPS data from csv table
    passbyTab = csv2table([filename, '.csv']);      
    passbyLat = cellfun(@str2num, passbyTab.lat);           % Latitude
    passbyLon = cellfun(@str2num, passbyTab.lon);           % Longitude
    passbySpeed = cellfun(@str2num, passbyTab.kph);         % Speed
    passbyTime = cellfun(@str2num, passbyTab.secs);         % Time
    
    % Interpolate GPS measures
    maxtQ = min([length(passbyTime), tAudio(end)]);
    tQ = 0 : mean(diff(tAudio)) : maxtQ; % Time vector with audio sampling frequency
    lonQ = interp1(passbyTime, passbyLon, tQ, 'linear');
    lonQ(isnan(lonQ)) = 0.0;
    latQ = interp1(passbyTime, passbyLat, tQ, 'linear');
    latQ(isnan(latQ)) = 0.0;
    vQ = interp1(passbyTime, passbySpeed, tQ, 'nearest');
    vQ(isnan(vQ)) = 0.0;
    
    % Compare vehicle and array GPS positions and get the minimum car-array
    % distance point
    [~, nCentralGPS] = getCentralPoint(latQ, lonQ);
    
    % Shift speed vector to align GPS and audio central points
    %[~, iCentral] = min(abs(nCentralGPS - tAudio));
    [~, nCentralAudio] = min(abs(tCentralAudio - tAudio));
    vQ = [ vQ zeros(1, length(tAudio) - length(tQ)) ];
    nDelay = nCentralAudio - nCentralGPS;
    vShifted = circshift(vQ, [0, nDelay]);                  % Circular shift
    % Ignore circular data
    if nDelay > 0
        vShifted = [zeros(1, nDelay) vShifted(nDelay+1 : end)];
    elseif nDelay < 0
        vShifted = [vShifted(1 : end + nDelay) zeros(1, abs(nDelay))];
    end
        
end

