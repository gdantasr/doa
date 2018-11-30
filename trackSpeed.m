function [ Vshifted ] = trackSpeed( filename )
%   Summary of this function goes here
%   Detailed explanation goes here

    % Speed tracking
    passbyTab = csv2table([filename, '.csv']);
    passbyLat = cellfun(@str2num, passbyTab.lat);
    passbyLon = cellfun(@str2num, passbyTab.lon);
    passbySpeed = cellfun(@str2num, passbyTab.kph);
    passbyTime = cellfun(@str2num, passbyTab.secs);
    [~, iMin] = getCentralPoint(passbyLat, passbyLon);   
    Vq = interp1(passbyTime, passbySpeed, t, 'v5cubic');
    Vq(isnan(Vq)) = 0.0;
    [~, iCentral] = min(abs(iMin - t));
    [~, nShift] = min(abs(t90 - t));
    Vshifted = circshift(Vq, [0, nShift - iCentral]);

end

