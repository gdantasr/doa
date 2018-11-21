function [dMin, iMin] = getCentralPoint(passbyFileName)
%   UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Get array position. Average from measurements
arrayFileName = 'pos_array.csv';
tab = csv2table(arrayFileName);
arrayLat = cellfun(@str2num, tab.lat);
arrayLon = cellfun(@str2num, tab.lon);
avgArrayLat = mean(arrayLat);
avgArrayLon = mean(arrayLon);
loc1 = [avgArrayLat avgArrayLon];

% Get vehicle positions
tab = csv2table(passbyFileName);
passbyLat = cellfun(@str2num, tab.lat);
passbyLon = cellfun(@str2num, tab.lon);

% Compute array-vehicle distante for each measure (1 per sec.)
d = zeros(1, length(passbyLat));
for i = 1 : length(passbyLat)
    loc2 = [passbyLat(i) passbyLon(i)];     % Get current vehicle lat./lon. 
    km = haversine(loc1, loc2);             % Lat./lon. coordinates to distance in kilometres
    d(i) = km*1000;                         % Distance in metres
end

[dMin, iMin] = min(d);

end