% Create comparative tables between DoA methods

% load('v_estimated.mat');
% load('wb_estimated.mat');
vReal = cell2mat(v);

carID = zeros(length(fileNames), 1);
vGroup = zeros(length(fileNames), 1);

% Define function for statistical calculation on data
mystats = @(x) [mean(x) std(x)];

for fileID = 1:length(fileNames)
    % Check car ID from filename
    carID(fileID) = strfind('bfjm', fileNames{fileID}(1));
    
    % Same speed grouping
    if v{fileID} >= 30 && v{fileID} < 40
        vGroup(fileID) = 1;
    elseif v{fileID} >= 40 && v{fileID} < 50
        vGroup(fileID) = 2;
    elseif v{fileID} >= 50 && v{fileID} < 60
        vGroup(fileID) = 3;
    else % v >=60
        vGroup(fileID) = 4;
    end
end
carGroup = findgroups(carID);  % Same car grouping
wbTable = zeros(5, 8);
vTable = zeros(5, 8);
for col = 1:5
    stat = splitapply(mystats, wb_estimated(:,col), carGroup)';
    stat = stat(:)';
    wbTable(col, :) = stat;
    vError = abs(v_estimated(:,col) - vReal);
    stat = splitapply(mystats, vError, vGroup)';
    stat = stat(:)';
    vTable(col, :) = stat;
end

methodTable = cell2table(upper(methodNames)', 'VariableNames', {'method'}); % Method names column
% Round value to 2 decimal places and add method names columns
wbTable = [methodTable table(round(wbTable*100)/100)];
vTable = [methodTable table(round(vTable*100)/100)];

writetable( wbTable, [filePath, 'wbTable.csv']);
writetable( vTable, [filePath, 'vTable.csv']);