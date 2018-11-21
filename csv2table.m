function tab = csv2table( fileName )
%   L� dados de arquivo .csv formatado em json e retorna uma tabela 
%   com os dados

fileID = fopen(fileName); % Abre arquivo
columnID = fgetl(fileID); % L� primeira linha contendo o identificador das colunas
columnID = strsplit(columnID,', ');
nCols = length(columnID);
format = repmat('%s', [1 nCols]);

data = textscan(fileID, format, 'delimiter',','); % L� restante do arquivo. Cada coluna do arquivo � guardada como uma c�lula em data

% Guarda cada coluna do arquivo em uma coluna de dataCell
dataCell = cell(length(data{1,1}), length(data) ); 
for col = 1:size(data, 2)
    dataCell(:, col) = data{1,col};
end

tab = cell2table(dataCell,...
    'VariableNames',columnID);

% v = cellfun(@str2num, T.kph);
% t = cellfun(@str2num, T.secs);
% ac = gradient(v);
% plot(t, v, 'b', t, ac, 'r')

end

