filename = 'Solo_IQMFeatures.xlsx';
T = readtable(filename);

% Clean column names
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

% Select relevant numeric features
exclude = {'Ranking', 'Score', 'NormalizedScore', 'Quadrant', 'Question'};
numericVars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
numericVars = setdiff(numericVars, exclude);

% Unique questions
questions = unique(T.Question);
questions = questions(~isnan(questions));  % filter out ghost NaNs

% Build clean delta matrix
metricDiffs = [];
for i = 1:length(questions)
    qData = T(T.Question == questions(i), :);

    % 🔍 Filter: Only include questions where there are 4 Score entries,
    % not five.
    scoreValues = qData.Score;
    if sum(scoreValues ~= 0) == 4
        fprintf('❌ Question %d skipped: Reference Image Present\n', questions(i));
        continue;
    end

    best = qData(qData.Ranking == 1, :);
    worst = qData(qData.Ranking == 0, :); %just change this to get the base ref

    if ~isempty(best) && ~isempty(worst)
        delta = zeros(1, length(numericVars));
        for j = 1:length(numericVars)
            delta(j) = best.(numericVars{j}) - worst.(numericVars{j});
        end
        metricDiffs = [metricDiffs; delta];
    else
        fprintf('⚠️ Question %d skipped: missing best or worst\n', questions(i));
    end
end


% Create a table from differences
diffTable = array2table(metricDiffs, 'VariableNames', numericVars);

% Compute correlation matrix between feature deltas
correlationMatrix = corr(table2array(diffTable), 'Type', 'Pearson');

% Convert to table for clarity
correlationTable = array2table(correlationMatrix, 'VariableNames', numericVars, 'RowNames', numericVars);

% Display as heatmap
figure;
heatmap(numericVars, numericVars, correlationMatrix, 'Colormap', parula, 'ColorLimits', [-1 1]);
title('Correlation Between Feature Differences Based On Relative Difference of Metrics Between the Best Image and its Reference Image');

%%
% Set the folder containing your .fig files
folderPath = "C:\Users\Pamposh\Desktop\Wilkinson Project Stuff\LiveHologramsForProject\IQM_Stuff";  % <- CHANGE THIS

% Get all .fig files in the folder
figFiles = dir(fullfile(folderPath, '*.fig'));

% Loop through each .fig file
for k = 1:length(figFiles)
    figName = figFiles(k).name;
    figPath = fullfile(folderPath, figName);

    % Open the .fig file
    fig = openfig(figPath, 'invisible');

    % Create output filename (.jpg)
    [~, name, ~] = fileparts(figName);
    jpgName = fullfile(folderPath, [name, '.jpg']);

    % Save as JPEG
    saveas(fig, jpgName);

    % Close the figure
    close(fig);
end

disp('All .fig files converted to .jpg.');
