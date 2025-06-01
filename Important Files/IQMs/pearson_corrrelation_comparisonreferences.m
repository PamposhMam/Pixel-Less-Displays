filename = 'Comparison_Metrics.xlsx';  % replace with actual filename
T = readtable(filename);
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

% Define numeric feature columns (exclude metadata)
exclude = {'Ranking', 'Score', 'NormalizedScore', 'Question'};
numericVars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
numericVars = setdiff(numericVars, exclude);

% Get list of unique questions
questions = unique(T.Question);
questions = questions(~isnan(questions));

% Compute best-worst deltas per question
metricDiffs = [];
for i = 1:length(questions)
    qData = T(T.Question == questions(i), :);
    best = qData(qData.Ranking == 1, :);
    worst = qData(qData.Ranking == 4, :);

    if ~isempty(best) && ~isempty(worst)
        delta = zeros(1, length(numericVars));
        for j = 1:length(numericVars)
            delta(j) = best.(numericVars{j}) - worst.(numericVars{j});
        end
        metricDiffs = [metricDiffs; delta];
    end
end

% Create table of deltas
diffTable = array2table(metricDiffs, 'VariableNames', numericVars);

% Compute correlation matrix
correlationMatrix = corr(table2array(diffTable), 'Type', 'Pearson');
correlationTable = array2table(correlationMatrix, 'VariableNames', numericVars, 'RowNames', numericVars);

% Show heatmap
figure;
heatmap(numericVars, numericVars, correlationMatrix, 'Colormap', parula, 'ColorLimits', [-1 1]);
title('Correlation Between IQA Feature Differences (Best - Worst)');

%%
T = readtable(filename);
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

% Define numeric features (exclude identifiers/scores)
exclude = {'Ranking', 'Score', 'NormalizedScore', 'Question'};
numericVars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
numericVars = setdiff(numericVars, exclude);

% Create binary label: 1 if image is best (Ranking == 1), else 0
T.WasBest = double(T.Ranking == 1);

% Compute Pearson correlation between WasBest and each metric
correlationResults = table;
for i = 1:length(numericVars)
    metric = numericVars{i};
    x = T.(metric);
    y = T.WasBest;
    if all(~isnan(x)) && all(~isnan(y))
        r = corr(x, y, 'Type', 'Pearson');
        correlationResults = [correlationResults;
            table(string(metric), r, 'VariableNames', {'Metric', 'PearsonWithBest'})];
    end
end

% Sort by absolute correlation
correlationResults = sortrows(correlationResults, 'PearsonWithBest', 'descend');
disp(correlationResults);