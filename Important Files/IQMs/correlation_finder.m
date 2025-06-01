%% Load data correctly using readtable
T = readtable('Solo_IQMFeatures.xlsx');

% Extract feature labels (column headers from 7th column onward)
data_labels = T.Properties.VariableNames(7:end);

%% Extract data for each ranking (1 to 4)
data1 = T{T.Ranking == 1, 7:end};
data2 = T{T.Ranking == 2, 7:end};
data3 = T{T.Ranking == 3, 7:end};
data4 = T{T.Ranking == 4, 7:end};

datasets = {data1, data2, data3, data4};
dataset_names = {'Dataset1', 'Dataset2', 'Dataset3', 'Dataset4'};

%% Compute Weighted Dataset (5th dataset)
question_ids = unique(T.Question);
num_questions = numel(question_ids);
num_features = numel(T.Properties.VariableNames(7:end));

weighted_feature_data = NaN(num_questions, num_features);  % preallocate as NaN

for q = 1:num_questions
    qid = question_ids(q);
    
    % Rows for this question where NormalizedScore is finite (i.e. not NaN or Inf)
    rows = T.Question == qid & isfinite(T.NormalizedScore);
    scores = T.NormalizedScore(rows);
    features = T{rows, 7:end};

    % Check everything is numeric and compatible
    if isempty(scores) || all(scores == 0) || any(isnan(scores)) || size(features, 1) ~= numel(scores)
        warning('Skipping Question %d due to invalid scores/features', qid);
        continue;
    end

    % Safe weighted average
    weighted_feature_data(q, :) = sum(features .* scores, 1) ./ sum(scores);
end


datasets{end+1} = weighted_feature_data;
dataset_names{end+1} = 'WeightedDataset';

%% Compute relative datasets (Ref vs Best/Worst, Best vs Worst)
ref_vs_best   = zeros(num_questions, num_features);
ref_vs_worst  = zeros(num_questions, num_features);
best_vs_worst = zeros(num_questions, num_features);

for q = 1:num_questions
    qid = question_ids(q);
    group = T(T.Question == qid, :);
    scores = group.Score;
features = group{:, 7:end};

ref_idx = find(scores == 0);
if numel(ref_idx) ~= 1
    warning('Question %d: invalid reference count (%d)', qid, numel(ref_idx));
    continue;
end

score_candidates = scores;
score_candidates(ref_idx) = NaN;

[~, best_idx] = max(score_candidates);
[~, worst_idx] = min(score_candidates);

if any(isnan([best_idx, worst_idx]))
    warning('Question %d: could not find valid best/worst', qid);
    continue;
end

ref_features   = features(ref_idx, :);
best_features  = features(best_idx, :);
worst_features = features(worst_idx, :);

end

datasets{end+1} = ref_vs_best;
dataset_names{end+1} = 'RefVsBest';
datasets{end+1} = ref_vs_worst;
dataset_names{end+1} = 'RefVsWorst';
datasets{end+1} = best_vs_worst;
dataset_names{end+1} = 'BestVsWorst';

%% Loop through and compute correlation & regression stats for feature pairs
for d = 5:length(datasets)
    current_data = datasets{d};
    current_name = dataset_names{d};
    num_features = size(current_data, 2);

    output_dir = fullfile(pwd, current_name);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    for i = 1:num_features-1
        for j = i+1:num_features
            x = current_data(:, i);
            y = current_data(:, j);

            % Remove invalid entries
            valid = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);
            x = x(valid);
            y = y(valid);

            if numel(x) > 2
                r = corr(x, y, 'Type', 'Pearson');
                p = polyfit(x, y, 1);
                y_fit = polyval(p, x);
                residuals = y - y_fit;

                stats = struct();
                stats.slope = p(1);
                stats.intercept = p(2);
                stats.mean_residual = mean(residuals);
                stats.std_residual = std(residuals);
                stats.correlation_coefficient = r;
            else
                r = NaN;
                y_fit = NaN(size(x));
            end

            % Plot
            figure;
            scatter(x, y, 'filled');
            hold on;
            plot(x, y_fit, 'r', 'LineWidth', 2);
            title(sprintf('%s vs %s (r=%.2f)', data_labels{i}, data_labels{j}, r));
            xlabel(data_labels{i});
            ylabel(data_labels{j});
            grid on;

            % Save plot
            plot_filename = fullfile(output_dir, sprintf('%s_vs_%s.png', data_labels{i}, data_labels{j}));
            saveas(gcf, plot_filename);
            close(gcf);

            % Save stats
            stats_filename = fullfile(output_dir, sprintf('%s_vs_%s_Stats.txt', data_labels{i}, data_labels{j}));
            fileID = fopen(stats_filename, 'w');
            fprintf(fileID, '%s vs %s\n', data_labels{i}, data_labels{j});
            fprintf(fileID, 'Correlation Coefficient (r): %.4f\n', stats.correlation_coefficient);
            fprintf(fileID, 'Linear Regression Slope: %.4f\n', stats.slope);
            fprintf(fileID, 'Linear Regression Intercept: %.4f\n', stats.intercept);
            fprintf(fileID, 'Mean Residual: %.4f\n', stats.mean_residual);
            fprintf(fileID, 'Standard Deviation of Residuals: %.4f\n', stats.std_residual);
            fclose(fileID);
        end
    end
end
