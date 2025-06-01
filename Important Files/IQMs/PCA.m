% Load data from Excel sheet
[data, txt, raw] = xlsread('Solo_IQMFeatures.xlsx'); % Replace with your file name

% Extract relevant columns and remove irrelevant ones (2nd and 4th columns)
data = data(:, [1, 2:end]); % Retain relevant features (ranking in the first column)
data_labels = raw(1, 2:end); % First row (excluding the first column) -> feature labels

% Check for unique rankings in the 3rd column
rankings = unique(data(:, 3)); % Extract unique ranking values from the 3rd column


 % Extract data for ranking 1
 data1 = data(data(:, 2) == 1, 4:end); % Exclude first three columns

 % Extract data for ranking 2
 data2 = data(data(:, 2) == 2, 4:end); % Exclude first three columns
    
% Extract data for ranking 3
 data3 = data(data(:, 2) == 3, 4:end); % Exclude first three columns
  
% Extract data for ranking 4
data4 = data(data(:, 2) == 4, 4:end); % Exclude first three columns


% Normalize the data (optional but recommended)
% Subtract the mean and divide by the standard deviation
data_norm1 = (data1 - mean(data1)) ./ std(data1);
data_norm2 = (data2 - mean(data2)) ./ std(data2);
data_norm3 = (data3 - mean(data3)) ./ std(data3);
data_norm4 = (data4 - mean(data4)) ./ std(data4);

% Perform PCA
[coeff1, score1, latent1, tsquared1, explained1] = pca(data_norm1);
[coeff2, score2, latent2, tsquared2, explained2] = pca(data_norm2);
[coeff3, score3, latent3, tsquared3, explained3] = pca(data_norm3);
[coeff4, score4, latent4, tsquared4, explained4] = pca(data_norm4);

% Display the percentage of variance explained by each principal component for all rankings
fprintf('Explained Variance (%%) - Ranking 1:\n');
disp(explained1);
fprintf('Explained Variance (%%) - Ranking 2:\n');
disp(explained2);
fprintf('Explained Variance (%%) - Ranking 3:\n');
disp(explained3);
fprintf('Explained Variance (%%) - Ranking 4:\n');
disp(explained4);

% Visualization and saving for Ranking 1
figure;
scatter(score1(:,1), score1(:,2), 'filled');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
title('PCA - Ranking 1');
grid on;
saveas(gcf, 'PCA_Ranking_1_Scatter.png'); % Save scatter plot

figure;
pareto(explained1);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot - Ranking 1');
saveas(gcf, 'PCA_Ranking_1_Scree.png'); % Save scree plot

% Visualization and saving for Ranking 2
figure;
scatter(score2(:,1), score2(:,2), 'filled');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
title('PCA - Ranking 2');
grid on;
saveas(gcf, 'PCA_Ranking_2_Scatter.png'); % Save scatter plot

figure;
pareto(explained2);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot - Ranking 2');
saveas(gcf, 'PCA_Ranking_2_Scree.png'); % Save scree plot

% Visualization and saving for Ranking 3
figure;
scatter(score3(:,1), score3(:,2), 'filled');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
title('PCA - Ranking 3');
grid on;
saveas(gcf, 'PCA_Ranking_3_Scatter.png'); % Save scatter plot

figure;
pareto(explained3);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot - Ranking 3');
saveas(gcf, 'PCA_Ranking_3_Scree.png'); % Save scree plot

% Visualization and saving for Ranking 4
figure;
scatter(score4(:,1), score4(:,2), 'filled');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
title('PCA - Ranking 4');
grid on;
saveas(gcf, 'PCA_Ranking_4_Scatter.png'); % Save scatter plot

figure;
pareto(explained4);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot - Ranking 4');
saveas(gcf, 'PCA_Ranking_4_Scree.png'); % Save scree plot


% Assuming data_labels are the feature names from your first row
% Get the top contributing features to each principal component
% Get the indices of the features with the highest absolute values for PC1 and PC2
[~, idx1] = sort(abs(coeff1(:,1)), 'descend');
[~, idx2] = sort(abs(coeff1(:,2)), 'descend');

% Create labels based on the most influential features
pc1_labels = data_labels(idx1(1:5)); % Top 5 features for PC1
pc2_labels = data_labels(idx2(1:5)); % Top 5 features for PC2

% Display the percentage of variance explained by each principal component for all rankings
fprintf('Explained Variance (%%) - Ranking 1:\n');
disp(explained1);
fprintf('Explained Variance (%%) - Ranking 2:\n');
disp(explained2);
fprintf('Explained Variance (%%) - Ranking 3:\n');
disp(explained3);
fprintf('Explained Variance (%%) - Ranking 4:\n');
disp(explained4);

% Visualization and saving for Ranking 1
figure;
scatter(score1(:,1), score1(:,2), 'filled');
xlabel(['Principal Component 1 (' strjoin(pc1_labels, ', ') ')']);
ylabel(['Principal Component 2 (' strjoin(pc2_labels, ', ') ')']);
title('PCA - Ranking 1');
grid on;
saveas(gcf, 'PCA_Ranking_1_Scatter.png'); % Save scatter plot

figure;
pareto(explained1);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot - Ranking 1');
saveas(gcf, 'PCA_Ranking_1_Scree.png'); % Save scree plot

% Visualization and saving for Ranking 2
figure;
scatter(score2(:,1), score2(:,2), 'filled');
xlabel(['Principal Component 1 (' strjoin(pc1_labels, ', ') ')']);
ylabel(['Principal Component 2 (' strjoin(pc2_labels, ', ') ')']);
title('PCA - Ranking 2');
grid on;
saveas(gcf, 'PCA_Ranking_2_Scatter.png'); % Save scatter plot

figure;
pareto(explained2);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot - Ranking 2');
saveas(gcf, 'PCA_Ranking_2_Scree.png'); % Save scree plot

% Visualization and saving for Ranking 3
figure;
scatter(score3(:,1), score3(:,2), 'filled');
xlabel(['Principal Component 1 (' strjoin(pc1_labels, ', ') ')']);
ylabel(['Principal Component 2 (' strjoin(pc2_labels, ', ') ')']);
title('PCA - Ranking 3');
grid on;
saveas(gcf, 'PCA_Ranking_3_Scatter.png'); % Save scatter plot

figure;
pareto(explained3);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot - Ranking 3');
saveas(gcf, 'PCA_Ranking_3_Scree.png'); % Save scree plot

% Visualization and saving for Ranking 4
figure;
scatter(score4(:,1), score4(:,2), 'filled');
xlabel(['Principal Component 1 (' strjoin(pc1_labels, ', ') ')']);
ylabel(['Principal Component 2 (' strjoin(pc2_labels, ', ') ')']);
title('PCA - Ranking 4');
grid on;
saveas(gcf, 'PCA_Ranking_4_Scatter.png'); % Save scatter plot

figure;
pareto(explained4);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot - Ranking 4');
saveas(gcf, 'PCA_Ranking_4_Scree.png'); % Save scree plot
