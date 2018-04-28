% This file produces the result for Breast Cancer ER- repurposing
% System requirement: Matlab 2016 or later, Microsoft Office 365 or later, 
% which can read Excel file in .xlsx format
% the see the paper Nguyen et al, 'DeCoST: A New Approach in Drug
% Repurposing from Control System Theory' for more mathematical description



clear all
clc

% import the gene information
[num, text, temp] = xlsread('GeneInfo.xlsx');
geneName = text(2:size(text, 1), 1); % gene list
Expression = num; % gene expression

% import and preprocess the gene-gene interaction (from a pathway or
% interaction databases
[num, text, temp] = xlsread('GeneInteraction.xlsx');
text(1, :) = [];
adj = zeros(length(geneName)); % this is the matrix of gene-gene interaction, 
                               % which is the important part of the
                               % mathematical model
for i = 1:size(text, 1)
    [~, index1] = ismember(text{i, 1}, geneName);
    [~, index2] = ismember(text{i, 2}, geneName);
    if strcmp(text{i, 3}, 'Activation') == 1
        adj(index1, index2) = 1;
    elseif strcmp(text{i, 3}, 'Inhibition') == 1
        adj(index1, index2) = -1;
    end
end

%adj = adj';
for i = 1:size(adj, 1)
    sumCol = sum(adj(:, i));
    if sumCol > 0
        adj(:, i) = adj(:, i) / sumCol;
    end
end

% import the drug information
[num, text, temp] = xlsread('DrugInfo.xlsx');
drugName = text(2:size(text, 1), 1);
drugCategory = text(2:size(text, 1), 2);

% load the drug-vector interaction
drugVector = zeros(length(geneName), length(drugName)); % this matrix has the drug-gene interactions
[num, text, temp] = xlsread('DrugGeneInteraction.xlsx');
text(1, :) = [];
for i = 1:size(text, 1)
    [~, index1] = ismember(text{i, 2}, geneName);
    [~, index2] = ismember(text{i, 1}, drugName);
    if strcmp(text{i, 3}, 'Activation') == 1
        drugVector(index1, index2) = 1;
    elseif strcmp(text{i, 3}, 'Inhibition') == 1
        drugVector(index1, index2) = -1;
    end
end

% begin computing drug therapeutic score
B = eye(size(adj, 1)); % identity matrix B means Bu = u. See equation 2
Q = eye(size(adj, 1)); % Q and R matrix are needed in DARE algorithm. Setting them to identity matrix would lead to equation 7 and 8
R = eye(size(adj, 1));
[P,L,G, report] = dare(adj, B, Q, R);% P is the result of DARE algorithm, equation 9 
F = ( ( R + B'*P*B ) ^(-1) ) * B' * P * adj;
u= -F*Expression;

numData = 100;
x_t = ones(length(Expression), numData);
u_t = ones(length(Expression), numData);
x_t(:, 1) = Expression;
u_t(:, 1) = u;
for t = 2:numData
    x_t(:, t) = adj*x_t(:, t-1) + u_t(:, t-1);
    u_t(:, t) = -F*x_t(:, t);
end

%u(find(find(u>-0.5&u<0.49))) = 0;
u = sign(u); % u is the hypo-control

drugRankingScore = zeros(length(drugName), 1); % compute Td score
for i = 1:length(drugName)
    numOverlap = dot( abs(u), abs(sign(drugVector(:, i))) );
    if numOverlap > 0
        drugRankingScore(i) = dot(u, sign(drugVector(:, i)) ) / numOverlap;
    end
end


%% draw the result
close all;
zeroScore = find(drugRankingScore == 0);
% boxplot of the Td score
greenID = find ( strcmp(drugCategory, 'D1') == 1); % green for D1 (approved drugs) group
greenID = setdiff(greenID, zeroScore);
redID = find( strcmp(drugCategory, 'D2') == 1);% green for D2 (rejected/withdrawn drugs) group
redID = setdiff(redID, zeroScore);
yellowID = find( strcmp(drugCategory, 'D3') == 1);% yellow for D3 (candidates of repurposing) group
yellowID = setdiff(yellowID, zeroScore);

redGroup = drugRankingScore(redID);
yellowGroup = drugRankingScore(yellowID);
greenGroup = drugRankingScore(greenID);

hold on;
for i = 1:length(greenGroup)
    plot(1, greenGroup(i), 'g.', 'MarkerSize',16);
end
plot([0.8 1.2], [median(greenGroup) median(greenGroup)], '-g', 'LineWidth', 4);

for i = 1:length(redGroup)
    plot(2, redGroup(i), 'r.', 'MarkerSize',16);
end
plot([1.8 2.2], [median(redGroup) median(redGroup)], '-r', 'LineWidth', 4);

for i = 1:length(yellowGroup)
    plot(3, yellowGroup(i), 'y.', 'MarkerSize',16);
end
plot([2.8 3.2], [median(yellowGroup) median(yellowGroup)], '-y', 'LineWidth', 4);

set(gca,'XTick',[1:1:3]);
labelList = {'D1'; 'D2' ; 'D3'};
set(gca,'XTickLabel',labelList);


axis([0 4 -1 1]);
hold off;

% draw the AUC
label = [ones(length(greenGroup), 1); -ones(length(redGroup), 1)];
score = [greenGroup; redGroup];
[X,Y,T,AUC] = perfcurve(label,score,1);
figure, plot(X,Y)
legend(['AUC: ', num2str(AUC)]);
xlabel('False positive rate');
ylabel('True positive rate');
title('ROC for classification between D1 and D2 group');
