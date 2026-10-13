% Caterina Pedersini, April 2024
% cpedersini86@gmail.com

% LME at ROI-level, to assess the effect of group on each 
% structural measure


clearvars;
close all;
clc

% load dependencies
githubDir = uigetdir(pwd, 'Select the base directory'); % define here the baseline Directory
addpath(genpath(fullfile(githubDir, 'fdr_bh'))); % https://github.com/cvnlab/cvncode
addpath(genpath(fullfile(githubDir, 'cvncode'))); % https://github.com/cvnlab/cvncode
addpath(genpath(fullfile(githubDir, 'knkutils'))); % https://github.com/cvnlab/knkutils
addpath(genpath(fullfile(githubDir, 'gifti'))); % https://github.com/gllmflndn/gifti

if isfolder('/Applications/freesurfer/7.4.1/bin')
    freesurferDir = '/Applications/freesurfer/7.4.1/';
    setenv('FREESURFER_HOME', '/Applications/freesurfer/7.4.1');
    PATH = getenv('PATH'); setenv('PATH', ['/opt/anaconda3/bin:/usr/local/bin:/usr/local/fsl/bin:/Applications/freesurfer/7.4.1/bin:' PATH]);
    addpath(genpath(fullfile(freesurferDir, 'matlab')));
elseif isfolder('/Applications/freesurfer/7.3.2/bin')
    freesurferDir = '/Applications/freesurfer/7.3.2/';
    setenv('FREESURFER_HOME', '/Applications/freesurfer/7.3.2');
    PATH = getenv('PATH'); setenv('PATH', ['/opt/anaconda3/bin:/usr/local/bin:/usr/local/fsl/bin:/Applications/freesurfer/7.3.2/bin:' PATH]);
    addpath(genpath(fullfile(freesurferDir, 'matlab')));
elseif isfolder('/Applications/freesurfer/7.2.0/bin')
    freesurferDir = '/Applications/freesurfer/7.2.0/';
    setenv('FREESURFER_HOME', '/Applications/freesurfer/7.2.0');
    PATH = getenv('PATH'); setenv('PATH', ['/opt/anaconda3/bin:/usr/local/bin:/usr/local/fsl/bin:/Applications/freesurfer/7.2.0/bin:' PATH]);
    addpath(genpath(fullfile(freesurferDir, 'matlab')));
elseif isfolder('/Applications/freesurfer/bin')
    freesurferDir = '/Applications/freesurfer';
    setenv('FREESURFER_HOME', '/Applications/freesurfer');
    PATH = getenv('PATH'); setenv('PATH', ['/opt/anaconda3/bin:/usr/local/bin:/usr/local/fsl/bin:/Applications/freesurfer/bin:' PATH]);
    addpath(genpath(fullfile(freesurferDir, 'matlab')));
else
    error('Cannot find freesurfer binary in or near /Applications/freesurfer')
end


%% load paths

baseDir = uigetdir(pwd, 'Select the base directory'); % define here the baseline Directory
projectDir = fullfile(baseDir, 'derivatives/freesurfer');    
outputDir = fullfile(projectDir, 'structural_measures/DesikanKilliany');

setenv('SUBJECTS_DIR', projectDir); 

step_surface = 0.5; 
smoothing = 0;
hemi = {'lh' 'rh'};
indx_samples = 172; % last index of controls
condition = 2; % 1 -> raw data; 2 -> z-scores separately for left and right hemisphere
output_plot = 1;



%% Load data and define paramethers

load(fullfile(outputDir, ['ct_parcels_' num2str(smoothing) 'FWHM']));
load(fullfile(outputDir, ['curv_parcels_' num2str(smoothing) 'FWHM']));
load(fullfile(outputDir, ['t1_parcels_' num2str(smoothing) 'FWHM']));

subjects = importdata(fullfile(projectDir, 'subjects_lme.txt'));        
Controls = readtable(fullfile(baseDir, 'derivatives/subjects/Controls.xlsx'));
Patients = readtable(fullfile(baseDir, 'derivatives/subjects/Patients.xlsx'), 'Sheet', 1);
Patients(56:end,:)=[];

subject_all = 1; % 1 -> subjects all; 2 -> only pre-surgery; 3 -> only post-surgery
if subject_all ~= 1
    x = 10; % number of days as a limit to include subjects in the subjects_all = 0 
end

measure = {'t1Intensity_0.5', 'thickness', 'curv'};
fixed_effect = 'Group';

% define the labels
labels = {'bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal', 'cuneus', ...
    'entorhinal', 'frontalpole', 'fusiform', 'inferiorparietal', 'inferiortemporal', 'insula', ...
    'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal', ...
    'lingual', 'medialorbitofrontal', 'middletemporal', 'paracentral', 'parahippocampal', ...
    'parsopercularis', 'parsorbitalis', 'parstriangularis', 'pericalcarine', 'postcentral', ...
    'posteriorcingulate', 'precentral', 'precuneus', 'rostralanteriorcingulate', 'rostralmiddlefrontal', ...
    'superiorfrontal', 'superiorparietal', 'superiortemporal', 'supramarginal', 'temporalpole', 'transversetemporal'};


if subject_all == 1
    % 1 -> controls; 0 -> patients
    group = [zeros(indx_samples,1); ones(length(subjects)-indx_samples,1)];
    group_boxplot = [ones(indx_samples,1); zeros(length(subjects)-indx_samples,1)]; % 1 -> controls; 0 -> patients
    age = [Controls.Age; Patients.AgeAtMeasurementYears];
    gender = [Controls.Gender; Patients.Gender];
    
elseif subject_all == 2 % timepoints within 10 days from the surgery
    ind = Patients.DaysSinceSurgery<=x;
    ind_all = [ones(172,1); ind];
    subjects = subjects(ind_all==1);

    % 1 -> controls; 0 -> patients
    group = [zeros(indx_samples,1); ones(sum(ind),1)];
    group_boxplot = [ones(indx_samples,1); zeros(sum(ind),1)]; % 1 -> controls; 0 -> patients
    Patients = Patients(Patients.DaysSinceSurgery <= x, :);
    age = [Controls.Age; Patients.AgeAtMeasurementYears];
    gender = [Controls.Gender; Patients.Gender];
    
elseif subject_all == 3 % timepoints > 10 days from surgery
    ind = Patients.DaysSinceSurgery>x;
    ind_all = [ones(172,1); ind];
    subjects = subjects(ind_all==1);

    % 1 -> controls; 0 -> patients
    group = [zeros(indx_samples,1); ones(sum(ind),1)];
    group_boxplot = [ones(indx_samples,1); zeros(sum(ind),1)]; % 1 -> controls; 0 -> patients
    Patients = Patients(Patients.DaysSinceSurgery > x, :);
    age = [Controls.Age; Patients.AgeAtMeasurementYears];
    gender = [Controls.Gender; Patients.Gender];
end


%% Run the LME at ROI-level to assess the effect of group

% Create the resultsDir
if subject_all == 1 % all subjects
    if condition == 1 
        resultsDir = fullfile(outputDir, 'figures', ['smoothing_' num2str(smoothing) 'FWHM'], 'Group_rawdata');   
    else
        resultsDir = fullfile(outputDir, 'figures', ['smoothing_' num2str(smoothing) 'FWHM'], 'Group_zscores');   
    end
elseif subject_all == 2
    if condition == 1 % 1 -> choose the z-values; 0 -> choose the original values
        resultsDir = fullfile(outputDir, 'figures', ['smoothing_' num2str(smoothing) 'FWHM'], [num2str(x) '_Group_rawdata'], 'presurgery');   
    else
        resultsDir = fullfile(outputDir, 'figures', ['smoothing_' num2str(smoothing) 'FWHM'], [num2str(x) '_Group_zscores'], 'presurgery');   
    end
elseif subject_all == 3
    if condition == 1 % 1 -> choose the z-values; 0 -> choose the original values
        resultsDir = fullfile(outputDir, 'figures', ['smoothing_' num2str(smoothing) 'FWHM'], [num2str(x) '_Group_rawdata'], 'postsurgery');   
    else
        resultsDir = fullfile(outputDir, 'figures', ['smoothing_' num2str(smoothing) 'FWHM'], [num2str(x) '_Group_zscores'], 'postsurgery');   
    end
end
mkdir(resultsDir)


for measure_i = 1:length(measure)
    
    for label_i = 1:length(labels)

        T1 = table2array(t1_parcels(:,label_i));
        CT = table2array(ct_parcels(:,label_i));
        Curv = table2array(curv_parcels(:,label_i));

        if subject_all ~= 1
            T1 = T1(ind_all==1,:);
            CT = CT(ind_all==1,:);
            Curv = Curv(ind_all==1,:);
        end            
        T = table(subjects, T1, CT, Curv, group, group_boxplot, age, gender, 'VariableNames', {'Subjects', 'T1', 'CT', 'Curv', 'Group', 'Group_Boxplot', 'AgeAtMeasurement', 'Gender'});

    
        %% Run the LME

        if measure_i == 1 % T1 
            response_variable = 'T1';
        elseif measure_i == 2 % CT
            response_variable = 'CT';
        elseif measure_i == 3 % CT
            response_variable = 'Curv';
        end

        lme{label_i} = fitlme(T, [response_variable '  ~ Group + AgeAtMeasurement + Gender + (1|Subjects)']);
        anova(lme{label_i})

        if output_plot == 1
            % Boxchart - effect of group
            % subplot(1,3,1)
            figure('Units', 'inches', 'Position', [1, 1, 4, 6]); % [left, bottom, width, height]
            boxchart(T.(response_variable),'GroupByColor',T.Group_Boxplot)
            legend('Patients', 'Controls', 'Location', 'southoutside');
            ylabel(response_variable, 'FontSize', 18); %ylim([-1 1]);
            xlabel('Group', 'FontSize', 18);
            title([labels{label_i}], 'FontSize', 18)
            saveas(gca, fullfile(resultsDir, [response_variable '_Cluster_' labels{label_i} '_EffectGroup.bmp']));

            % Scatterplot - effect of age
            group1_age = T.AgeAtMeasurement(T.Group_Boxplot == 0);
            group1_var = T.(response_variable)(T.Group_Boxplot == 0);
            group2_age = T.AgeAtMeasurement(T.Group_Boxplot == 1);
            group2_var = T.(response_variable)(T.Group_Boxplot == 1);
            
            figure;
            scatter(group1_age, group1_var, 'b', 'filled'); % Group 1
            hold on;
            scatter(group2_age, group2_var, 'r', 'filled'); % Group 2
            hold on;
            % Extract fixed effects coefficients from lme
            coefficients = lme{label_i}.Coefficients.Estimate;
            age_index = find(strcmp(lme{label_i}.CoefficientNames, 'AgeAtMeasurement'));
            intercept_index = find(strcmp(lme{label_i}.CoefficientNames, '(Intercept)'));
            % Extract slope and intercept
            slope = coefficients(age_index);
            intercept = coefficients(intercept_index) + coefficients(age_index); % Add fixed effect of Group
            % Plot line with slope and intercept
            x_values = 0:25; % Use min and max age from both groups
            y_values = slope * x_values + intercept;
            plot(x_values, y_values, 'k-', 'LineWidth', 2);
            hold off;
            xlabel('Age', 'FontSize', 18);
            ylabel(response_variable, 'FontSize', 18);%ylim([-1 1]);
            legend('Patients', 'Controls', 'Location', 'Best');
            title([labels{label_i}], 'FontSize', 18)
            saveas(gca, fullfile(resultsDir, [response_variable '_Cluster_' labels{label_i} '_Age.bmp']));

            % Boxchart - effect of gender
            %boxchart(T.(response_variable),'GroupByColor',T.Gender)
            %legend('Female', 'Male', 'Location', 'Best');
            %ylabel(response_variable, 'FontSize', 18); %ylim([-1 1]);
            %xlabel('Gender', 'FontSize', 18);
            %title([labels{label_i}], 'FontSize', 18)

        end
    end

    for ii = 1:length(labels)

        betas(1,ii) = lme{ii}.Coefficients.Estimate(2);
        betas(2,ii) = lme{ii}.Coefficients.Estimate(3);
        betas(3,ii) = lme{ii}.Coefficients.Estimate(4);

        SE(1,ii) = lme{ii}.Coefficients.SE(2);
        SE(2,ii) = lme{ii}.Coefficients.SE(3);
        SE(3,ii) = lme{ii}.Coefficients.SE(4);

        pvals(1,ii) = lme{ii}.Coefficients.pValue(2);
        pvals(2,ii) = lme{ii}.Coefficients.pValue(3);
        pvals(3,ii) = lme{ii}.Coefficients.pValue(4);

    end

    % Significant level (adjust as needed)
    alpha = 0.01;
    % apply FDR correction
    [~,~,~,fdr_p] = fdr_bh(reshape(pvals,[],1), alpha, 'pdep',true);
    fdr_p = reshape(fdr_p,3, []);

    ROI = 1:length(labels); % ROI values
    ROIs1 = ROI(1,(fdr_p(1, :) < alpha));
    Betas1 = betas(1,(fdr_p(1, :) < alpha));
    SE1 = SE(1,(fdr_p(1, :) < alpha));
    ROIs2 = ROI(1,(fdr_p(1, :) >= alpha));
    Betas2 = betas(1,(fdr_p(1, :) >= alpha));
    SE2 = SE(1,(fdr_p(1, :) >= alpha));

    % Plotting with error bars
    figure('Position', [100, 100, 1200, 600]); % Width = 1200, Height = 600 for a rectangular shape
    errorbar(ROIs1, Betas1, SE1, 'r', 'LineStyle', 'none', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 2, 'MarkerEdgeColor','r',...
                  'MarkerFaceColor','r');
    hold on;
    errorbar(ROIs2, Betas2, SE2, 'b', 'LineStyle', 'none', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 2, 'MarkerEdgeColor','b',...
                  'MarkerFaceColor','b');
    yline(0)
    ylabel([response_variable ' Betas'], 'FontSize', 18);
    ylim([-0.6 0.6])
    xlim([0 size(labels,2)+1])
    title(['Effect of ' fixed_effect], 'FontSize', 18)
    xticks(1:length(labels));
    xticklabels(labels);
    xtickangle(45);
    % Adjust the font size of the x-tick labels
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 20; % Set the font size as desired, e.g., 14

    % Add dashed lines for each label
    for i = 1:length(labels)
        xline(i, '--k', 'LineWidth', 0.5); % Dashed vertical line at each x-tick position
    end

    % Adjust the figure's position to ensure labels are clear
    set(gcf, 'PaperPositionMode', 'auto'); % This ensures that the figure is saved correctly
    exportgraphics(gcf, fullfile(resultsDir, ['Effect of' fixed_effect '_' response_variable '_ROIs.pdf']), 'ContentType', 'vector');

end

close all;
