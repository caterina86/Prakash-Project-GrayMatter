
% Caterina Pedersini, September 2024
% cpedersini86@gmail.com

% Plot the raw structural measures, to visually check for differences
% between groups

clearvars;
close all;
clc

% load dependencies
githubDir = uigetdir(pwd, 'Select the base directory'); % define here the baseline Directory
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
resultsDir = fullfile(projectDir, 'structural_measures/WholeBrain');
setenv('SUBJECTS_DIR', projectDir); 

subjects = importdata(fullfile(projectDir, 'subjects.txt'));
smoothing = [0 5 10];
hemi = {'lh' 'rh'};
indx_samples = 172;


%% Cortical thickness and curvature

% load data
for smooth_i = 1:length(smoothing)

    load(fullfile(resultsDir, ['CorticalThickness_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    load(fullfile(resultsDir, ['Curvature_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    load(fullfile(resultsDir, ['T1_intensities0.5_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    
    for sub_i = 1:length(subjects) 
        ct_all(:,sub_i) = [(ct{sub_i}(:,1)); (ct{sub_i}(:,2))];
        ct_all(ct_all==0) = NaN;
        cv_all(:,sub_i) = [(curv{sub_i}(:,1)); (curv{sub_i}(:,2))];
        t1_all(:,sub_i) = [(int{sub_i}(:,1)); (int{sub_i}(:,2))];
    end

    % plot the distribution
    group_data_ct = [mean(ct_all(:, 1:indx_samples), 2), mean(ct_all(:, indx_samples+1:end), 2)];
    group_data_cv = [mean(cv_all(:, 1:indx_samples), 2), mean(cv_all(:, indx_samples+1:end), 2)];
    group_data_t1 = [mean(t1_all(:, 1:indx_samples), 2), mean(t1_all(:, indx_samples+1:end), 2)];


    figure;
    % Plot for the first dataset (group_data_ct)
    subplot(1, 3, 1);
    violinplot(group_data_ct); % Include group labels for correct plotting
    hold on;
    % Customize the appearance of the violin plots
    set(gca, 'LineWidth', 2); % Thicker outline for the violin plots
    set(gca, 'XTick', []); % Removes the x-axis labels
    title({'Average', 'Cortical Thickness'}, 'FontSize', 18);

    subplot(1, 3, 2);
    violinplot(group_data_cv); % Include group labels for correct plotting
    hold on;
    % Customize the appearance of the violin plots
    set(gca, 'LineWidth', 2); % Thicker outline for the violin plots
    set(gca, 'XTick', []); % Removes the x-axis labels
    title({'Average', 'Curvature'}, 'FontSize', 18);

    subplot(1, 3, 3);
    violinplot(group_data_t1); % Include group labels for correct plotting
    hold on;
    % Customize the appearance of the violin plots
    set(gca, 'LineWidth', 2); % Thicker outline for the violin plots
    set(gca, 'XTick', []); % Removes the x-axis labels
    title({'Average', 'T1-w signal'}, 'FontSize', 18);
    exportgraphics(gcf, fullfile(resultsDir, 'figures', ['violinplot_rawData_' num2str(smoothing(smooth_i)) 'FWHM.pdf']), 'ContentType', 'vector')

end

close all;
