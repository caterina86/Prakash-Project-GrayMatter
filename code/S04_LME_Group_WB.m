% Caterina Pedersini, April 2024
% cpedersini86@gmail.com

% LME at vertex-level, to assess the effect of group on each 
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
resultsDir = fullfile(projectDir, 'structural_measures/WholeBrain');
smoothing = 0;
outputDir = fullfile(projectDir, 'structural_measures/LME_WholeBrain_Allsubjects', [num2str(smoothing) 'FWHM']);
mkdir(outputDir)

setenv('SUBJECTS_DIR', projectDir); 

subjects = importdata(fullfile(projectDir, 'subjects_lme.txt'));
hemi = {'lh' 'rh'};
subject_all = 1;
indx_samples = 172;
measure = {'t1Intensity_0.5', 'thickness', 'curv'};
contrast = 'CB-HC';
surface = [5 6];


%% Load data 
load(fullfile(resultsDir, ['residuals_zintra_int_' num2str(smoothing) 'FWHM.mat']));
load(fullfile(resultsDir, ['CT_zscore_' num2str(smoothing) 'FWHM.mat']));
load(fullfile(resultsDir, ['Curv_zscore_' num2str(smoothing) 'FWHM.mat']));
  
Controls = readtable(fullfile(baseDir, 'derivatives/subjects/Controls.xlsx'));
Patients = readtable(fullfile(baseDir, 'derivatives/subjects/Patients.xlsx'), 'Sheet', 1);
Patients(56:end,:)=[];

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

    ROI_T1 = ROI_T1(ind_all==1,:);
    ROI_CT = ROI_CT(ind_all==1,:);
    ROI_Curv = ROI_Curv(ind_all==1,:);

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

    ROI_T1 = ROI_T1(ind_all==1,:);
    ROI_CT = ROI_CT(ind_all==1,:);
    ROI_Curv = ROI_Curv(ind_all==1,:);

    % 1 -> controls; 0 -> patients
    group = [zeros(indx_samples,1); ones(sum(ind),1)];
    group_boxplot = [ones(indx_samples,1); zeros(sum(ind),1)]; % 1 -> controls; 0 -> patients
    Patients = Patients(Patients.DaysSinceSurgery > x, :);
    age = [Controls.Age; Patients.AgeAtMeasurementYears];
    gender = [Controls.Gender; Patients.Gender];
end



%% Run the LME at vertex-level

for measure_i = 1:length(measure)
    
    disp(measure{measure_i})

    numVertices = size(residuals_int_T1_intensity{1,1}, 1);
    numSubjects = length(subjects);

    % Preallocate memory
    T1 = nan(numSubjects, 1);
    CT = nan(numSubjects, 1);
    Curv = nan(numSubjects, 1);
    lme = cell(numVertices, 1);
    pvals = nan(numVertices, length(measure));
    betas = nan(numVertices, length(measure));

    file_betas = fullfile(outputDir, ['betas_' measure{measure_i} '_' num2str(smoothing) 'FWHM.mat']);
    file_pvals = fullfile(outputDir, ['pvals_' measure{measure_i} '_' num2str(smoothing) 'FWHM.mat']);

    if ~exist(file_betas, 'file') && ~exist(file_pvals, 'file')

        for vertex_i = 1:numVertices

            % Extract data for all subjects
            for sub_i = 1:numSubjects
                T1(sub_i) = residuals_int_T1_intensity{sub_i}(vertex_i, 1);
                CT(sub_i) = z_intra_ct{sub_i}(vertex_i, 1);
                Curv(sub_i) = z_intra_curv{sub_i}(vertex_i, 1);
            end

            % Check for NaN values
            if sum(isnan(T1)) > 1 || sum(isnan(CT)) > 1 || sum(isnan(Curv)) > 1
                lme{vertex_i} = [];
            else
                % Create table
                T = table(subjects, T1, CT, Curv, group, group_boxplot, age, gender, ...
                          'VariableNames', {'Subjects', 'T1', 'CT', 'Curv', 'Group', 'Group_Boxplot', 'AgeAtMeasurement', 'Gender'});

                % Select response variable
                switch measure_i
                    case 1
                        response_variable = 'T1';
                    case 2
                        response_variable = 'CT';
                    case 3
                        response_variable = 'Curv';
                end

                % Fit linear mixed-effects model
                lme{vertex_i} = fitlme(T, [response_variable ' ~ Group + AgeAtMeasurement + Gender + (1|Subjects) ']);

                % extract the random effect
                % randomEffectsTable = randomEffects(lme{vertex_i});

            end

            if isempty(lme{vertex_i})
                pvals(vertex_i, 1) = NaN;
                betas(vertex_i, 1) = NaN;
            else
                pvals(vertex_i, 1) = lme{vertex_i}.Coefficients.pValue(2);
                betas(vertex_i, 1) = lme{vertex_i}.Coefficients.Estimate(2);
            end
        end

        save(file_betas, 'betas', '-v7.3');
        save(file_pvals, 'pvals', '-v7.3');

    else

        disp('LME already performed');
        load(file_betas)
        load(file_pvals)

    end
end

close all;
