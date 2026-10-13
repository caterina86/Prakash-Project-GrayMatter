
% Caterina Pedersini, April 2024
% cpedersini86@gmail.com

% LME at vertex-level, to assess the effect of surgery on each structural
% measure, in the group of congenital cataract-reversal individuals

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
smoothing = 5;
longitudinal = 0;
surface = [5 6];

outputDir = fullfile(projectDir, 'structural_measures/LME_WholeBrain_AllPatients', [num2str(smoothing) 'FWHM']);
mkdir(outputDir)

% define parameters:
measure = {'t1Intensity_0.5', 'thickness', 'curv'};
contrast = 'slope';
hemi= {'lh' 'rh'};
first_patient_idx = 173;



%% Define parameters and load data:

load(fullfile(resultsDir, ['residuals_zintra_int_' num2str(smoothing) 'FWHM.mat']));
load(fullfile(resultsDir, ['CT_zscore_' num2str(smoothing) 'FWHM.mat']));
load(fullfile(resultsDir,['Curv_zscore_' num2str(smoothing) 'FWHM.mat']));
% select only patients
ct = z_intra_ct(first_patient_idx:end);
curv = z_intra_curv(first_patient_idx:end);
t1 = residuals_zintra_int_T1_intensity(first_patient_idx:end);


Patients = readtable(fullfile(baseDir, 'derivatives/subjects/Patients.xlsx'), 'Sheet', 1);
Patients(56:end,:)=[];
ind = find(Patients.DaysSinceSurgery < 0);
Patients.DaysSinceSurgery(ind) = 1;


%% Run the LME at vertex-level

for measure_i = 1:length(measure)

    disp(measure{measure_i})
    numVertices = size(t1{1,1}, 1);
    numSubjects = size(Patients,1);

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
                T1(sub_i) = t1{sub_i}(vertex_i, 1);
                CT(sub_i) = ct{sub_i}(vertex_i, 1);
                Curv(sub_i) = curv{sub_i}(vertex_i, 1);
            end

            % Check for NaN values
            if sum(isnan(T1)) > 1 || sum(isnan(CT)) > 1 || sum(isnan(Curv)) > 1
                lme{vertex_i} = [];
            else    

                % Create the input table and select only subjects with longitudinal data
                Input = table(Patients.Subjects, T1, CT, Curv,log10(Patients.DaysSinceSurgery), Patients.AgeAtMeasurementYears, Patients.Gender, 'VariableNames', {'Subjects', 'T1', 'CT', 'Curv','TimeSurgery', 'AgeAtMeasurement','Gender'});        
                
                if longitudinal == 1
                    Input.Subjects = categorical(Input.Subjects);
                    % Count occurrences of each subject
                    subject_counts = tabulate(Input.Subjects);
                    % Extract subjects that occur more than once
                    repeated_subjects = subject_counts(:,1); % Subjects
                    counts = subject_counts(:,2); % Counts
                    repeated_subjects = repeated_subjects(cell2mat(counts) > 1);
                    Input = Input(ismember(Input.Subjects, repeated_subjects), :);
                end
                
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
                lme{vertex_i} = fitlme(Input, [response_variable ' ~ TimeSurgery * AgeAtMeasurement + Gender + (1|Subjects) ']);

            end

            if isempty(lme{vertex_i})
                pvals(vertex_i, 1) = NaN;
                betas(vertex_i, 1) = NaN;
                pvals(vertex_i, 2) = NaN;
                betas(vertex_i, 2) = NaN;
            else
                pvals(vertex_i, 1) = lme{vertex_i}.Coefficients.pValue(2); % time since surgery
                betas(vertex_i, 1) = lme{vertex_i}.Coefficients.Estimate(2); % time since surgery
                pvals(vertex_i, 2) = lme{vertex_i}.Coefficients.pValue(3); % age
                betas(vertex_i, 2) = lme{vertex_i}.Coefficients.Estimate(3); % age
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
