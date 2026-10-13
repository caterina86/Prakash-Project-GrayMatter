
% Caterina Pedersini, April 2024
% cpedersini86@gmail.com

% Script to extract the residuals of the T1-w signal 

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

subjects = (importdata(fullfile(projectDir, 'subjects.txt')));
smoothing = [0 5 10];
hemi = {'lh' 'rh'};
indx_samples = 172; % last index of controls


%% Extract the residuals of T1-w signal using z-scores

for smooth_i = 1:length(smoothing)
    
    load(fullfile(resultsDir, ['T1_intensities_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    load(fullfile(resultsDir, ['CT_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    load(fullfile(resultsDir, ['Curv_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    t1 = z_intra_t1;
    cortical_thickness = z_intra_ct;
    curvature = z_intra_curv;  

    for sub_i = 1:length(subjects)

        T1_intensity = t1{sub_i}(:,1);
        CT = cortical_thickness{sub_i};       
        CV = curvature{sub_i};

        % Generate additional features for quadratic terms
        cortical_thickness_sq = CT.^2;
        curvature_sq = CV.^2;

        T = table(T1_intensity, CT, CV, cortical_thickness_sq, curvature_sq, 'VariableNames', {'T1', 'CT', 'CV', 'CT2', 'CV2'});

        lm = fitlme(T, 'T1 ~ CT + CV + CT*CV + CT2 + CV2');

        % Extract residuals
        residuals_T1_intensity{sub_i}(:,1) = lm.Residuals.Raw;      
        residuals_int_T1_intensity{sub_i}(:,1) = lm.Residuals.Raw + lm.Coefficients.Estimate(1);

    end

    save(fullfile(resultsDir, ['residuals_zintra_int_' num2str(smoothing(smooth_i)) 'FWHM.mat']), 'residuals_int_T1_intensity', '-v7.3');



    %% Create the map of the midsurface

    res = 0; % 1 -> only residuals; 0 -> residuals + intercept

    for sub_i = 1:length(subjects) 
        if res == 1 % only residuals
            residuals(:,sub_i) = residuals_T1_intensity{sub_i}(:,1);
        else % residuals + intercept
            residuals(:,sub_i) = residuals_int_T1_intensity{sub_i}(:,1);
        end
    end
    mean_t1_intensity_patients = mean(residuals(:,(indx_samples+1):end),2);
    mean_t1_intensity_controls = mean(residuals(:,1:indx_samples),2);

end



%% Extract the residuals of T1-w signal using raw data

for smooth_i = 1:length(smoothing)
    
    load(fullfile(resultsDir, ['T1_intensities_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    load(fullfile(resultsDir, ['CorticalThickness_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    load(fullfile(resultsDir, ['Curvature_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    t1 = t1intensity;
    for sub_i = 1:length(subjects)
        cortical_thickness{1,sub_i} = [ct{sub_i}(:,1);ct{sub_i}(:,2)];
        curvature{1,sub_i} = [curv{sub_i}(:,1);curv{sub_i}(:,2)];
    end
   
    for sub_i = 1:length(subjects)

        T1_intensity = t1{sub_i}(:,1);
        CT = cortical_thickness{sub_i};       
        CV = curvature{sub_i};

        % Generate additional features for quadratic terms
        cortical_thickness_sq = CT.^2;
        curvature_sq = CV.^2;

        T = table(T1_intensity, CT, CV, cortical_thickness_sq, curvature_sq, 'VariableNames', {'T1', 'CT', 'CV', 'CT2', 'CV2'});

        lm = fitlme(T, 'T1 ~ CT + CV + CT*CV + CT2 + CV2');

        % Extract residuals
        residuals_T1_intensity{sub_i}(:,1) = lm.Residuals.Raw;      
        residuals_int_T1_intensity{sub_i}(:,1) = lm.Residuals.Raw + lm.Coefficients.Estimate(1);

    end
    
    save(fullfile(resultsDir, ['residuals_int_' num2str(smoothing(smooth_i)) 'FWHM.mat']), 'residuals_int_T1_intensity', '-v7.3');

    res = 0; % 1 -> only residuals; 0 -> residuals + intercept

    for sub_i = 1:length(subjects) 
        if res == 1 % only residuals
            residuals(:,sub_i) = residuals_T1_intensity{sub_i}(:,1);
        else % residuals + intercept
            residuals(:,sub_i) = residuals_int_T1_intensity{sub_i}(:,1);
        end
    end
    mean_t1_intensity_patients = mean(residuals(:,(indx_samples+1):end),2);
    mean_t1_intensity_controls = mean(residuals(:,1:indx_samples),2);

end

