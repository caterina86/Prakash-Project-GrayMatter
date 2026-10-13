% Caterina Pedersini, April 2024
% cpedersini86@gmail.com

% LME at ROI-level, to assess the effect of surgery on each structural
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
outputDir = fullfile(projectDir, 'structural_measures/DesikanKilliany');
setenv('SUBJECTS_DIR', projectDir); 

smoothing = 0;
hemi = {'lh' 'rh'};
indx_samples = 172; % last index of controls
output_plot = 0;
effect_i = 2; % 2 -> time since surgery; 3 -> age at measurement


%% Load data and define paramethers

load(fullfile(outputDir, ['ct_parcels_' num2str(smoothing) 'FWHM']));
load(fullfile(outputDir, ['curv_parcels_' num2str(smoothing) 'FWHM']));
load(fullfile(outputDir, ['t1_parcels_' num2str(smoothing) 'FWHM']));

subjects = importdata(fullfile(projectDir, 'subjects_lme.txt'));        
Controls = readtable(fullfile(baseDir, 'derivatives/subjects/Controls.xlsx'));
Patients = readtable(fullfile(baseDir, 'derivatives/subjects/Patients.xlsx'), 'Sheet', 1);
Patients(56:end,:)=[];
ind = find(Patients.DaysSinceSurgery < 0);
Patients.DaysSinceSurgery(ind) = 1;

labels = {'bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal', 'cuneus', ...
    'entorhinal', 'frontalpole', 'fusiform', 'inferiorparietal', 'inferiortemporal', 'insula', ...
    'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal', ...
    'lingual', 'medialorbitofrontal', 'middletemporal', 'paracentral', 'parahippocampal', ...
    'parsopercularis', 'parsorbitalis', 'parstriangularis', 'pericalcarine', 'postcentral', ...
    'posteriorcingulate', 'precentral', 'precuneus', 'rostralanteriorcingulate', 'rostralmiddlefrontal', ...
    'superiorfrontal', 'superiorparietal', 'superiortemporal', 'supramarginal', 'temporalpole', 'transversetemporal'};


% Create the resultsDir
resultsDir = fullfile(outputDir, 'figures', ['smoothing_' num2str(smoothing) 'FWHM'], 'TimeSurgery_zscores');   
mkdir(resultsDir)



%% Run the LME to assess the fixed-effect of time since surgery
    
for label_i = 1:length(labels)

    T1 = table2array(t1_parcels(indx_samples+1:end,label_i));
    CT = table2array(ct_parcels(indx_samples+1:end,label_i));
    Curv = table2array(curv_parcels(indx_samples+1:end,label_i));
    P = table(Patients.Subjects, T1, CT, Curv,log10(Patients.DaysSinceSurgery), Patients.AgeAtMeasurementYears, Patients.AgeAtSurgeryYears, Patients.Gender, 'VariableNames', {'Subjects', 'T1', 'CT', 'Curv','TimeSurgery', 'AgeAtMeasurement','AgeAtSurgery', 'Gender'});


    % Define response variables
    response_variables = {'T1', 'CT', 'Curv'};

    % Create figure
    figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);

    % Loop over response variables
    for i = 1:numel(response_variables)

        dep_var = P.(response_variables{i});

        if sum(dep_var) == 0
        else
            lme1{i,label_i} = fitlme(P, [response_variables{i} ' ~ AgeAtMeasurement + Gender + (1|Subjects)']); % fit the null model (without the fixed-effect)
            lme2{i,label_i} = fitlme(P, [response_variables{i} ' ~ TimeSurgery * AgeAtMeasurement + Gender + (1|Subjects)']); % fit the alternative model
            bic1{i,label_i} = lme1{i,label_i}.ModelCriterion.BIC;
            bic2{i,label_i} = lme2{i,label_i}.ModelCriterion.BIC;

            logB12 = 0.5 * (bic1{i,label_i} - bic2{i,label_i});
            BF{i,label_i} = exp(logB12);

            if BF{i,label_i} >= 3
                fprintf('Strong evidence in favor of lme2 for %s %s\n', response_variables{i}, labels{label_i});
            elseif BF{i,label_i} <= 0.33
                fprintf('Evidence against the lme2 for %s %s\n', response_variables{i}, labels{label_i});
            elseif BF{i,label_i} > 0.33 && BF{i,label_i} < 3
                fprintf('No evidence for %s %s\n', response_variables{i}, labels{label_i});
            end
            
            if output_plot == 1
                subplot(1,3,i);
                h = gscatter(P.TimeSurgery, dep_var, P.Subjects, [], [], [], 15); % 15 is the new marker size
                hold on;

                % Iterate over each subject
                uniqueSubjects = unique(P.Subjects);
                for iii = 1:numel(uniqueSubjects)
                    subject = uniqueSubjects(iii);
                    % Extract data points for the current subject
                    subjectDataIdx = strcmp(cellstr(P.Subjects), cellstr(subject));
                    subjectTime = P.TimeSurgery(subjectDataIdx);
                    subjectT1 = dep_var(subjectDataIdx);

                    % Find the corresponding color from gscatter (h is an array of line objects)
                    subjectColor = h(iii).Color; 

                    % Plot the lines with the same color as the dots, increase line width
                    plot(subjectTime, subjectT1, 'Color', subjectColor, 'LineWidth', 2); % Set LineWidth to 2
                end

                % Extract the beta value (slope) from your linear mixed model
                beta = lme2{i,label_i}.Coefficients.Estimate(effect_i); % Replace with correct indexing

                % Create a range of x values for plotting the regression line (adjust as needed)
                xRange = linspace(min(P.TimeSurgery), max(P.TimeSurgery), 100); % 100 points between min and max

                % Calculate the y-intercept (optional, assume 0 if not applicable)
                intercept = lme2{i,label_i}.Coefficients.Estimate(1); % If intercept is known, replace 0 with the actual value
                % Calculate the corresponding y values for the regression line
                yFit = intercept + beta * xRange;

                % Plot the dashed line (beta line)
                plot(xRange, yFit, '--k', 'LineWidth', 2); % '--k' specifies a black dashed line

                % Customize the plot labels and increase font size
                xlabel('Log Time Since Surgery', 'FontSize', 18);
                ylabel(response_variables{i}, 'FontSize', 18);

                % Set axis tick label font size
                set(gca, 'FontSize', 16); % Set FontSize for axis numbers (ticks)

                % Ensure everything displays properly
                hold off;
                % title(labels{label_i});
                legend('off'); % Adjust the legend location as needed
                hold off; % Release the plot hold
                saveas(gca, fullfile(resultsDir, ['ROI_' labels{label_i} '_EffectSurgery.bmp']));
            end
        end
    end

end
    
% Correct p-values based on multiple comparisons
for ii = 1:length(labels)
    for measure = 1:3
        if isempty(lme2{measure,ii})
            betas(measure,ii) = NaN;
            SE(measure,ii) = NaN;
            pvals(measure,ii) = NaN;
        else
            betas(measure,ii) = lme2{measure,ii}.Coefficients.Estimate(effect_i);
            SE(measure,ii) = lme2{measure,ii}.Coefficients.SE(effect_i);
            pvals(measure,ii) = lme2{measure,ii}.Coefficients.pValue(effect_i);
        end   
    end
end    

alpha = 0.01;
for measure = 1:3
    % apply FDR correction
    [~,~,~,fdr_p] = fdr_bh(reshape(pvals(measure,:),[],1), alpha, 'pdep',true);
    fdr_p_all(measure,:) = reshape(fdr_p,1,[]);
end

min(min(fdr_p_all))

close all;


