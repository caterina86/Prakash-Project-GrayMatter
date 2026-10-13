
% Caterina Pedersini, April 2024
% cpedersini86@gmail.com

% Script to calculate the intra-subject z-score for the T1 intensity, CT
% and Curvature, save them and plot them on the surface (save the .mgh file)

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
addpath(genpath(fullfile(baseDir, 'code'))); 
projectDir = fullfile(baseDir, 'derivatives/freesurfer');
resultsDir = fullfile(projectDir, 'structural_measures/WholeBrain');
mkdir(fullfile(resultsDir,'figures/'));
setenv('SUBJECTS_DIR', projectDir); 

subjects = importdata(fullfile(projectDir, 'subjects.txt'));
smoothing = [0 5 10];
hemi = {'lh' 'rh'};
indx_samples = 172;
condition = 0; % 1 -> whole brain z-scores; 0 -> z-scores separately for left and right hemisphere
surface_i = 0.5; % which surface to take -> midsurface
surface_plot = 1; % 1 -> create the .mgh file of the surface output; 0 -> NOT create the .mgh files



%% Convert the t1-intensity maps in intrasubject z-scores
 
for smooth_i = 1:length(smoothing)
    
    disp(['Smoothing ' num2str(smoothing(smooth_i)) ' FWHM']);

    if condition == 1 % if z-score extracted from the whole brain (left and right hemisphere together)

        file_t1 = fullfile(resultsDir, ['T1_intensities_zscoreWB_' num2str(smoothing(smooth_i)) 'FWHM.mat']);
        
        if exist(file_t1, 'file') ~= 2 
            
            load(fullfile(resultsDir, ['T1_intensities' num2str(surface_i) '_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
            
            for sub_i = 1:length(subjects)
                int_sub = reshape(int{sub_i},[],1);
                mean_sub = nanmean(int_sub,1);
                std_sub = nanstd(int_sub,1);

                z_intra_t1{sub_i}(:,1) = ((int_sub - mean_sub)./std_sub);
                t1intensity{sub_i}(:,1) = int_sub;
            end
            
            save(file_t1, 'z_intra_t1', '-v7.3');
            save(fullfile(resultsDir, ['T1_intensities_' num2str(smoothing(smooth_i)) 'FWHM.mat']), 't1intensity', '-v7.3');

        end
        
    else % if z-score calculated separately for the 2 hemispheres
        file_t1 = fullfile(resultsDir, ['T1_intensities_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']);
        
        if exist(file_t1, 'file') ~= 2 
            
            load(fullfile(resultsDir, ['T1_intensities' num2str(surface_i) '_' num2str(smoothing(smooth_i)) 'FWHM.mat']));

            for sub_i = 1:length(subjects)
                for hemi_i = 1:length(hemi)
                    int_sub = int{sub_i}(:,hemi_i);
                    mean_sub_t1(sub_i) = nanmean(int_sub,1);
                    std_sub_t1(sub_i) = nanstd(int_sub,1);

                    z_intra(:,hemi_i) = ((int_sub - mean_sub_t1(sub_i))./std_sub_t1(sub_i));
                end

                z_intra_t1{sub_i}(:,1) = [z_intra(:,1); z_intra(:,2)];
                t1intensity{sub_i}(:,1) = reshape(int{sub_i},[],1);
            end
            
            save(file_t1, 'z_intra_t1', '-v7.3');
            save(fullfile(resultsDir, ['T1_intensities_' num2str(smoothing(smooth_i)) 'FWHM.mat']), 't1intensity', '-v7.3');
        end
    end

    clear z_intra t1intensity z_intra_t1

end        




%% Save the output of the midsurface for visualization on the surface

roiLabel = readtable(fullfile(projectDir, 'fsaverage/atlasmgz/Glasser2016_ColorLUT.txt'));
roiIdx = get_roi('fsaverage'); % two columnn vector, one for left hemi and one for right hemi
max(roiIdx); % has to be 180 for Glasser atlas

% extract cortical thickness values for the ROIs where we observe a difference in the vertex-wise analysis
list = roiLabel.LabelName_(2:end);

patients = [0 1];
surface = [5 6];

for smooth_i = 1:length(smoothing)

    load(fullfile(resultsDir, ['T1_intensities_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']));

    % extract the t1-w intensity for all subjects
    for sub_i = 1:length(subjects) 
        t1_intensity_zscore(:,sub_i) = (z_intra_t1{sub_i}(:,1));
    end
    mean_t1_intensity_patients = mean(t1_intensity_zscore(:,(indx_samples+1):end),2);
    mean_t1_intensity_controls = mean(t1_intensity_zscore(:,1:indx_samples),2);

    % create a mask to elimiate all the values outside the Glasser Atlas
    hSize = get_surfsize('fsaverage'); % get fsnative hemi size
    hSizeIdx=[1,hSize(1);hSize(1)+1,(hSize(1) + hSize(2))];  % get hemi index where the maximum is the number of vertices with a functional activation (pRF)
    hh = {'lh','rh'}; % two hemispheres, left and right
    
    for iROI = 1:max(max(roiIdx))
        % create an index with the vertices that belong to a specific ROI
        indx = roiIdx(:,1)==iROI;
        roiIndices(:, iROI) = indx;
    end    
    
    % Create a vector with the same length as the number of vertices
    mergedVector_lh = zeros(size(roiIdx, 1), 1);
    mergedVector_lh(any(roiIndices, 2)) = 1;    
    
    clear roiIndices
    
    for iROI = 1:max(max(roiIdx))
        % create an index with the vertices that belong to a specific ROI
        indx = roiIdx(:,2)==iROI;
        roiIndices(:, iROI) = indx;
    end    
    
    % Create a vector with the same length as the number of vertices
    mergedVector_rh = zeros(size(roiIdx, 1), 1);
    % Set entries to 1 if any ROI index is true for that vertex
    mergedVector_rh(any(roiIndices, 2)) = 1;    
    
    mergedVector = [mergedVector_lh; mergedVector_rh];

    if surface_plot == 1

        for group_i = 1:length(patients)

            for surface_ii = 1:length(surface)

                if group_i == 1
                    vals = mean_t1_intensity_patients;
                else
                    vals = mean_t1_intensity_controls;
                end

                valsalpha = ~isnan(vals);
                vals(isnan(vals)) = 0;
                vals(mergedVector == 0) = nan;
                [~,~, rgbimg] = cvnlookup('fsaverage',surface(surface_ii),vals,[-2 2],(jet),[],[],[],{'rgbnan',0,'overlayalpha',valsalpha});
                a = imagesc(rgbimg); axis image tight; axis off;
                title('Mean T1 intensity')
                colormap((jet)) % color map
                caxis([-2 2]) % range of the colorbar    
                cb = colorbar('SouthOutside');
                ylabel(cb,'T1 intensity intra-subject z-scores','Rotation',360);
                set(gca,'FontSize',20) 

                if group_i == 1
                    if surface(surface_ii) == 6
                        exportgraphics(gca, fullfile(resultsDir,'figures/', ['Mean_T1Intensity_zscores_' num2str(surface_i) '_' num2str(smoothing(smooth_i)) 'FWHM_patients_lateral.pdf']));
                    else
                        exportgraphics(gca, fullfile(resultsDir,'figures/', ['Mean_T1Intensity_zscores_' num2str(surface_i) '_' num2str(smoothing(smooth_i)) 'FWHM_patients_medial.pdf']));
                    end
                else 
                    if surface(surface_ii) == 6
                        exportgraphics(gca, fullfile(resultsDir,'figures/', ['Mean_T1Intensity_zscores_' num2str(surface_i) '_' num2str(smoothing(smooth_i)) 'FWHM_controls_lateral.pdf']));
                    else
                        exportgraphics(gca, fullfile(resultsDir,'figures/', ['Mean_T1Intensity_zscores_' num2str(surface_i) '_' num2str(smoothing(smooth_i)) 'FWHM_controls_medial.pdf']));
                    end
                end
            end
        end
    end
        
    clear z_intra

end

close all;


%% Convert Cortical thickness and curvature in intrasubject z-scores

for smooth_i = 1:length(smoothing)
    
    disp(['Smoothing ' num2str(smoothing(smooth_i)) ' FWHM']);
    
    if condition == 1 % if z-score extracted from the whole brain (left and right hemisphere together)
        
        file_ct = fullfile(resultsDir, ['CT_zscoreWB_' num2str(smoothing(smooth_i)) 'FWHM.mat']);
        file_cv = fullfile(resultsDir, ['Curv_zscoreWB_' num2str(smoothing(smooth_i)) 'FWHM.mat']);
        
        if exist(file_ct, 'file') ~= 2 && exist(file_cv, 'file') ~= 2
        
            load(fullfile(resultsDir, ['CorticalThickness_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
            load(fullfile(resultsDir, ['Curvature_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    
            for sub_i = 1:length(subjects)
                ct_sub = reshape(ct{sub_i},[],1);
                ct_sub(ct_sub==0) = NaN;
                mean_sub_ct(sub_i) = nanmean(ct_sub,1);
                std_sub = nanstd(ct_sub,1);
                z_intra_ct{sub_i}(:,1) = ((ct_sub - mean_sub_ct(sub_i))./std_sub);

                curv_sub = reshape(curv{sub_i},[],1);
                curv_sub(ct_sub==0) = NaN;
                mean_sub_curv(sub_i) = nanmean(curv_sub,1);
                std_sub = nanstd(curv_sub,1);
                z_intra_curv{sub_i}(:,1) = ((curv_sub - mean_sub_curv(sub_i))./std_sub);
            end
            
            save(fullfile(resultsDir, ['CT_zscoreWB_' num2str(smoothing(smooth_i)) 'FWHM.mat']), 'z_intra_ct', '-v7.3');
            save(fullfile(resultsDir, ['Curv_zscoreWB_' num2str(smoothing(smooth_i)) 'FWHM.mat']), 'z_intra_curv', '-v7.3');

        end
    else
        
        file_ct = fullfile(resultsDir, ['CT_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']);
        file_cv = fullfile(resultsDir, ['Curv_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']);
       
        if exist(file_ct, 'file') ~= 2 && exist(file_cv, 'file') ~= 2
            
            load(fullfile(resultsDir, ['CorticalThickness_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
            load(fullfile(resultsDir, ['Curvature_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
            
            for sub_i = 1:length(subjects)
            
                for hemi_i = 1:length(hemi)
                    ct_sub = ct{sub_i}(:,hemi_i);
                    ct_sub(ct_sub==0) = NaN;
                    mean_sub_ct(sub_i) = nanmean(ct_sub,1);
                    std_sub_ct(sub_i) = nanstd(ct_sub,1);
                    z_intra(:,hemi_i) = ((ct_sub - mean_sub_ct(sub_i))./std_sub_ct(sub_i));
                end
                        
                z_intra_ct{sub_i}(:,1) = [z_intra(:,1); z_intra(:,2)];
                        
                clear z_intra;
            
                for hemi_i = 1:length(hemi)
                    curv_sub = abs(curv{sub_i}(:,hemi_i));
                    curv_sub(curv_sub==0) = NaN;
                    mean_sub_curv(sub_i) = nanmean(curv_sub,1);
                    std_sub_curv(sub_i) = nanstd(curv_sub,1);
                    z_intra(:,hemi_i) = ((curv_sub - mean_sub_curv(sub_i))./std_sub_curv(sub_i));
                end
                
                z_intra_curv{sub_i}(:,1) = [z_intra(:,1); z_intra(:,2)];

            end
            
            save(fullfile(resultsDir, ['CT_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']), 'z_intra_ct', '-v7.3');
            save(fullfile(resultsDir, ['Curv_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']), 'z_intra_curv', '-v7.3');
           
            clear z_intra;
        end
    end

end


%% Distribution of controls and patients' standardized data

for smooth_i = 1:length(smoothing)

    disp(['Violin plot - Smoothing ' num2str(smoothing(smooth_i)) ' FWHM']);
    
    % load the data
    load(fullfile(resultsDir, ['T1_intensities_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']));    
    load(fullfile(resultsDir, ['CT_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']));
    load(fullfile(resultsDir, ['Curv_zscore_' num2str(smoothing(smooth_i)) 'FWHM.mat']));

    % create the histogram of z-scores
    for sub_i = 1:length(subjects)
        t1_all(:,sub_i) = z_intra_t1{:,sub_i}(:,1); 
        ct_all(:,sub_i) = z_intra_ct{:,sub_i}(:,1); 
        cv_all(:,sub_i) = z_intra_curv{:,sub_i}(:,1); 
    end

    group_data_ct = [mean(ct_all(:, 1:172), 2), mean(ct_all(:, 173:end), 2)];
    group_data_cv = [mean(cv_all(:, 1:172), 2), mean(cv_all(:, 173:end), 2)];
    group_data_t1 = [mean(t1_all(:, 1:172), 2), mean(t1_all(:, 173:end), 2)];
    group_labels = [repmat({'Sighted'}, length(ct_all(:, 1:172)), 1); repmat({'Sight-Recovery'}, length(ct_all(:, 173:end)), 1)];

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

    exportgraphics(gcf, fullfile(resultsDir, 'figures', ['violinplot_zscores_' num2str(smoothing(smooth_i)) 'FWHM.pdf']), 'ContentType', 'vector')

    clear t1_all ct_all cv_all
    
end

close all;
