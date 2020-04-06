% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

remove_outliers = 'yes';

prefix = ''; % 'streamlinecount25920-'
blprojectid = ['/' prefix 'proj-5a74d1e26ed91402ce400cca/'];

w_measures = {'fa', 'ad', 'md', 'rd'};

blprojectid = 'proj-5e61139282b37f2cfe8fdb28';

% Set working directories.
rootDir = '/Volumes/240/spade/';

% Read in behavioral data.
beh_data_in_tbl = readtable([rootDir 'supportFiles/SPADE_demographics.csv'], 'TreatAsEmpty', {'.', 'na'});

% Parse table into array and header to make things easier.
% beh_data_in = table2array(beh_data_in_tbl);
beh_data_in_header = beh_data_in_tbl.Properties.VariableNames;

% Identify outliers to be removed.
% outlier = [128 315 318];

for w = 1:length(w_measures)
    
    wm_measure = w_measures{w};
    
    observation_count = 0;
    
    % Each session one at a time.
    for s = 1:2
        
        %% TRACTOGRAPHY.
        
        % Get contents of the directory where the tract measures for this subject are stored.
        grp_contents = dir(fullfile(rootDir, blprojectid));
        
        % Remove the '.' and '..' files.
        grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
        
        % Keep only names that are subject folders.
        grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's' & arrayfun(@(x) x.name(end), grp_contents) == num2str(s));
        
        % Load in each tract's tractography measures for this subject.
        for i = 1:size(grp_contents, 1)
            
            % Make sure behavioral data exist.
            idx = find(beh_data_in_tbl.No == str2num(grp_contents(i).name(5:7)));
            
            if exist('idx', 'var') & ~isempty(idx)
                
                observation_count = observation_count + 1;
                
                % Grab subID.
                sub(observation_count) = str2num(grp_contents(i).name(5:7));
                
                % Grab session.
                ses(observation_count) = s;
                
                % Grab other behavioral variables from beh file.
                % Get group.
                dance_group(observation_count) = beh_data_in_tbl.DanceLevelCode(idx);
                
                % Get age.
                age(observation_count) = beh_data_in_tbl.Age(idx);
                
                % Get practice hours.
                hours(observation_count) = beh_data_in_tbl.TotalPracticeHours(idx);
                
                % Get sex.
                sex(observation_count) = beh_data_in_tbl.gender(idx);
                
                % Display current sub ID.
                disp(grp_contents(i).name)
                
                % Get contents of the directory where the tract measures for this subject are stored.
                sub_contents_tractprofiles = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '/dt-neuro-tractprofile*/profiles/*.csv'));
                
                % Remove the '.' and '..' files.
                sub_contents_tractprofiles = sub_contents_tractprofiles(arrayfun(@(x) x.name(1), sub_contents_tractprofiles) ~= '.');
                
                for j = 1:size(sub_contents_tractprofiles)
                    
                    % Preallocate based on number of subjects(size(grp_contents)) and number of tracts (size(sub_contents...)).
                    if i == 1 && j == 1
                        
                        tract = {}; m = NaN(size(grp_contents, 1), size(sub_contents_tractprofiles, 1));
                        
                    end
                    
                    % Read in data for this subject and this tract.
                    data_temp = readtable([sub_contents_tractprofiles(j).folder filesep sub_contents_tractprofiles(j).name]);
                    
                    % Get middle 80%.
                    start = size(data_temp, 1)*.1;
                    stop = size(data_temp, 1)*.9;
                    
                    % Read in mean WM measure.
                    if strcmp(wm_measure, 'ad')
                        
                        m(observation_count, j) = nanmean(data_temp.ad_1(start:stop));
                        
                    elseif strcmp(wm_measure, 'fa')
                        
                        m(observation_count, j) = nanmean(data_temp.fa_1(start:stop));
                        
                    elseif strcmp(wm_measure, 'md')
                        
                        m(observation_count, j) = nanmean(data_temp.md_1(start:stop));
                        
                    elseif strcmp(wm_measure, 'rd')
                        
                        m(observation_count, j) = nanmean(data_temp.rd_1(start:stop));
                        
                    end
                    
                    % Grab tract name for grouping variable.
                    tract{observation_count, j} = sub_contents_tractprofiles(j).name(1:end-13);
                    
                    clear data_temp
                    
                end % end j
                
            end % end exist
            
        end % end i
        
    end % end s
    
    % Find empty cells.
    t = find(cellfun(@isempty,tract));
    
    % Enter 'empty' in empty cells.
    tract(t) = {'empty'};
    
    % Get a list of unique tract names.
    list_tract = unique(tract);
    
    % Get WM measurements for each tract (reorganizing so that each column is a tract).
    for k = 1:size(list_tract, 1)
        
        % Select the wm_measurements for this tract from each subject.
        temp = m.*strcmp(tract, list_tract{k});
        
        % Convert all zeros to NaN;
        temp(temp == 0) = NaN;
        
        % Get the mean of the wm_measure for this tract (take sum because don't want to include zero columns; only one value of interest per row).
        wm(:, k) = nansum(temp, 2);
        
        clear temp
        
    end % end k
    
    % Convert all zeros to NaN;
    wm(wm == 0) = NaN;
    
    % Remove 'empty' column from data and header and append subID.
    wm = cat(2, transpose(sub), wm(:, find(~all(isnan(wm), 1))));
    wm_header = [{'subID'}, transpose(list_tract(~strcmp(list_tract, 'empty')))];
    
    % Create grouping and behavioral vectors.
    %     beh = cat(2, beh_data_in_tbl.No, beh_data_in_tbl.DanceLevelCode, beh_data_in_tbl.Age, beh_data_in_tbl.TotalPracticeHours, beh_data_in_tbl.gender);
    
    % Select only subjects who appear in both WM and BEH.
    % Concatenate into one data array and one header array.
    % Remove redundant subID columns.
    data_all = cat(2, sub', ses', dance_group', age', hours', sex', wm(:, find(strcmp(wm_header, 'subID'))+1:end));
    data_all_header = [{'subID',  'session', 'group', 'cov_age', 'cov_practicehours', 'cov_sex'}, wm_header{find(strcmp(wm_header, 'subID'))+1:end}];
    
    % Remove outliers.
    if strcmp(remove_outliers, 'yes') && exist('outlier')
        
        % Get index for outliers to be removed.
        idx_outlier = ismember(data_all(:, find(strcmp(data_all_header, {'subID'}))), outlier);
        
        % Remove outliers.
        data_all = data_all(~idx_outlier, :);
        
    end
    
    data_tbl = array2table(data_all, 'VariableNames', data_all_header);
    
    % Save all variables.
    save([rootDir 'supportFiles/spade_data_' wm_measure '.mat'], 'data_tbl')
    
    % Write out table.
    writetable(data_tbl, fullfile(rootDir, 'supportFiles', ['spade_data_' wm_measure '.csv']));
    
    % Reset for next loop.
    clearvars -except w rootDir beh_data_in_tbl beh_data_in_header beh_data_in blprojectid remove_outliers w_measures outlier
    
end