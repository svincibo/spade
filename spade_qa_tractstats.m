% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

remove_outliers = 'no';

stats = 'nfibers'; % aka number of streamlines

blprojectid = 'proj-5e61139282b37f2cfe8fdb28';

% Set working directories.
rootDir = '/Volumes/240/spade/';

% Read in behavioral data.
beh_data_in_tbl = readtable([rootDir 'supportFiles/SPADE_demographics.csv'], 'TreatAsEmpty', {'.', 'na'});

% Parse table into array and header to make things easier.
beh_data_in = table2array(beh_data_in_tbl);
beh_data_in_header = beh_data_in_tbl.Properties.VariableNames;

% Identify outliers to be removed.
% outlier = [108 128 203 315 316 318];

%% TRACTOGRAPHY.

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
for i = 1:size(grp_contents, 1)
    
    % Grab subID.
    sub(i) = str2num(grp_contents(i).name(5:7));
    
    % Gather session, for ease.
    ses(i) = str2num(grp_contents(i).name(end));
    
    % Display current session and sub ID.
    disp(['Session ' s' ', ' grp_contents(i).name])
    
    % Get contents of the directory where the tract measures for this subject are stored.
    sub_contents_tractstats = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, '/dt-neuro-tractmeasures*/*.csv'));
    
    % Remove the '.' and '..' files.
    sub_contents_tractstats = sub_contents_tractstats(arrayfun(@(x) x.name(1), sub_contents_tractstats) ~= '.');
    
    % Read in data for this subject and this tract.
    data_tbl_in = readtable(fullfile(sub_contents_tractstats.folder, sub_contents_tractstats.name));
    
    % Convert into header for ease.
    data_all_in_header = data_tbl_in.TractName;
    
    % Get index matrices for hypothesis-driven grouping of WM tracts.
    for k = 1:length(data_all_in_header)
        
        % Indices of horizontal tracts.
        toi_idx(k) = strcmp(data_all_in_header{k}, 'leftSLF1And2') || strcmp(data_all_in_header{k}, 'rightSLF1And2') ...
            || strcmp(data_all_in_header{k}, 'leftIFOF') || strcmp(data_all_in_header{k}, 'rightIFOF') ...
            || strcmp(data_all_in_header{k}, 'leftILF') || strcmp(data_all_in_header{k}, 'rightILF') ...
            || strcmp(data_all_in_header{k}, 'leftArc') || strcmp(data_all_in_header{k}, 'rightArc') ...
            || strcmp(data_all_in_header{k}, 'leftSLF3') || strcmp(data_all_in_header{k}, 'rightSLF3') ...
            || strcmp(data_all_in_header{k}, 'leftAslant') || strcmp(data_all_in_header{k}, 'rightAslant') ...
            || strcmp(data_all_in_header{k}, 'leftTPC') || strcmp(data_all_in_header{k}, 'rightTPC') ...
            || strcmp(data_all_in_header{k}, 'leftpArc') || strcmp(data_all_in_header{k}, 'rightpArc') ...
            || strcmp(data_all_in_header{k}, 'leftMDLFspl') || strcmp(data_all_in_header{k}, 'rightMDLFspl') ...
            || strcmp(data_all_in_header{k}, 'leftVOF') || strcmp(data_all_in_header{k}, 'rightVOF') ...
            || strcmp(data_all_in_header{k}, 'leftMDLFang') || strcmp(data_all_in_header{k}, 'rightMDLFang');
        
    end % end k
    
    % Get positions of tracts of interest in data_tbl.
    toi = find(toi_idx == 1);
    
    for t = 1:length(toi)
        
        % Read in mean stat.
        m(i, t) = data_tbl_in.StreamlineCount(toi(t));
        
        % Grab tract name for grouping variable.
        tract{i, t} = char(data_tbl_in.TractName(toi(t)));
        
    end % end t
    
    clear data_temp toi toi_idx
    
end % end i

% Find empty cells.
idx = find(cellfun(@isempty,tract));

% Enter 'empty' in empty cells.
tract(idx) = {'empty'};

% Enter NaN for m in empty cells.
m(idx) = NaN;

% Get a list of unique tract names.
list_tract = tract(1, :);
list_tract = list_tract(~strcmp(list_tract, 'empty'));

% Determine which subIDs appear in both WM and BEH.
sub_tract_beh = intersect(sub, beh_data_in_tbl.SubjectID);

% Get indices of subjects who appear in both WM and BEH.
sub_idx_wm = ismember(sub, sub_tract_beh);
sub_idx_beh = ismember(beh_data_in_tbl.SubjectID, sub_tract_beh);

% Select only subjects who appear in both WM and BEH.
% Concatenate into one data array and one header array.
% Remove redundant subID columns.
data_out = cat(2, sub', ses', beh_data_in_tbl.group_age(sub_idx_beh, :), beh_data_in_tbl.Age_months(sub_idx_beh, :), m(sub_idx_wm, :));
data_out_header = [{'subID', 'session', 'gp_age', 'cov_age'}, list_tract];

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(data_out(:, find(strcmp(data_out_header, {'subID'}))), outlier);
    
    % Remove outliers.
    data_out = data_out(~idx_outlier, :);
    
end

data_tbl = array2table(data_out, 'VariableNames', data_out_header);

% Save all variables.
save([rootDir 'supportFiles/spade_data_' stat '.mat'], 'data_tbl')









