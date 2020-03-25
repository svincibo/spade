% This script reads in streamline count values (i.e., nfibers) and 
% checks that the number of streamlines in each tract is correlated across
% subjects between sessions and checks that there are no significant
% differences in streamline count within tracts between groups at either session. 
% A box plot is provided to help identify unusually low streamline counts within 
% a particular tract and to compare sessions. 

clear all; close all; clc
format shortG

remove_outliers = 'no';
% Identify outliers to be removed.
% outlier = [107];% subID of outliers to be removed

stat = 'nfibers';

blprojectid = 'proj-5e61139282b37f2cfe8fdb28';

% Set working directories.
rootDir = '/Volumes/240/spade/';

% Read in behavioral data.
beh_data_in_tbl = readtable([rootDir 'supportFiles/SPADE_demographics.csv'], 'TreatAsEmpty', {'.', 'na'});

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
sub_count = 0;
for t = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behaviora/demographic data.
    if ~isempty(find((beh_data_in_tbl.No == str2num(grp_contents(t).name(5:7)))))
        
        % Display current sub ID.
        disp(grp_contents(t).name)
        
        % Update subject counter for when not all subjects are used/needed.
        sub_count = sub_count + 1;
            
            % Get contents of the directory where the tract measures for this subject are stored.
            sub_contents_tractstats = dir(fullfile(grp_contents(t).folder, grp_contents(t).name, '/dt-neuro-tractmeasures*/*.csv'));
            
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
            
            for t2 = 1:length(toi)
                
                % Read in mean stat.
                m(t, t2) = data_tbl_in.StreamlineCount(toi(t2));
                
                % Grab tract name for grouping variable.
                tract{t, t2} = char(data_tbl_in.TractName(toi(t2)));
                
            end % end t
            
            % Grab subID.
            sub(sub_count) = str2num(grp_contents(t).name(5:7));
            
            % Gather session, for ease.
            ses(sub_count) = str2num(grp_contents(t).name(end));
            
            % Get age group.
            group(sub_count) = beh_data_in_tbl.DanceLevelCode(find((beh_data_in_tbl.No == str2num(grp_contents(t).name(5:7)))));
            
            % Get age in months.
            age(sub_count) = beh_data_in_tbl.Age(find((beh_data_in_tbl.No == str2num(grp_contents(t).name(5:7)))));
            
            clear data_temp toi toi_idx
        
    end % end if exist
    
end % end t

% Find empty cells.
idx = find(cellfun(@isempty,tract));

% Enter 'empty' in empty cells.
tract(idx) = {'empty'};

% Enter NaN for m in empty cells.
m(idx) = NaN;

% Get a list of unique tract names.
list_tract = tract(1, :);
list_tract = list_tract(~strcmp(list_tract, 'empty'));

% % Determine which subIDs appear in both WM and BEH.
sub_tract_beh = intersect(sub, beh_data_in_tbl.No);
% 
% % Get indices of subjects who appear in both WM and BEH.
sub_idx_wm = ismember(sub, sub_tract_beh);
% sub_idx_beh = ismember(beh_data_in_tbl.No, sub_tract_beh);

% Select only subjects who appear in both WM and BEH.
% Concatenate into one data array and one header array.
% Remove redundant subID columns.
data_out = cat(2, sub', ses', group', age', m(sub_idx_wm, :));
data_out_header = [{'subID', 'ses', 'group', 'cov_age'}, list_tract];

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

% Write out table.
writetable(data_tbl, fullfile(rootDir, 'supportFiles', ['spade_data_' stat '.csv']));

% Group differences test
d = m(sub_idx_wm, :);
for tn = 1:size(d, 2)

    disp(list_tract{tn});
    
    disp('Are there streamline count differences between groups at session 1?')
    tracttotest = d(ses == 1, tn);
    grouptotest = group(ses == 1)';
    [p, tableout, stats] = anova1(tracttotest, grouptotest, 'off');
    disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
    
    disp('Are there streamline count differences between groups at session 2?')
    tracttotest = d(ses == 2, tn);
    grouptotest = group(ses == 2)';
    [p, tableout, stats] = anova1(tracttotest, grouptotest, 'off');
    disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
    
    disp('Are there streamline count differences between sessions in Experts?')
    tracttotest = d(group == 1, tn);
    sestotest = ses(group == 1)';
    [p, tableout, stats] = anova1(tracttotest, sestotest, 'off');
    disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
    
    disp('Are there streamline count differences between sessions in Beginners?')
    tracttotest = d(group == 1, tn);
    sestotest = ses(group == 1)';
    [p, tableout, stats] = anova1(tracttotest, sestotest, 'off');
    disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
    
    disp('Are there streamline count differences between sessions in Controls?')
    tracttotest = d(group == 1, tn);
    sestotest = ses(group == 1)';
    [p, tableout, stats] = anova1(tracttotest, sestotest, 'off');
    disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
    
end % end tn

% Check correlation between session 1 and session 2.
d = m(sub_idx_wm, :);
for tn = 1:size(d, 2)
    
    [rho(tn), p(tn)] = corr(d(ses==1, tn), d(ses == 2, tn), 'type', 'Pearson', 'rows', 'complete', 'tail', 'both');
    
    disp(['Between-session correlation for ' list_tract{tn} ': rho = ' num2str(rho(tn)) ', p = ' num2str(p(tn)) '.']);
    
end

% Look at boxplot of streamline count.
prepbox = NaN(size(d, 1)/2, size(d, 2));
prepbox(:, 1:2:size(d, 2)*2) = d(ses==1, :);
prepbox(:, 2:2:size(d, 2)*2) = d(ses==2, :);
labelsforx = cell(size(d, 2)*2, 1);
labelsforx(1:2:size(d, 2)*2) = strcat('s1-', list_tract);
labelsforx(2:2:size(d, 2)*2) = strcat('s2-', list_tract);
figure(2)
boxplot(prepbox, 'colors', 'rb', 'Labels', labelsforx)
xtickangle(90)
ylabel('Streamline Count')
xlabel('Each Tract at Each Session')




