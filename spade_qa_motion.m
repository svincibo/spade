clear all; close all; clc
format shortG

expert_color = [0.6350 0.0780 0.1840]; %red
beginner_color = [0 0.4470 0.7410]; %blue
control_color = [0.41176 0.41176 0.41176]; %gray

blprojectid = 'proj-5e61139282b37f2cfe8fdb28';

% Set working directories.
rootDir = '/Volumes/240/spade/';
% addpath(genpath(fullfile(rootDir, 'proj-5e5672430f7fa65e1d3c9621')));

remove_outliers = 'no';
outlier = [];

%%%%%%%%%%%%%%% TESTING AREA %%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% For each subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behaviora/demographic data.
    if ~isempty(find((beh_data_in_tbl.No == str2num(grp_contents(s).name(5:7)))))
        
        % Display current sub ID.
        disp(grp_contents(s).name)
        
        % Get contents of the directory where the SNR values for this subject are stored.
        sub_contents_snr = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/*product.json'));
        % Remove the '.' and '..' files.
        sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
        
        if ~isempty(sub_contents_snr)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get contents of the directory where the SNR values for this subject are stored.
            sub_contents_motion = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-neuro-dtiinit*/*ecXform.mat'));
            % Remove the '.' and '..' files.
            sub_contents_motion = sub_contents_motion(arrayfun(@(x) x.name(1), sub_contents_motion) ~= '.');
            
            % Get motion parameters for this subject.
            load(fullfile(sub_contents_motion.folder,sub_contents_motion.name));
            
            % Select only the translation/rotation parameters.
            mot = vertcat(xform(:).ecParams);
            mot = mot(:, 1:6); % xyz translation xyz rotation
            
            % Motion parameters represent the translation/rotation that must occur
            % to return the image at timepoint tn to the place that it was at timepoint
            % t0. Thus, to calculate instantaneouse parameters, we need a moving
            % difference.
            for m = 1:size(mot, 2)
                
                % Get moving difference for each parameter. Append row of zeros for t1, per convention (Power et al., 2014).
                % This step creates the fd timecourses for each motion parameter for each run in the motion struct and represents instantaneous movement.
                movingdifference(:, m) = [0 ; diff(mot(:, m), 1, 1)]';
                
            end
            
            % Get an overall fd for all 6 parameters for each run.
            % This step creates the fd summary statistic for all 6 parameters for each timepoint in a run for each run (e.g., scalar FD timecourse).
            motion(sub_count, :) = sum(abs(movingdifference), 2)';
            
            % Get subID.
            subID(sub_count) = str2num(grp_contents(s).name(5:7));
            
            % Get session.
            ses(sub_count) = str2num(grp_contents(s).name(end));
            
            % Get ylabel.
            lab{sub_count} = grp_contents(s).name;
            
            % Get training group.
            group(sub_count) = beh_data_in_tbl.DanceLevelCode(find((beh_data_in_tbl.No == str2num(grp_contents(s).name(5:7)))));
            
        end % if ~isempty
        
    end % end if exist
    
end % end s

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(subID, outlier);
    
    % Remove outliers.
    subID = subID(~idx_outlier);
    motion = motion(~idx_outlier, :);
    group = group(~idx_outlier);
    ses = ses(~idx_outlier);
    
    % Set figure note.
    ttlstr = 'Motion outlier removed.';
    
else
    
    % Set figure note.
    ttlstr = 'Motion outlier retained.';
    
end

meanmotion = mean(motion, 2);

% Write out table for anova.
t_out = array2table(cat(2, subID', ses', group', meanmotion), 'VariableNames', {'subID', 'session', 'group', 'fd'});
writetable(t_out, fullfile(rootDir, 'supportFiles', 'spade_data_motion.csv'));

disp('Check for group differences in FD.')
[~, tableout, ~] = anova1(meanmotion, group', 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])

% Visualize: group differences
figure(2)
hold on;

capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;
xtickvalues = [1 2 3];
alphablend = .8;

b1 = bar(1, nanmean(meanmotion(group == 1)), 'FaceColor', control_color, 'EdgeColor', control_color, 'FaceAlpha', alphablend);
plot([1 1], [nanmean(meanmotion(group == 1)) - nanstd(meanmotion(group == 1)) nanmean(meanmotion(group == 1)) + nanstd(meanmotion(group == 1))], 'Color', control_color)
b2 = bar(2, nanmean(meanmotion(group == 2)), 'FaceColor', beginner_color, 'EdgeColor', beginner_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(meanmotion(group == 2)) - nanstd(meanmotion(group == 2)) nanmean(meanmotion(group == 2)) + nanstd(meanmotion(group == 2))], 'Color', beginner_color)
b3 = bar(3, nanmean(meanmotion(group == 3)), 'FaceColor', expert_color, 'EdgeColor', expert_color, 'FaceAlpha', alphablend);
plot([3 3], [nanmean(meanmotion(group == 3)) - nanstd(meanmotion(group == 3)) nanmean(meanmotion(group == 3)) + nanstd(meanmotion(group == 3))], 'Color', expert_color)

% xlim_lo = min(age)-1;
% xlim_hi = max(age)+1;
ylim_lo = 0;
ylim_hi = max(meanmotion)+0.5;
    
% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 3.5];
xax.TickValues = [1 2 3];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'Control', 'Beginner', 'Experienced'};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), '', num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'Mean Framewise Displacement (FD)';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_barplot_motion'), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', 'plot_barplot_motion'), '-depsc')

hold off;

% % Manually record outliers. Include observations with unusually high FD and any above 2mm. 
% % (0 indicates no outliers)
% outliers.motion = 90;
% % outliers.motion = subID(meanmotion>2);
% 
% save('devti_remove_motionoutliers.mat', 'outliers')
