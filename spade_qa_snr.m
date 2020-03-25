% This script reads in SNR values and plots them according to session (pre-
% vs post-training) and group (expert=3, beginner=2, control=1).

clear all; close all; clc
format shortG

remove_outliers = 'no';
include = 'all';

blprojectid = 'proj-5e61139282b37f2cfe8fdb28';

% Set working directories.
rootDir = '/Volumes/240/spade/';

% Read in behavioral data.
beh_data_in_tbl = readtable([rootDir 'supportFiles/SPADE_demographics.csv'], 'TreatAsEmpty', {'.', 'na'});

% Identify outliers to be removed.
% outlier = [108 315 318];% 128 315 318];

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
        
        % Get contents of the directory where the SNR values for this subject are stored.
        sub_contents_snr = dir(fullfile(grp_contents(t).folder, grp_contents(t).name, '/dt-raw.tag-snr*/*product.json'));
        % Remove the '.' and '..' files.
        sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
        
        % Get SNR for this subject.
        data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
        
        for g = 1:size(data_snr_temp.SNRInB0_X_Y_Z, 1)
            
            get_temp(g) = str2num(data_snr_temp.SNRInB0_X_Y_Z{g});
            
        end
        
        % Get subID.
        subID(sub_count) = str2num(grp_contents(t).name(5:7));
        
        % Get session.
        ses(sub_count) = str2num(grp_contents(t).name(end));
        
        % Get SNR.
        snr(sub_count) = min(get_temp);
        
        % Get age group.
        group(sub_count) = beh_data_in_tbl.DanceLevelCode(find((beh_data_in_tbl.No == str2num(grp_contents(t).name(5:7)))));
        
        % Get age in months.
        age(sub_count) = beh_data_in_tbl.Age(find((beh_data_in_tbl.No == str2num(grp_contents(t).name(5:7)))));  
        
        clear data_snr_temp get_temp
        
    end % end if exist
    
end % end t

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(subID, outlier);
    
    % Remove outliers.
    subID = subID(~idx_outlier);
    snr = snr(~idx_outlier);
    group = group(~idx_outlier);
    age = age(~idx_outlier);
    ses = ses(~idx_outlier);
    
    % Set figure note.
    ttlstr = 'SNR outlier removed.';
    
else
    
    % Set figure note.
    ttlstr = 'SNR outlier retained.';
    
end

% Write out table for anova.
t_out = array2table(cat(2, subID', ses', group', age', snr'), 'VariableNames', {'subID', 'ses', 'group', 'cov_age', 'snr'});
writetable(t_out, fullfile(rootDir, 'supportFiles', ['spade_data_snr_' include '.csv']));

% Group differences test

disp('Are there SNR differences between groups at session 1?')
snrtotest = snr(ses == 1)';
grouptotest = group(ses == 1)';
[p, tableout, stats] = anova1(snrtotest, grouptotest, 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
clear grouptotest snrtotest

disp('Are there SNR differences between groups at session 2?')
snrtotest = snr(ses == 2)';
grouptotest = group(ses == 2)';
[p, tableout, stats] = anova1(snrtotest, grouptotest, 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
clear grouptotest snrtotest

disp('Are there SNR differences between sessions in Experts?')
snrtotest = snr(group == 3)';
sestotest = ses(group == 3)';
[p, tableout, stats] = anova1(snrtotest, sestotest, 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
clear sesttest snrtotest

disp('Are there SNR differences between sessions in Beginners?')
snrtotest = snr(group == 2)';
sestotest = ses(group == 2)';
[p, tableout, stats] = anova1(snrtotest, sestotest, 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
clear sesttest snrtotest

disp('Are there SNR differences between sessions in Controls?')
snrtotest = snr(group == 1)';
sestotest = ses(group == 1)';
[p, tableout, stats] = anova1(snrtotest, sestotest, 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
clear sesttest snrtotest

% Group differences plot
figure(1)
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
ylim_lo = 5;
ylim_hi = 16;

c_color = [0 0.4470 0.7410];
b_color = [0.4660 0.6740 0.1880];
e_color = [0.6350 0.0780 0.1840];

% Controls
b1 = bar(1, nanmean(snr(group == 1)), 'FaceColor', c_color, 'EdgeColor', c_color, 'FaceAlpha', alphablend);
plot([1 1], [nanmean(snr(group == 1)) - nanstd(snr(group == 1)) nanmean(snr(group == 1)) + nanstd(snr(group == 1))], 'Color', c_color)
% Beginners
b2 = bar(2, nanmean(snr(group == 2)), 'FaceColor', b_color, 'EdgeColor', b_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(snr(group == 2)) - nanstd(snr(group == 2)) nanmean(snr(group == 2)) + nanstd(snr(group == 2))], 'Color', b_color)
% Experts
b3 = bar(3, nanmean(snr(group == 3)), 'FaceColor', e_color, 'EdgeColor', e_color, 'FaceAlpha', alphablend);
plot([3 3], [nanmean(snr(group == 3)) - nanstd(snr(group == 3)) nanmean(snr(group == 3)) + nanstd(snr(group == 3))], 'Color', e_color)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 3.5];
xax.TickValues = [1 2 3];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'Controls', 'Beginners', 'Experts'};
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

a.YLabel.String = 'SNR';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', ['plot_barplot_snr_bygroup' include]), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', ['plot_barplot_snr_bygroup' include]), '-depsc')

hold off;

% Session differences plot
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
ylim_lo = 5;
ylim_hi = 16;

s1_color = [0.41176 0.41176 0.41176];
s2_color = [0.82745 0.82745 0.82745];

% Session 1
b1 = bar(1, nanmean(snr(ses == 1)), 'FaceColor', s1_color, 'EdgeColor', s1_color, 'FaceAlpha', alphablend);
plot([1 1], [nanmean(snr(ses == 1)) - nanstd(snr(ses == 1)) nanmean(snr(ses == 1)) + nanstd(snr(ses == 1))], 'Color', s1_color)
% Session 2
b2 = bar(2, nanmean(snr(ses == 2)), 'FaceColor', s2_color, 'EdgeColor', s2_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(snr(ses == 2)) - nanstd(snr(ses == 2)) nanmean(snr(ses == 2)) + nanstd(snr(ses == 2))], 'Color', s2_color)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 2.5];
xax.TickValues = [1 2];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'Pretraining', 'Posttraining'};
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

a.YLabel.String = 'SNR';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', ['plot_barplot_snr_bysession_' include]), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', ['plot_barplot_snr_bysession_' include]), '-depsc')

hold off;

% Session x Group differences plot
figure(3)
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
alphablend = .9; alphablend2 = .4;
ylim_lo = 5;
ylim_hi = 16;

% Session 1
b1 = bar(1, nanmean(snr(ses == 1 & group == 1)), 'FaceColor', c_color, 'EdgeColor', c_color, 'FaceAlpha', alphablend, 'EdgeAlpha', alphablend2);
plot([1 1], [nanmean(snr(ses == 1 & group == 1)) - nanstd(snr(ses == 1 & group == 1)) nanmean(snr(ses == 1 & group == 1)) + nanstd(snr(ses == 1 & group == 1))], 'Color', c_color)
b2 = bar(2, nanmean(snr(ses == 1 & group == 2)), 'FaceColor', b_color, 'EdgeColor', b_color, 'FaceAlpha', alphablend, 'EdgeAlpha', alphablend2);
plot([2 2], [nanmean(snr(ses == 1 & group == 2)) - nanstd(snr(ses == 1 & group == 2)) nanmean(snr(ses == 1 & group == 2)) + nanstd(snr(ses == 1 & group == 2))], 'Color', b_color)
b3 = bar(3, nanmean(snr(ses == 1 & group == 3)), 'FaceColor', e_color, 'EdgeColor', e_color, 'FaceAlpha', alphablend, 'EdgeAlpha', alphablend2);
plot([3 3], [nanmean(snr(ses == 1 & group == 3)) - nanstd(snr(ses == 1 & group == 3)) nanmean(snr(ses == 1 & group == 3)) + nanstd(snr(ses == 1 & group == 3))], 'Color', e_color)
% Session 2
b4 = bar(4, nanmean(snr(ses == 2 & group == 1)), 'FaceColor', c_color, 'EdgeColor', c_color, 'FaceAlpha', alphablend2, 'EdgeAlpha', alphablend2);
plot([4 4], [nanmean(snr(ses == 2 & group == 1)) - nanstd(snr(ses == 2 & group == 1)) nanmean(snr(ses == 2 & group == 1)) + nanstd(snr(ses == 2 & group == 1))], 'Color', c_color)
b5 = bar(5, nanmean(snr(ses == 2 & group == 2)), 'FaceColor', b_color, 'EdgeColor', b_color, 'FaceAlpha', alphablend2, 'EdgeAlpha', alphablend2);
plot([5 5], [nanmean(snr(ses == 2 & group == 2)) - nanstd(snr(ses == 2 & group == 2)) nanmean(snr(ses == 2 & group == 2)) + nanstd(snr(ses == 2 & group == 2))], 'Color', b_color)
b6 = bar(6, nanmean(snr(ses == 2 & group == 3)), 'FaceColor', e_color, 'EdgeColor', e_color, 'FaceAlpha', alphablend2, 'EdgeAlpha', alphablend2);
plot([6 6], [nanmean(snr(ses == 2 & group == 3)) - nanstd(snr(ses == 2 & group == 3)) nanmean(snr(ses == 2 & group == 3)) + nanstd(snr(ses == 2 & group == 3))], 'Color', e_color)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 6.5];
xax.TickValues = [1 2 3 4 5 6];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'Pre, Cont.', 'Pre, Beg.', 'Pre, Exp.', ...
    'Post, Cont.', 'Post, Beg.', 'Post, Exp.',};
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

a.YLabel.String = 'SNR';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', ['plot_barplot_snr_bysessionxgroup_' include]), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', ['plot_barplot_snr_bysessionxgroup_' include]), '-depsc')

hold off;