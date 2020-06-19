% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/spade/';

% Get bl project foldername.
blprojectid = 'proj-5e61139282b37f2cfe8fdb28';

w_measures = {'fa', 'md'};

exp_color = [0.6350 0.0780 0.1840]; %red
beg_color = [0 0.4470 0.7410]; %blue
con_color = [0.41176 0.41176 0.41176]; %gray

hold on;
linewidth = 1.5;
linestyle = 'none';
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
xticklength = 0;

save_figures = 'yes';

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'yes';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed - e.g., outlier = [108 126 212 214 318];  
    outlier = [];
    
else
    
    outlier = [];
    
end

% Read in behavioral data.
beh_data_in_tbl = readtable([rootDir 'supportFiles/SPADE_demographics.csv'], 'TreatAsEmpty', {'.', 'na'});

figcount = 0;
for w = 2:length(w_measures)
    
    wm_measure = w_measures{w};
    
    if strcmp(wm_measure, 'fa')
        ylim_lo = 0.20; ylim_hi = 0.70;
        ylim_diff_lo = -0.25; ylim_diff_hi = 0.25;
    elseif strcmp(wm_measure, 'md')
        ylim_lo = 0.60; ylim_hi = 1.20;
        ylim_diff_lo = -0.25; ylim_diff_hi = 0.25;
    end
    
    %% TRACTOGRAPHY.
    
    % Get contents of the directory where the tract measures for this subject are stored.
    grp_contents = dir(fullfile(rootDir, blprojectid));
    
    % Remove the '.' and '..' files.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
    
    % Keep only names that are subject folders.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');
    
    % Load in each tract's tractography measures for this subject.
    sub_count = 0;
    for i = 1:size(grp_contents, 1)
        
        % Only collect values for subjects that have both MRI and behaviora/demographic data.
        if ~isempty(find((beh_data_in_tbl.No == str2num(grp_contents(i).name(5:7)))))
            
            % Display current sub ID.
            disp(grp_contents(i).name)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get contents of the directory where the tract measures for this subject are stored.
            sub_contents_tractprofiles = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, 'dt-neuro-tractprofile*', 'profiles', '*.csv'));
            
            % Remove the '.' and '..' files.
            sub_contents_tractprofiles = sub_contents_tractprofiles(arrayfun(@(x) x.name(1), sub_contents_tractprofiles) ~= '.');
            
            for j = 1:size(sub_contents_tractprofiles)
                
                % Preallocate based on number of subjects(size(grp_contents)) and number of tracts (size(sub_contents...)).
                if i == 1 && j == 1
                    
                    tract = {}; m = NaN(size(grp_contents, 1), size(sub_contents_tractprofiles, 1));
                    
                end
                
                % Read in data for this subject and this tract.
                data_temp = readtable(fullfile(sub_contents_tractprofiles(j).folder, sub_contents_tractprofiles(j).name));
                
                % Get middle 80%.
                start = size(data_temp, 1)*.1;
                stop = size(data_temp, 1)*.9;
                
                % Read in mean WM measure.
                if strcmp(wm_measure, 'ad')
                    
                    m_wm(:, j, sub_count) = data_temp.ad_1(start:stop);
                    sd_wm(:, j, sub_count) = data_temp.ad_2(start:stop);
                    
                elseif strcmp(wm_measure, 'fa')
                    
                    m_wm(:, j, sub_count) = data_temp.fa_1(start:stop);
                    sd_wm(:, j, sub_count) = data_temp.fa_2(start:stop);
                    
                elseif strcmp(wm_measure, 'md')
                    
                    m_wm(:, j, sub_count) = data_temp.md_1(start:stop);
                    sd_wm(:, j, sub_count) = data_temp.md_2(start:stop);
                    
                    
                end
                
                % Grab tract name for grouping variable.
                tract(:, j, sub_count) = repmat({sub_contents_tractprofiles(j).name(1:end-13)}, 161, 1);
                
                % Grab subID.
                sub(:, j, sub_count) = repmat(str2num(grp_contents(i).name(5:7)), 161, 1);
                
                % Gather session, for ease.
                ses(:, j, sub_count) = repmat(str2num(grp_contents(i).name(end)), 161, 1);
                
                % Get exp group.
                group(sub_count) = beh_data_in_tbl.DanceLevelCode(find((beh_data_in_tbl.No == str2num(grp_contents(i).name(5:7)))));
                
                % Get age in months.
                age(sub_count) = beh_data_in_tbl.Age(find((beh_data_in_tbl.No == str2num(grp_contents(i).name(5:7)))));
                
                clear data_temp
                
            end % end if exist
            
        end % sub_contents
        
    end % group_contents
    
    % Find empty cells and fill with 'empty'.
    t = find(cellfun(@isempty,tract));
    tract(t) = {'empty'};
    
    % Get a list of unique tract names.
    list_tract = unique(tract);
    
    % Get a list of unique sub IDs.
    subID = unique(sub);
    
    % Plot tract profiles for each tract.
    for k = 1:size(list_tract, 1)
        
        % Only plot for tracts of interest and non-empty tracts.
        if strcmp(list_tract{k}, 'leftSLF1And2') || strcmp(list_tract{k}, 'rightSLF1And2') ...
                || strcmp(list_tract{k}, 'leftIFOF') || strcmp(list_tract{k}, 'rightIFOF') ...
                || strcmp(list_tract{k}, 'leftILF') || strcmp(list_tract{k}, 'rightILF') ...
                || strcmp(list_tract{k}, 'leftArc') || strcmp(list_tract{k}, 'rightArc') ...
                || strcmp(list_tract{k}, 'leftSLF3') || strcmp(list_tract{k}, 'rightSLF3') ...
                || strcmp(list_tract{k}, 'leftAslant') || strcmp(list_tract{k}, 'rightAslant') ...
                || strcmp(list_tract{k}, 'leftTPC') || strcmp(list_tract{k}, 'rightTPC') ...
                || strcmp(list_tract{k}, 'leftpArc') || strcmp(list_tract{k}, 'rightpArc') ...
                || strcmp(list_tract{k}, 'leftMDLFspl') || strcmp(list_tract{k}, 'rightMDLFspl') ...
                || strcmp(list_tract{k}, 'leftVOF') || strcmp(list_tract{k}, 'rightVOF') ...
                || strcmp(list_tract{k}, 'leftMDLFang') || strcmp(list_tract{k}, 'rightMDLFang') ...
                && ~strcmp(list_tract{k}, 'empty')
            
            % Find entries that are for this tract.
            t_idx = strcmp(tract, list_tract{k});
            
            % Open a new figure for this tract.
            figcount = figcount + 1;
            figure(figcount)
            
            exp_count = 0; beg_count = 0; con_count = 0;
            for s = 1:size(m_wm, 3)
                
                % Only include subjects who are not outliers.
                if sub(1, 1, s)~=0 && ~ismember(sub(1, 1, s), outlier)
                    
                    % Find entries that are for this subject.
                    s_idx = sub == sub(1, 1, s);
                    
                    for session = 1:2
                        
                        % Find entries that are for this session.
                        ses_idx = ses == session;
                        
                        % Subset the thing so that we only plot for this tract, subject, and session.
                        t_temp = m_wm(find(t_idx == 1 & s_idx == 1 & ses_idx == 1));
                        
                        if ~isempty(t_temp)
                            
                            % Different line styles for session 1 and session 2.
                            if session == 1
                                linestyle = '-';
                            elseif session == 2
                                linestyle = ':';
                            end
                            
                            % Code the plot for subject and keep data for inspection (yc, oc, a).
                            if group(s) == 3 % expert
                                
                                exp_count = exp_count + 1;

                                % Collect.
                                expert(:, exp_count) = t_temp;
                                expert_ses(:, exp_count) = session;
                                
                                
                            elseif group(s) == 2 % beginner
                                
                                beg_count = beg_count + 1;
                                
                                % Collect.
                                beg(:, beg_count) = t_temp;
                                beg_ses(:, beg_count) = session;
                                
                            else % control
                                
                                con_count = con_count + 1;
                                
                                % Collect.
                                con(:, con_count) = t_temp;
                                con_ses(:, con_count) = session;
                                
                            end
                            hold on;
                            
                            clear t_temp;
                            
                        end
                        
                    end
                    
                end % if subID{s} ~= 0
                
            end %if exist
            
            % Plot means and 95% confidence intervals (calculated from standard error: 1.96*SE). 
            subplot(1, 3, 1)
            xnew = expert(:, expert_ses == 1); c = exp_color;
            plot(nanmean(xnew, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', c(1:3))
            hi = nanmean(xnew, 2) + 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); lo = nanmean(xnew, 2) - 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); x = (1:size(nanmean(xnew, 2),1))';
            hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], c(1:3));
            set(hp1, 'facecolor', c(1:3), 'edgecolor', 'none', 'facealpha', .2);
            hold on;
            xnew = expert(:, expert_ses == 2);
            plot(nanmean(xnew, 2), 'LineWidth', 3, 'LineStyle', ':', 'Color', c(1:3))
            hi = nanmean(xnew, 2) + 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); lo = nanmean(xnew, 2) - 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); x = (1:size(nanmean(xnew, 2),1))';
            hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], c(1:3));
            set(hp1, 'facecolor', c(1:3), 'edgecolor', 'none', 'facealpha', .2);
            g = gca;
            g.YLabel.String = 'Fractional Anisotropy (FA)';
            g.YLabel.FontSize = fontsize;
            
            % xaxis
            xax = get(gca, 'xaxis');
            xax.Limits = [0 160];
            xax.TickValues = [0 80 160];
            xax.TickLabels = {'20', '100', '180'};
            xax.TickDirection = 'out';
            xax.TickLength = [xticklength xticklength];
            xax.FontName = fontname;
            xax.FontSize = fontsize;
            xax.FontAngle = fontangle;
            
            % yaxis
            yax = get(gca,'yaxis');
            yax.Limits = [ylim_lo ylim_hi];
            yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
            yax.TickDirection = 'out';
            yax.TickLabels = {num2str(ylim_lo, '%2.2f'), num2str((ylim_lo+ylim_hi)/2, '%2.2f'), num2str(ylim_hi, '%2.2f')};
            yax.FontName = fontname;
            yax.FontSize = fontsize;
            
            % general
            g = gca;
            box off
            g.XLabel.String = 'Location along tract';
            g.XLabel.FontSize = fontsize;
            g.XLabel.FontAngle = fontangle;
            pbaspect([1 1 1])
            title('Experts')
            hold off;
            
            subplot(1, 3, 2)
            xnew = beg(:, beg_ses == 1); c = beg_color;
            plot(nanmean(xnew, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', c(1:3))
            hi = nanmean(xnew, 2) + 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); lo = nanmean(xnew, 2) - 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); x = (1:size(nanmean(xnew, 2),1))';
            hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], c(1:3));
            set(hp1, 'facecolor', c(1:3), 'edgecolor', 'none', 'facealpha', .2);
            hold on;
            xnew = beg(:, beg_ses == 2);
            plot(nanmean(xnew, 2), 'LineWidth', 3, 'LineStyle', ':', 'Color', c(1:3))
            hi = nanmean(xnew, 2) + 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); lo = nanmean(xnew, 2) - 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); x = (1:size(nanmean(xnew, 2),1))';
            hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], c(1:3));
            set(hp1, 'facecolor', c(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            % xaxis
            xax = get(gca, 'xaxis');
            xax.Limits = [0 160];
            xax.TickValues = [0 80 160];
            xax.TickLabels = {'20', '100', '180'};
            xax.TickDirection = 'out';
            xax.TickLength = [xticklength xticklength];
            xax.FontName = fontname;
            xax.FontSize = fontsize;
            xax.FontAngle = fontangle;
            
            % yaxis
            yax = get(gca,'yaxis');
            yax.Limits = [ylim_lo ylim_hi];
            yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
            yax.TickDirection = 'out';
            yax.TickLabels = {num2str(ylim_lo, '%2.2f'), num2str((ylim_lo+ylim_hi)/2, '%2.2f'), num2str(ylim_hi, '%2.2f')};
            yax.FontName = fontname;
            yax.FontSize = fontsize;
            
            % general
            g = gca;
            box off
            g.XLabel.String = 'Location along tract';
            g.XLabel.FontSize = fontsize;
            g.XLabel.FontAngle = fontangle;
            pbaspect([1 1 1])
            title('Beginners')
            hold off;
            
            subplot(1, 3, 3)
            xnew = con(:, con_ses == 1); c = con_color;
            plot(nanmean(xnew, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', c(1:3))
            hi = nanmean(xnew, 2) + 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); lo = nanmean(xnew, 2) - 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); x = (1:size(nanmean(xnew, 2),1))';
            hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], c(1:3));
            set(hp1, 'facecolor', c(1:3), 'edgecolor', 'none', 'facealpha', .2);
            hold on;
            xnew = expert(:, con_ses == 2);
            plot(nanmean(xnew, 2), 'LineWidth', 3, 'LineStyle', ':', 'Color', c(1:3))
            hi = nanmean(xnew, 2) + 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); lo = nanmean(xnew, 2) - 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); x = (1:size(nanmean(xnew, 2),1))';
            hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], c(1:3));
            set(hp1, 'facecolor', c(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            % xaxis
            xax = get(gca, 'xaxis');
            xax.Limits = [0 160];
            xax.TickValues = [0 80 160];
            xax.TickLabels = {'20', '100', '180'};
            xax.TickDirection = 'out';
            xax.TickLength = [xticklength xticklength];
            xax.FontName = fontname;
            xax.FontSize = fontsize;
            xax.FontAngle = fontangle;
            
            % yaxis
            yax = get(gca,'yaxis');
            yax.Limits = [ylim_lo ylim_hi];
            yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
            yax.TickDirection = 'out';
            yax.TickLabels = {num2str(ylim_lo, '%2.2f'), num2str((ylim_lo+ylim_hi)/2, '%2.2f'), num2str(ylim_hi, '%2.2f')};
            yax.FontName = fontname;
            yax.FontSize = fontsize;
            
            % general
            g = gca;
            box off
            g.XLabel.String = 'Location along tract';
            g.XLabel.FontSize = fontsize;
            g.XLabel.FontAngle = fontangle;
            pbaspect([1 1 1])
            title('Controls')
            hold off;
            
            sgtitle(list_tract{k})
            
            print(fullfile(rootDir, 'plots', ['plot_tractprofiles_' wm_measure '_' list_tract{k}]), '-dpng')
            print(fullfile(rootDir, 'plots', 'eps', ['plot_tractprofiles_' wm_measure '_' list_tract{k}]), '-depsc')
            
            % Open a new figure for the mean plot.
            figcount = figcount + 1;
            figure(figcount)
            
            xnew1 = expert(:, expert_ses == 1); xnew2 = expert(:, expert_ses == 2); c = exp_color;
            if size(xnew1, 2) == size(xnew2, 2) % must have both session 1 and session 2 for all subs
                xnew = xnew2 - xnew1;
                plot(nanmean(xnew, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', c(1:3))
                hi = nanmean(xnew, 2) + 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); lo = nanmean(xnew, 2) - 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); x = (1:size(nanmean(xnew, 2),1))';
                hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], c(1:3));
                set(hp1, 'facecolor', c(1:3), 'edgecolor', 'none', 'facealpha', .2);
                hold on;
            elseif size(xnew1, 2) < size(xnew2, 2)
                disp(['Unable to plot session difference for ' wm_measure ' for ' list_tract{k} ' because at least one subject is missing tractprofile for session 1.'])
            else size(xnew1, 2) > size(xnew2, 2)
                disp(['Unable to plot session difference for ' wm_measure ' for ' list_tract{k} ' because at least one subject is missing tractprofile for session 2.'])
            end
            xnew1 = beg(:, beg_ses == 1); xnew2 = beg(:, beg_ses == 2); c = beg_color;
            if size(xnew1, 2) == size(xnew2, 2) % must have both session 1 and session 2 for all subs
                xnew = xnew2 - xnew1;
                plot(nanmean(xnew, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', c(1:3))
                hi = nanmean(xnew, 2) + 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); lo = nanmean(xnew, 2) - 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); x = (1:size(nanmean(xnew, 2),1))';
                hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], c(1:3));
                set(hp1, 'facecolor', c(1:3), 'edgecolor', 'none', 'facealpha', .2);
                hold on;
            elseif size(xnew1, 2) < size(xnew2, 2)
                disp(['Unable to plot session difference for ' wm_measure ' for ' list_tract{k} ' because at least one subject is missing tractprofile for session 1.'])
            else size(xnew1, 2) > size(xnew2, 2)
                disp(['Unable to plot session difference for ' wm_measure ' for ' list_tract{k} ' because at least one subject is missing tractprofile for session 2.'])
            end
            xnew1 = con(:, con_ses == 1); xnew2 = con(:, con_ses == 2); c = con_color;
            if size(xnew1, 2) == size(xnew2, 2) % must have both session 1 and session 2 for all subs
                xnew = xnew2 - xnew1;
                plot(nanmean(xnew, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', c(1:3))
                hi = nanmean(xnew, 2) + 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); lo = nanmean(xnew, 2) - 1.96*nanstd(xnew, 0, 2)/sqrt(size(~isnan(xnew), 2)); x = (1:size(nanmean(xnew, 2),1))';
                hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], c(1:3));
                set(hp1, 'facecolor', c(1:3), 'edgecolor', 'none', 'facealpha', .2);
                hold on;
            elseif size(xnew1, 2) < size(xnew2, 2)
                disp(['Unable to plot session difference for ' wm_measure ' for ' list_tract{k} ' because at least one subject is missing tractprofile for session 1.'])
            else size(xnew1, 2) > size(xnew2, 2)
                disp(['Unable to plot session difference for ' wm_measure ' for ' list_tract{k} ' because at least one subject is missing tractprofile for session 2.'])
            end
            
            % Add reference line at y=0.
            plot([0 180], [0 0], 'k')
            
            % xaxis
            xax = get(gca, 'xaxis');
            xax.Limits = [0 160];
            xax.TickValues = [0 80 160];
            xax.TickLabels = {'20', '100', '180'};
            xax.TickDirection = 'out';
            xax.TickLength = [xticklength xticklength];
            xax.FontName = fontname;
            xax.FontSize = fontsize;
            xax.FontAngle = fontangle;
            
            % yaxis
            yax = get(gca,'yaxis');
            yax.Limits = [ylim_diff_lo ylim_diff_hi];
            yax.TickValues = [ylim_diff_lo (ylim_diff_lo+ylim_diff_hi)/2 ylim_diff_hi];
            yax.TickDirection = 'out';
            yax.TickLabels = {num2str(ylim_diff_lo, '%2.2f'), num2str((ylim_diff_lo+ylim_diff_hi)/2, '%2.2f'), num2str(ylim_diff_hi, '%2.2f')};
            yax.FontName = fontname;
            yax.FontSize = fontsize;
            
            % general
            g = gca;
            box off
            g.XLabel.String = 'Location along tract';
            g.XLabel.FontSize = fontsize;
            g.XLabel.FontAngle = fontangle;
            g.YLabel.String = {'Difference in Fractional Anisotropy (FA)'; '(Session 2 - Session 1)'};
            g.YLabel.FontSize = fontsize;
            
            pbaspect([1 1 1])
            title(list_tract{k})
            
            print(fullfile(rootDir, 'plots', ['plot_tractprofiles_diff_' wm_measure '_' list_tract{k}]), '-dpng')
            print(fullfile(rootDir, 'plots', 'eps', ['plot_tractprofiles_diff_' wm_measure '_' list_tract{k}]), '-depsc')
            
            hold off;
            
        end % if toi
        
    end %tract
    
end %w




