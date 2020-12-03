% This script reads in FA and MD measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App).
% It also reads in behavioral data (e.g., experimental group) collected as part of the
% spade study.
%it runs twice, one for producing the plots for mena Fa and difference (S2-S1) and one for mean MD and
%Md difference (S2-S1)

clear all; close all; clc
format shortG

% Set working directories.
% rootDir = '/Users/EPL_VR_02/OneDrive - University of Cyprus/PhD_SPADE/SPADE_BLData/';
rootDir = '/Volumes/240/spade/';

% Get bl project foldername.
% blprojectid = 'SPADE_TractProfiles';
blprojectid = 'proj-5e61139282b37f2cfe8fdb28';

w_measures = {'fa', 'md'};

exp_color = [0.6350 0.0780 0.1840]; %red
beg_color = [0 0.4470 0.7410]; %blue
con_color = [0.41176 0.41176 0.41176]; %gray

hold on;
linewidth = 1.5;
linestyle = 'none';
fontname = 'Times New Roman';
fontsize = 10;
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
for w = 1%:length(w_measures)

    wm_measure = w_measures{w};

    if strcmp(wm_measure, 'fa')
        ylim_lo = 0.20; ylim_hi = 0.70; %start the yaxis numbering from 0.20 to 0.70. the FA has a range from 0-1.
        ylim_diff_lo = -0.25; ylim_diff_hi = 0.25; % do a step on the axis every .25.
        ylabel = 'Fractional Anisotropy (FA)';
        ylabel_diff = 'Difference in Fractional Anisotropy (FA)';
    elseif strcmp(wm_measure, 'md')
        ylim_lo = 0.10; ylim_hi = 1.20;
        ylim_diff_lo = -0.25; ylim_diff_hi = 0.25;
        ylab = 'Mean Diffusivity (MD)';
        ylab_diff = 'Difference in Mean Diffusivity (MD)';
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

        % Only collect values for subjects that have both MRI and behavioral/demographic data.
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
                if i == 1 && j == 1 %logical indexing, means all true

                    tract = {};

                end

                % Read in data for this subject and this tract.
                data_temp = readtable(fullfile(sub_contents_tractprofiles(j).folder, sub_contents_tractprofiles(j).name));

                % Get middle 80%.
                start = size(data_temp, 1)*.1;
                stop = size(data_temp, 1)*.9;

                % Read in mean WM measure.
                if strcmp(wm_measure, 'fa')

                    m_wm(:, j, sub_count) = data_temp.fa_mean(start:stop);
                    sd_wm(:, j, sub_count) = data_temp.fa_sd(start:stop);

                elseif strcmp(wm_measure, 'md')

                    m_wm(:, j, sub_count) = data_temp.md_mean(start:stop);
                    sd_wm(:, j, sub_count) = data_temp.md_sd(start:stop);

                end

                % Grab tract name for grouping variable.
                %note: repmat fuction stands for 'Repeat copies of array'

                tract(:, j, sub_count) = repmat({sub_contents_tractprofiles(j).name(1:end-13)}, 161, 1);

                % Grab subID.
                sub(:, j, sub_count) = repmat(str2num(grp_contents(i).name(5:7)), 161, 1);

                % Gather session, for ease.
                ses(:, j, sub_count) = repmat(str2num(grp_contents(i).name(end)), 161, 1);

                % Get exp group.
                group(sub_count) = beh_data_in_tbl.DanceLevelCode(find((beh_data_in_tbl.No == str2num(grp_contents(i).name(5:7)))));

                % Get age in months.
                age(sub_count) = beh_data_in_tbl.Age(find((beh_data_in_tbl.No == str2num(grp_contents(i).name(5:7)))));

                %clear data_temp

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
    subID = subID(subID ~= 0);

% Select only tracts that we care about.
t_count = 0;
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
            || strcmp(list_tract{k}, 'leftCST') || strcmp(list_tract{k}, 'rightCST') ...
            && ~strcmp(list_tract{k}, 'empty')
        
        % keep tract name for later
        t_count = t_count + 1;
        tractname{t_count} = list_tract{k};
        
        % Find entries that are for this tract.
        %note: t_idx stands for tract index
        t_idx = strcmp(tract, list_tract{k});
        
        %             % Open a new figure for this tract.
        %             figcount = figcount + 1;
        %             figure(figcount)
        
        count = 0; exp_count = 0; beg_count = 0; con_count = 0;
        for s = 1:length(subID)
            
            % Only include subjects who are not outliers.
            if ~ismember(subID(s), outlier)
                
                % Find entries that are for this subject.%note:s_idx
                % stands for subject index
                s_idx = sub == subID(s);
                
                for session = 1:2
                    
                    % Find entries that are for this session.
                    ses_idx = ses == session;
                    
                    % Subset the thing so that we only plot for this tract, subject, and session.
                    %note: t_idx ==1 means that its the index for this
                    %tract when this is true equals to 1 (logical) and it changes
                    %everytime based on the loop eg., for session 2, sub 114, rightVOF
                    
                    t_temp = m_wm(t_idx == 1 & s_idx == 1 & ses_idx == 1);
                    
                    if isempty(t_temp)
                        
                        t_temp = NaN(161, 1);
                        
                    end
                    
                    count = count + 1;
                    %
                    %                         % Different line styles for session 1 and session 2.
                    %                         if session == 1
                    %                             linestyle = '-';
                    %                         elseif session == 2
                    %                             linestyle = ':';
                    %                         end
                    
                    % Code the plot for subject and keep data for inspection (yc, oc, a).
                    if group(count) == 3 % expert
                        
                        exp_count = exp_count + 1;
                        
                        % Collect.
                        expert(:, exp_count) = t_temp;
                        expert_ses(exp_count) = session;
                        
                        
                    elseif group(count) == 2 % beginner
                        
                        beg_count = beg_count + 1;
                        
                        % Collect.
                        beg(:, beg_count) = t_temp;
                        beg_ses(beg_count) = session;
                        
                    elseif group(count) == 1 % control
                        
                        con_count = con_count + 1;
                        
                        % Collect.
                        con(:, con_count) = t_temp;
                        con_ses(con_count) = session;
                        
                    end
                    
                    %                         end % if ~isempty
                    
                    clear t_temp;
                    
                end % for session
                
            end % if ~ismember(subID(s), outlier)
            
        end % for s
        
        %             figcount = figcount + 1;
        %             figure(figcount)
        %
        %             disp(list_tract{k})
        %note: 'xnew1' variable is the mean for session 1 and 'xnew2'
        %variable is the mean for session 2', 'xnew' is the difference
        %between the two sessions
        
        xnew1 = expert(:, expert_ses == 1); xnew2 = expert(:, expert_ses == 2);
        xnew = xnew2 - xnew1;
        %         s1_exp(:, k) = mean(xnew1, 1)'; s2_exp(:, k) = mean(xnew2, 1)';
        tdiff_exp(:, k) = mean(xnew, 1)';
        
        xnew1 = beg(:, beg_ses == 1); xnew2 = beg(:, beg_ses == 2);
        xnew = xnew2 - xnew1;
        %         s1_beg(:, k) = mean(xnew1, 1)'; s2_beg(:, k) = mean(xnew2, 1)';
        tdiff_beg(:, k) = mean(xnew, 1)';
        
        xnew1 = con(:, con_ses == 1); xnew2 = con(:, con_ses == 2);
        xnew = xnew2 - xnew1;
        xnew = xnew(:, 1:end-1); % because I think the final con subject is bad
        %         s1_con(:, k) = mean(xnew1, 1)'; s2_con(:, k) = mean(xnew2, 1)';
        tdiff_con(:, k) = mean(xnew, 1)';
        
    end % if toi
    
end %tract
%
% end %w

% Come back here -- this is a potential source of error in the future.
tdiff_exp( :, ~any(tdiff_exp,1) ) = [];  %columns
tdiff_beg( :, ~any(tdiff_beg,1) ) = [];  %columns
tdiff_con( :, ~any(tdiff_con,1) ) = [];  %columns

diff_exp = nanmean(tdiff_exp, 1);
diff_beg = nanmean(tdiff_beg, 1);
diff_con = nanmean(tdiff_con, 1);

figure(1)

subplot(3, 1, 1)
% clear tidx
[~, tidx2] = sort(diff_beg, 'descend');

imagesc(cat(1, diff_exp(tidx2), diff_beg(tidx2), diff_con(tidx2)));

% xaxis
xax = get(gca, 'xaxis');
%         xax.Limits = [0.5 3.5];
xax.TickValues = 1:size(tdiff_exp, 2);
xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
xlabels = tractname(tidx2);
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
%     xax.FontName = fontname;
%     xax.FontSize = fontsize;
xax.TickLabelRotation = 45;

yax = get(gca, 'yaxis');
yax.Limits = [0.5 3.5];
yax.TickValues = [1 2 3];
yax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
ylabels = {'Experienced' ; 'Beginner'; 'Controls'};
ylabels = cellfun(@(x) strrep(x, ',', '\newline'), ylabels, 'UniformOutput', false);
yax.TickLabels = ylabels;

colorbar; caxis([-.005, .005]);
%     xax.FontName = fontname;
%     xax.FontSize = fontsize;
%

title('Sorted for Beginner')
pbaspect([5 1 1])

subplot(3, 1, 2)
% clear tidx2;
[~, tidx3] = sort(diff_exp, 'descend');

imagesc(cat(1, diff_exp(tidx3), diff_beg(tidx3), diff_con(tidx3)));

% xaxis
xax = get(gca, 'xaxis');
%         xax.Limits = [0.5 3.5];
xax.TickValues = 1:size(tdiff_exp, 2);
xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
xlabels = tractname(tidx3);
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
%     xax.FontName = fontname;
%     xax.FontSize = fontsize;
xax.TickLabelRotation = 45;

yax = get(gca, 'yaxis');
yax.Limits = [0.5 3.5];
yax.TickValues = [1 2 3];
yax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
ylabels = {'Experienced' ; 'Beginner'; 'Controls'};
ylabels = cellfun(@(x) strrep(x, ',', '\newline'), ylabels, 'UniformOutput', false);
yax.TickLabels = ylabels;

colorbar; caxis([-.005, .005]);
%     xax.FontName = fontname;
%     xax.FontSize = fontsize;
%

title('Sorted for Expert')
pbaspect([5 1 1])

subplot(3, 1, 3)
% clear tidx2;
[~, tidx4] = sort(diff_con, 'descend');

imagesc(cat(1, diff_exp(tidx4), diff_beg(tidx4), diff_con(tidx4)));

% xaxis
xax = get(gca, 'xaxis');
%         xax.Limits = [0.5 3.5];
xax.TickValues = 1:size(tdiff_exp, 2);
xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
xlabels = tractname(tidx4);
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
%     xax.FontName = fontname;
%     xax.FontSize = fontsize;
xax.TickLabelRotation = 45;

yax = get(gca, 'yaxis');
yax.Limits = [0.5 3.5];
yax.TickValues = [1 2 3];
yax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
ylabels = {'Experienced' ; 'Beginner'; 'Controls'};
ylabels = cellfun(@(x) strrep(x, ',', '\newline'), ylabels, 'UniformOutput', false);
yax.TickLabels = ylabels;

colorbar; caxis([-.005, .005]);
%     xax.FontName = fontname;
%     xax.FontSize = fontsize;
%

title('Sorted for Control')
pbaspect([5 1 1])

print(fullfile(rootDir, 'plots', ['plot_diffmatrix_' wm_measure ]), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', ['plot_diffmatrix_' wm_measure]), '-depsc')

end