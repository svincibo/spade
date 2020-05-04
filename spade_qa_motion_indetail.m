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
            lab{sub_count} = grp_contents(s).name(5:7);
            
            % Get training group.
            group(sub_count) = beh_data_in_tbl.DanceLevelCode(find((beh_data_in_tbl.No == str2num(grp_contents(s).name(5:7)))));
            
        end % if ~isempty
        
    end % end if exist
    
end % end s

meanmotion = mean(motion, 2);
sd = std(motion, 0, 2);

% FD plot
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
xticklength = 0;
alphablend = .8;

gscatter(meanmotion(ses == 1), 1:2:2*length(meanmotion(ses == 1)), group(ses == 1), [control_color; beginner_color; expert_color], '.', 30)

hold on;
for p = 1:length(meanmotion)
    
    if group(p) == 1 && ses(p) == 1
        
        plot([meanmotion(p) - abs(sd(p)) meanmotion(p) + abs(sd(p))], [p p], 'Color', control_color)
        
    elseif group(p) == 1 && ses(p) == 2
        
        plot([meanmotion(p) - abs(sd(p)) meanmotion(p) + abs(sd(p))], [p p], 'Color', control_color, 'LineStyle', ':')
        
    elseif group(p) == 2 && ses(p) == 1
        
        plot([meanmotion(p) - abs(sd(p)) meanmotion(p) + abs(sd(p))], [p p], 'Color', beginner_color)
        
    elseif group(p) == 2 && ses(p) == 2
        
        plot([meanmotion(p) - abs(sd(p)) meanmotion(p) + abs(sd(p))], [p p], 'Color', beginner_color, 'LineStyle', ':')
        
    elseif group(p) == 3 && ses(p) == 1
        
        plot([meanmotion(p) - abs(sd(p)) meanmotion(p) + abs(sd(p))], [p p], 'Color', expert_color)
        
    elseif group(p) == 3 && ses(p) == 2
        
        plot([meanmotion(p) - abs(sd(p)) meanmotion(p) + abs(sd(p))], [p p], 'Color', expert_color, 'LineStyle', ':')
        
    end
    
end
gscatter(meanmotion(ses == 2), (1:2:2*length(meanmotion(ses == 2)))+1, group(ses == 2), [control_color; beginner_color; expert_color], 'o', 8)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 1];
xax.TickValues = 0:0.25:1.25;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(lab)+0.5];
yax.TickValues = 1.5:2:ceil(length(lab));
yax.TickDirection = 'out';
yax.TickLabels = lab(ses == 1);
yax.FontName = fontname;
yax.FontSize = 8;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

legend({'Controls', 'Beginners', 'Experts', 'Session 1', 'Session 2'}, 'Location', 'northeast');
legend box off

a.XLabel.String = 'Framewise Displacement (FD)';
a.XLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_fd'), '-dpng')

hold off;