% This script reads in SNR values and plots them according to session (pre-
% vs post-training) and group (expert=3, beginner=2, control=1).

clear all; close all; clc
format shortG

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
            
            % Get SNR for this subject.
            data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
            
            % Get SNR in b0 images.
            b0(sub_count) = str2num(data_snr_temp.SNRInB0_X_Y_Z{1});
            
            % Get mean SNR in X, Y, and Z directions.
            m(sub_count) = mean([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get standard deviation of SNR in X, Y, and Z directions.
            sd(sub_count) = std([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get subID.
            subID(sub_count) = str2num(grp_contents(s).name(5:7));
            
            % Get session.
            ses(sub_count) = str2num(grp_contents(s).name(end));
            
            % Get ylabel.
            lab{sub_count} = grp_contents(s).name;
            
            % Get training group.
            group(sub_count) = beh_data_in_tbl.DanceLevelCode(find((beh_data_in_tbl.No == str2num(grp_contents(s).name(5:7)))));
            
            clear data_snr_temp get_temp
            
        end % if ~isempty
        
    end % end if exist
    
end % end s

% SNR plot
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

c_color = [0 0.4470 0.7410];
b_color = [0.4660 0.6740 0.1880];
e_color = [0.6350 0.0780 0.1840];

gscatter(b0, 1:length(b0), group, [c_color; b_color; e_color], 'x', 10)
hold on;
for p = 1:length(m)
    
    if group(p) == 1 && ses(p) == 1
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', c_color)
        
    elseif group(p) == 1 && ses(p) == 2
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', c_color, 'LineStyle', ':')
        
    elseif group(p) == 2 && ses(p) == 1
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', b_color)
        
    elseif group(p) == 2 && ses(p) == 2
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', b_color, 'LineStyle', ':')
        
    elseif group(p) == 3 && ses(p) == 1
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', e_color)
        
    elseif group(p) == 3 && ses(p) == 2
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', e_color, 'LineStyle', ':')
        
    end
    
end
gscatter(m, 1:length(m), group, [c_color; b_color; e_color], '.', 30)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [min([b0 m + sd])-5 max([b0 m + sd])+5];
xax.TickValues = floor(min([b0 m + sd])):10:ceil(max([b0 m + sd]));
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(lab)+0.5];
yax.TickValues = 1:length(lab);
yax.TickDirection = 'out';
yax.TickLabels = lab;
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

legend({'Controls', 'Beginners', 'Experts', 'Session 1', 'Session 2'}, 'Location', 'northeast');
legend box off

a.XLabel.String = 'SNR';
a.XLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_snr'), '-dpng')

hold off;