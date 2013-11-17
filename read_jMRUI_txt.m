function [spectra, Peak_ratios, Metabolites] = MRS_read_jMRUI_txt_output_v3()
% Load jMRUI output (.txt) data and interpret the format
% of the saved data. If multiple acquisitions were acquired, it will detect
% values outside a specified range and provide a sorted list of spectra
% based upon the number of outliers. This is a simple way of identifying
% artifacts in the spectra.
%
% If the multiple spectra come from a Diff acquisition, then the peaks are
% identified within the specified range.

% igreenhouse@berkeley.edu August, 2013
clear all
%% Edit These Parameters

% xrange = 1300:2000; % 1300:2000 is specific to 1536 zero fill from UC Berkeley DICOM
% xrange = 1600:2300; % 1600:2300 is specific to 2048 zero fill from UC Berkeley DICOM
xrange = 1325:1925; % 1300:2000 is specific to 2048 zero fill from UC Berkeley RDA

std_thresh = 2; % number of std for identifying outliers

NAA_peak = 1384; % index for NAA peak (approximate is fine)

dif_minimum_peak_height = 300; % 300 based upon pilot testing at UC Berkeley. this may need to be changed for frontal voxels (try 200)

xlims = [1200 2200]; % x-values for outliers plots

%% Select txt file
[filename, filepath] = uigetfile('.txt');
%% Read Header
fid = fopen(fullfile(filepath, filename),'r');  % Open text file


InputText=textscan(fid,'%s',19,'delimiter','\n');  % Read strings delimited
                                                  % by a carriage return                                                  
Header=InputText{1};
disp(Header);
Columns = textscan(fid, '%s',4, 'delimiter','\t');
ColumnNames = Columns{1};
disp(ColumnNames);
skip_lines = textscan(fid,'%s', 2, 'delimiter','\n');
junk = skip_lines{1};
disp(junk);
%% Read Data
measure = 1;
while (~feof(fid))
    C = textscan(fid,'%n%n%n%n');    
    data(measure,:) = C{3};
    skip_lines = textscan(fid,'%s', 1, 'delimiter','\n');
    junk = skip_lines{1};
    disp(junk);
    measure = measure+1;    
end

data = data';

%% Determine Vector Size and Number of Spectra
vector_size = size(data,1); % vector size
num_spectra = size(data,2); % Number of Spectra

switch num_spectra; 
    case 3, 
        disp('3 Spectra: probably On, Off, Diff'), 
    otherwise
        disp(sprintf('%d spectra', num_spectra));
end

%% categorize spectra
a = 1; b = 1; c = 1;
for i = 1:size(data,2),
    if data(NAA_peak,i)<-2000,
        spectra.diff(:,a) = data(:,i);
        a = a + 1;
    elseif data(NAA_peak,i)>2000,
        spectra.off(:,b) = data(:,i);
        b = b + 1;
    else
        spectra.on(:,c) = data(:,i);
        c = c + 1;
    end
end
%% Plot all data if single measure from scanner

switch num_spectra,
    case 3,

        figure1 = figure();
%         on_off_diff_spectra = figure(1);
        subplot(4,1,1);
        plot([spectra.diff(xrange) spectra.off(xrange), spectra.on(xrange)]);
        title('All');
        legend('Diff', 'Off', 'On');
        
        set(gca, 'FontWeight', 'Bold');
        xlim([0 600]);

        subplot(4,1,2);
        on_and_off_spectra = spectra.on + spectra.off;
        plot(on_and_off_spectra(xrange));
        title('Summed On & Off Spectra');
        set(gca, 'FontWeight', 'Bold');
        xlim([0 600]);

        subplot(4,1,3);
        plot(spectra.diff(xrange));
        title('Diff');
        set(gca, 'FontWeight', 'Bold');
        xlim([0 600]);
    
    
    %% Identify peaks in Diff spectra
        mean_diff = spectra.diff(xrange);
        [pks, locs] = findpeaks(mean_diff, 'minpeakheight', dif_minimum_peak_height);                  
        for i = 1:length(locs),
            plotline(locs(i),'v');
        end
        if pks,
            disp('Peaks have been identified.');
        end            
    
        sorted_num_of_outliers = 'Single On, Off, Diff Spectra';

%% Review data acquired with multiple measures
    
    otherwise
        figure1 = figure(1);        
        
%         spectra.stdev_diff = std(spectra.diff'); % calculate std of diff
%         spectra % diff isn't needed for outlier detection since it is
%         calculated from the Off and On resonance samples
        spectra.stdev_off = std(spectra.off'); % calculate std of off spectra
        spectra.stdev_on = std(spectra.on'); % calculate std of off spectra

        % preallocate variables
%         outliers(size(data,1), size(data,2)) = zeros;
%         diff_outliers(size(spectra.diff,1), size(spectra.diff,2)) = zeros;
        
        spectra.on_outliers(size(spectra.on,1), size(spectra.on,2)) = zeros;
        spectra.off_outliers(size(spectra.off,1), size(spectra.off,2)) = zeros;

        
%         number_of_outliers(1:size(data,2),1) = 1:size(data,2); % number of spectral measures
%         number_diff_outliers(1:size(spectra.diff,2),1) = 1:size(spectra.diff,2); % number of diff measures
        spectra.number_on_outliers(1:size(spectra.on,2),1) = 1:size(spectra.on,2); % number of diff measures        
        spectra.number_off_outliers(1:size(spectra.off,2),1) = 1:size(spectra.off,2); % number of diff measures
        
        for sample = 1:size(spectra.on,2),    % step through each sample
            for i = 1:length(spectra.on),     % step through each datapoint
                if abs(spectra.on(i,sample)-mean(spectra.on(i,:))) > std_thresh*spectra.stdev_on(i); % identify datapoints outside threshold
                    spectra.on_outliers(i,sample) = spectra.on(i,sample); % record outlying data
                end
                if abs(spectra.off(i,sample)-mean(spectra.off(i,:))) > std_thresh*spectra.stdev_off(i); % identify datapoints outside threshold
                    spectra.off_outliers(i,sample) = spectra.off(i,sample); % record outlying data
                end
            end
            spectra.number_on_outliers(sample,2) = length(find(spectra.on_outliers(:,sample))); % record number of outliers for each acquisition
            spectra.number_off_outliers(sample,2) = length(find(spectra.off_outliers(:,sample))); % record number of outliers for each acquisition
        end

        % Sort according to number of statistical outliers
        spectra.sorted_num_off_outliers = sortrows(spectra.number_off_outliers,-2);
        spectra.sorted_num_on_outliers = sortrows(spectra.number_on_outliers,-2);
%% Plot Spectra        
        subplot(2,1,1);
        plot(spectra.off(:,spectra.sorted_num_off_outliers(:,1)));
        legend(num2str(spectra.sorted_num_off_outliers(:,1)));
        title_text = sprintf('Off: Acquisitions Ordered Most -> Least Variable About Mean');
        Title(title_text);
        xlim(xlims);        
        
        subplot(2,1,2);
        plot(spectra.on(:,spectra.sorted_num_on_outliers(:,1)));
        legend(num2str(spectra.sorted_num_on_outliers(:,1)));
        title_text = sprintf('On: Acquisitions Ordered Most -> Least Variable About Mean');
        Title(title_text);
        xlim(xlims);
        
        disp('Visually inspect spectra to identify outliers, and enter any spectra you wish to ignore in the analysis.');
%% Input Spectra to Trimming
        trim_off_spectra = input('Enter number(s) of Off spectra to trim (e.g. [1 3 7], else 0): ');
        trim_on_spectra = input('Enter number(s) of On spectra to trim (e.g. [1 3 7], else 0): ');
    
%% Create indices for Trimming
        for i = 1:length(spectra.sorted_num_off_outliers),
            if find(trim_off_spectra==spectra.sorted_num_off_outliers(i,1))
                spectra.trim_off(i,1) = 0;
            else spectra.trim_off(i,1) = 1;
            end
        end
        for i = 1:length(spectra.sorted_num_on_outliers),
            if find(trim_on_spectra==spectra.sorted_num_on_outliers(i,1))
                spectra.trim_on(i,1) = 0;
            else spectra.trim_on(i,1) = 1;
            end
        end
%% Plot Trimmed Spectra       
        figure2 = figure(2);
        subplot(2,1,1);
        plot(spectra.off(:,spectra.sorted_num_off_outliers(:,1) & spectra.trim_off));
        legend(num2str(spectra.sorted_num_off_outliers(find(spectra.trim_off))));
        title_text = sprintf('Off: Trimmed ');
        Title(title_text);
        xlim(xlims);        
        
        subplot(2,1,2);
        plot(spectra.on(:,spectra.sorted_num_on_outliers(:,1) & spectra.trim_on));
        legend(num2str(spectra.sorted_num_on_outliers(find(spectra.trim_on))));
        title_text = sprintf('On: Acquisitions Ordered Most -> Least Variable About Mean');
        Title(title_text);
        xlim(xlims);
        
        disp('Compare with earlier plots to determine if outliers were properly removed. Hit any key to proceed.');
        pause

%% Calculate Mean Difference and Summed On & Off from Trimmed Spectra
        spectra.trimmed_mean_off = mean(spectra.off(:,spectra.sorted_num_off_outliers(:,1) & spectra.trim_off),2);
        spectra.trimmed_mean_on = mean(spectra.on(:,spectra.sorted_num_on_outliers(:,1) & spectra.trim_on),2);
        spectra.trimmed_mean_diff = spectra.trimmed_mean_on - spectra.trimmed_mean_off;
%         spectra.trimmed_sum_On_w_Off = spectra.trimmed_mean_on + spectra.trimmed_mean_off;
         spectra.trimmed_sum_On_w_Off = mean([spectra.trimmed_mean_off'; spectra.trimmed_mean_on'])'; %possibly use average instead of sum;

%%      Surface Plot of mean +/- stdev
%         xlims = [0 600]; 
%         figure3 = figure();
%         errorSurface(xrange-(xrange(1)-1),mean(data(xrange,:),2)', stdev_data(xrange)');
%         hold
%         plot(mean(data(xrange,:)')','k');
%         title_text = sprintf('%s Mean (+/- Stdev)', spectra_type);
%         Title(title_text);
%         xlim(xlims);
%         worst_sample = sorted_num_of_outliers(1,1);
%         plot(data(xrange,worst_sample),'r');

%% Plot and Identify peaks in Diff Spectrum
        figure3 = figure(3);
        
        plot(spectra.trimmed_mean_diff(xrange));
        hold
        
        title_text = 'Diff: Mean';
        Title(title_text);
        set(gca, 'FontWeight', 'Bold');
        xlim([0 600]);
        
        mean_diff = spectra.trimmed_mean_diff(xrange);
        [pks, locs] = findpeaks(mean_diff, 'minpeakheight', dif_minimum_peak_height);                  
        for i = 1:length(locs),
            plotline(locs(i),'v');
        end
        if pks,
            disp('Peaks have been identified.');
        end
end
%% Calculate Glx reference point and adjusted xrange for data       
    Glx_peak1 = locs(end);
    Glx_peak2 = locs(end-1);
    Glx_peak_separation = Glx_peak1-Glx_peak2;
    half_Glx_peak_separation = round(Glx_peak_separation/2);
    if mod(Glx_peak_separation,2),
        Glx_reference_point = Glx_peak1 - half_Glx_peak_separation;
    elseif ((Glx_peak1 - half_Glx_peak_separation)~=(Glx_peak2 + (half_Glx_peak_separation-1))) &...
            mean_diff(Glx_peak1 - half_Glx_peak_separation) < ...
            mean_diff(Glx_peak2 + (half_Glx_peak_separation-1));
        Glx_reference_point = Glx_peak1 - half_Glx_peak_separation;
    elseif ((Glx_peak1 - half_Glx_peak_separation)~=(Glx_peak2 + (half_Glx_peak_separation-1))) &...
            mean_diff(Glx_peak1 - half_Glx_peak_separation) > ...
            mean_diff(Glx_peak2 + (half_Glx_peak_separation-1));
        Glx_reference_point = Glx_peak2 + half_Glx_peak_separation;
    end

    limits_ref_to_Glx = [Glx_reference_point - 500 Glx_reference_point + 100];
    adjusted_xrange = (xrange(1) + limits_ref_to_Glx(1)):(xrange(1) + limits_ref_to_Glx(2));
%% Plot and Calculate GABA and Glx Area Under Peaks from Difference Spectra
    % NOTE: These ranges were derived from the UC Davis Template
    % provided by Rick Maddock (rjmaddock@ucdavis.edu)
 if num_spectra > 3,
    figure4 = figure(4);

    subplot(2,1,1);
    plot(spectra.trimmed_mean_diff(adjusted_xrange));
    hold

    title_text = 'Diff: Mean';
    Title(title_text);
    set(gca, 'FontWeight', 'Bold');
    xlim([0 600]);

    adjusted_mean_diff = spectra.trimmed_mean_diff(adjusted_xrange);
    spectra.adjusted_trimmed_mean_diff = adjusted_mean_diff;
 elseif num_spectra==3,
     adjusted_mean_diff = spectra.diff(adjusted_xrange);
     plot(adjusted_mean_diff);
     set(gca, 'FontWeight', 'Bold');
     xlim([0 600]);
     hold
 end

%         [pks, locs] = findpeaks(spectra.trimmed_mean_diff(adjusted_xrange), 'minpeakheight', dif_minimum_peak_height);                  
%         for i = 1:length(locs),
%             plotline(locs(i),'v');
%         end

    % Cr

%             GABA_peak_range = locs(4)-round(size(data,1)/120):locs(4)+round(size(data,1)/100);
%             GABA_baseline.left = locs(4)-round(size(data,1)/30):locs(4)-round(size(data,1)/60);            
%             GABA_baseline.right = locs(4)+round(size(data,1)/100):locs(4)+round(size(data,1)/60);

    GABA_peak_range = 273:348;
    GABA_baseline.left = 251:272;
    GABA_baseline.right = 386:436;

    plot(GABA_peak_range, adjusted_mean_diff(GABA_peak_range,:)','g', 'LineWidth',2);
    plot(GABA_baseline.left, adjusted_mean_diff(GABA_baseline.left,:)','r', 'LineWidth',2);
    plot(GABA_baseline.right, adjusted_mean_diff(GABA_baseline.right,:)','r', 'LineWidth',2);
    GABA.peak_area = sum(adjusted_mean_diff(GABA_peak_range));
    GABA.baseline_left_area = mean(adjusted_mean_diff(GABA_baseline.left)); % note mean instead of sum
    GABA.baseline_right_area = mean(adjusted_mean_diff(GABA_baseline.right)); % note mean instead of sum

    % Glx
    Glx_peak_range = 478:528;
    Glx_baseline.left = 473:477;
    Glx_baseline.right = 529:533;

    plot(Glx_peak_range, adjusted_mean_diff(Glx_peak_range,:)','g', 'LineWidth',2);
    plot(Glx_baseline.left, adjusted_mean_diff(Glx_baseline.left,:)','r', 'LineWidth',2);
    plot(Glx_baseline.right, adjusted_mean_diff(Glx_baseline.right,:)','r', 'LineWidth',2);
    Glx.peak_area = sum(adjusted_mean_diff(Glx_peak_range));
    Glx.baseline_left_area = sum(adjusted_mean_diff(Glx_baseline.left));
    Glx.baseline_right_area = sum(adjusted_mean_diff(Glx_baseline.right));

    legend('raw', 'peak', 'baseline', 'Location', 'SouthEast');
    legend boxoff;

%% Plot and Calculate Creatine Area Under Peak from Summed Spectra        
 if num_spectra > 3,
    subplot(2,1,2);

    adjusted_sum = spectra.trimmed_sum_On_w_Off(adjusted_xrange);
    spectra.adjusted_trimmed_sum_On_w_Off = adjusted_sum;
 elseif num_spectra == 3,
     subplot(4,1,4);
     adjusted_sum = on_and_off_spectra(adjusted_xrange);
 end
    plot(adjusted_sum);
    hold
    title_text = 'On & Off: Summed';
    Title(title_text);
    set(gca, 'FontWeight', 'Bold');
    xlim([0 600]);



    % Creatine
    Cr_peak_range = 295:337;
    Cr_baseline.left = 290:294;
    Cr_baseline.right = 338:342;

    plot(Cr_peak_range, adjusted_sum(Cr_peak_range,:)','g', 'LineWidth',2);
    plot(Cr_baseline.left, adjusted_sum(Cr_baseline.left,:)','r', 'LineWidth',2);
    plot(Cr_baseline.right, adjusted_sum(Cr_baseline.right,:)','r', 'LineWidth',2);
    Cr.peak_area = sum(adjusted_sum(Cr_peak_range));
    Cr.baseline_left_area = sum(adjusted_sum(Cr_baseline.left));
    Cr.baseline_right_area = sum(adjusted_sum(Cr_baseline.right));

%% Calculate NAA Area Under Peak from Off Spectra
 if num_spectra > 3,  
    adjusted_off = spectra.trimmed_mean_off(adjusted_xrange);
    spectra.adjusted_trimmed_mean_Off = adjusted_off;
 elseif num_spectra==3,
    adjusted_off = spectra.off(adjusted_xrange);
 end

    NAA_peak_range = 25:97;
    NAA_baseline.left = 15:24;
    NAA_baseline.right = 98:107;

    Cr.Off_peak_area = sum(adjusted_off(Cr_peak_range));
    Cr.Off_baseline_left_area = sum(adjusted_off(Cr_baseline.left));
    Cr.Off_baseline_right_area = sum(adjusted_off(Cr_baseline.right));

    NAA.peak_area = sum(adjusted_off(NAA_peak_range));
    NAA.baseline_left_area = sum(adjusted_off(NAA_baseline.left));
    NAA.baseline_right_area = sum(adjusted_off(NAA_baseline.right));

%% Calculate adjusted peak ratios
    % NOTE: These calculations were derived from the UC Davis Template
    % provided by Rick Maddock (rjmaddock@ucdavis.edu)
    GABA.adjusted_peak = GABA.peak_area-76*(((49/150)*(GABA.baseline_right_area-GABA.baseline_left_area))+GABA.baseline_left_area);
    Glx.adjusted_peak = Glx.peak_area-51*((Glx.baseline_left_area+Glx.baseline_right_area)/10);
    Cr.adjusted_peak = Cr.peak_area-43*((Cr.baseline_left_area+Cr.baseline_right_area)/10);
    Cr.Off_adjusted_peak = Cr.Off_peak_area-43*((Cr.Off_baseline_left_area+Cr.Off_baseline_right_area)/10);
    NAA.adjusted_peak = NAA.peak_area-73*((NAA.baseline_left_area+NAA.baseline_right_area)/20);

    Metabolites = struct('GABA', GABA, 'Glx', Glx, 'Cr', Cr, 'NAA', NAA);

    Peak_ratios.GABA_to_Cr = GABA.adjusted_peak/Cr.adjusted_peak;
    Peak_ratios.Glx_to_Cr = Glx.adjusted_peak/Cr.adjusted_peak;
    Peak_ratios.NAA_to_Cr = NAA.adjusted_peak/Cr.Off_adjusted_peak;

%%
    disp('Peak area and baseline regions have been calculated along with ratios. Hit any key to proceed.');
    pause
%% Save data
d = date;
outname = sprintf('S##_analyzed_%d_samples_%s', num_spectra, d);
uisave({'spectra' 'Peak_ratios', 'Metabolites', 'filename'}, outname);
        

end
