%% December 2019: Programma om mijn power spectra te verwerken met FOOOF (Haller et al. 2018)

% The program takes as an input .dat files containing the power spectrum of
% LFP recordings calculated with ZebraExplore

% The output of the program gives measures of the goodness of fit of the
% exponential function, the slope and offset of the exponential, and
% detected 'oscillations' on top of the 1/f exponential function. Finally,
% it calculates entropy as in Foffani et al. 2007. 

% The output is copied to the clipboard for further processing in excel. 

%%  Input information

% Important: you should have created a row vector with all the names of the
% paths where the PWS files are located for this animal

foldertosave = 'F:\Data to be analyzed\LFP P25 with GABA\191015 number 2\191015\fooof1';

input_paths = ['F:\Data to be analyzed\LFP P25 with GABA\191015 number 2\191015\Dir_10\PWS.dat'
 ];

% settings of fooof object
settings = struct();
settings.peak_threshold = 3;
settings.max_n_peaks = Inf;
settings.min_peak_amplitude = 0;
settings.peak_width_limits = [0.5, 200];

% Import the data
[numberdirs, temp] = size(input_paths);
out = struct(); 
for dir = 1:numberdirs
    data = importdata(input_paths(dir, :));
    [lengthfreqs, columnsdata] = size(data);
    freqs = data(:,1)';  
    nCh = 2;
    pws = zeros(lengthfreqs, nCh);
    
    % Run the fooof functions on the two frequency bands

    for i = 1:nCh
        pws(:,i) = data(:, i+1)';  
        f_range = [0.5, 25];
        fooof_results.deltabeta(i) = fooof(freqs, pws(:,i), f_range, settings);
        f_range = [40, 100];
        fooof_results.gamma(i) = fooof(freqs, pws(:,i), f_range, settings);
        f_range = [30, 50];
        fooof_results.ei(i) = fooof(freqs, pws(:,i), f_range, settings);
    end

    %% plotting to check  % this is useful if you suspect fitting is not correct. I used it to optimise FOOOF settings.
    %
    f_range = [0.5, 25];
    exp_ch1 = 10.^fooof_results.deltabeta(1).background_params(1)*(1./(0+freqs.^fooof_results.deltabeta(1).background_params(2)));
    exp_ch2 = 10.^fooof_results(2).background_params(1)*(1./(0+freqs.^fooof_results(2).background_params(2)));
    [value_start start] = min(abs(freqs-f_range(1)));
    [value_stop stop] = min(abs(freqs-f_range(2)));
     
    plotch1_log = figure('Name', 'Ch1-logscale');
    plot(log10(freqs(start:stop)), log10(exp_ch1(start:stop)), 'r', log10(freqs(start:stop)), log10(pws(start:stop, 1)), 'g');
    
    plotch1 = figure('Name', 'Ch1');
    plot(freqs(start:stop), exp_ch1(start:stop), 'r', freqs(start:stop), pws(start:stop, 1), 'g');
     
    plotch2_log = figure('Name', 'Ch2-logscale');
    plot(log10(freqs(start:stop)), log10(exp_ch2(start:stop)), 'r', log10(freqs(start:stop)), log10(pws(start:stop, 2)), 'g');
     
    plotch2 = figure('Name', 'Ch2');
    plot(freqs(start:stop), exp_ch2(start:stop), 'r', freqs(start:stop), pws(start:stop, 2), 'g');


    %% Process output: 1) save as a struct
    
    out(dir).pws.ch1 = pws(:,1);
    out(dir).pws.ch2 = pws(:,2);
    out(dir).freqs = freqs';
    out(dir).db.ch1 = fooof_results.deltabeta(1);
    out(dir).db.ch2 = fooof_results.deltabeta(2);
    out(dir).g.ch1 = fooof_results.gamma(1);
    out(dir).g.ch2 = fooof_results.gamma(2);
    out(dir).ei.ch1 = fooof_results.ei(1);
    out(dir).ei.ch2 = fooof_results.ei(2);
    
end

%% Process output: 2) plot the histogram of the frequencies

% First create a row vector with all peak frequencies

peaks_freq_ch1 = [out(1).db.ch1.peak_params(:,1); out(1).g.ch1.peak_params(:,1)];
for c = 2:numberdirs
    pfc1 = [out(c).db.ch1.peak_params(:,1); out(c).g.ch1.peak_params(:,1)];
    peaks_freq_ch1 = [peaks_freq_ch1; pfc1];
end
peaks_freq_ch2 = [out(1).db.ch2.peak_params(:,1); out(1).g.ch2.peak_params(:,1)];
for c = 2:numberdirs
    pfc2 = [out(c).db.ch2.peak_params(:,1); out(c).g.ch2.peak_params(:,1)];
    peaks_freq_ch2 = [peaks_freq_ch2; pfc2];
end

% then plot it all in a histogram

hist = figure('Name', 'Hist_Freqs');
histch1 = histogram(peaks_freq_ch1, 50, 'FaceColor', 'r');
hold on 
histch2 = histogram(peaks_freq_ch2, 50, 'FaceColor', 'k');

%% Create output: the oscillations

% I will process those based on the average of the entire animal: look at
% the resulting histogram. From this I will find the peaks: their central
% frequencies, their width, and their count. 

% Then as a better output for amplitude I will include the amplitude values
% found by fooof within the peak boundaries and will average this

% So to start: let's find the peaks in the histogram

% for ch1:

% First let's find the frequency values at the center of each bin
bincentersch1 = histch1.BinEdges(1:end-1)+(histch1.BinWidth/2);

% Create a new histogram count file with added values before and after the
% histogram. This is needed since the findpeaks function otherwise ignores 
% peaks located on the edges. I did this by simply mirroring the first 3
% values of the histogram and pasting it at the start of the new histogram
% file, while mirroring the end of the histogram and pasting it at the end

newhistch1 = [histch1.Values(4:-1:2), histch1.Values, histch1.Values(end-1:-1:end-3)];

% Then we also need to extend the bincenters array
partone_bincentersch1 = bincentersch1-(3*histch1.BinWidth):histch1.BinWidth:bincentersch1-histch1.BinWidth;
partthree_bincentersch1 = bincentersch1(end)+histch1.BinWidth:histch1.BinWidth:bincentersch1(end)+(3*histch1.BinWidth);
newbincentersch1 = [partone_bincentersch1, bincentersch1, partthree_bincentersch1];

% then we can find the peaks and their amplitudes, central frequency, and
% widths
[peakampch1, peakfreqch1, peakwidthch1] = findpeaks(newhistch1, newbincentersch1, 'MinPeakHeight', 5); 
    % add a minimal height of the peaks to really include only the important ones
    % Note that the width in this case is the width at half the prominence
    % of the peak (so how much the peak exceeds from the next trough)
    
% Finally, if we accidentally created a peak outside of the real range it has to be removed:
nonvalidsch1 = peakfreqch1 < bincentersch1(1)|peakfreqch1 > bincentersch1(end);
peakampch1(nonvalidsch1) = [];
peakfreqch1(nonvalidsch1) = [];
peakwidthch1(nonvalidsch1) = [];

% Save the results
results(1).oscillation_amp_histcount = peakampch1;
results(1).oscillation_freq = peakfreqch1;
results(1).oscillation_width = peakwidthch1;

% Now we do the same for channel 2
bincentersch2 = histch2.BinEdges(1:end-1)+(histch2.BinWidth/2);
newhistch2 = [histch2.Values(4:-1:2), histch2.Values, histch2.Values(end-1:-1:end-3)];
partone_bincentersch2 = bincentersch2-(3*histch2.BinWidth):histch2.BinWidth:bincentersch2-histch2.BinWidth;
partthree_bincentersch2 = bincentersch2(end)+histch2.BinWidth:histch2.BinWidth:bincentersch2(end)+(3*histch2.BinWidth);
newbincentersch2 = [partone_bincentersch2, bincentersch2, partthree_bincentersch2];
[peakampch2, peakfreqch2, peakwidthch2] = findpeaks(newhistch2, newbincentersch2, 'MinPeakHeight', 5); 
nonvalidsch2 = peakfreqch2 < bincentersch2(1)|peakfreqch2 > bincentersch2(end);
peakampch2(nonvalidsch2) = [];
peakfreqch2(nonvalidsch2) = [];
peakwidthch2(nonvalidsch2) = [];

% Save the results
results(2).oscillation_amp_histcount = peakampch2;
results(2).oscillation_freq = peakfreqch2;
results(2).oscillation_width = peakwidthch2;

% The final processing involves including the average 'real amplitudes', frequency, and SD (two times standard
% deviation of the estimated Gaussian of fooof) in the three band passes of my interest

% For Ch1

% First initiate all the values
delta_amp_ch1 = 0;
beta_amp_ch1 = 0;
gamma_amp_ch1 = 0;

delta_freq_ch1 = 0;
beta_freq_ch1 = 0;
gamma_freq_ch1 = 0;

delta_sd_ch1 = 0;
beta_sd_ch1 = 0;
gamma_sd_ch1 = 0;

count_delta_ch1 = 0;
count_beta_ch1 = 0;
count_gamma_ch1 = 0;


% Loop through all directories, and then all the peaks in each directory,
% sum up all amplitudes, SDs, and frequencies and count them
for dir = 1:numberdirs
    tempdb = size(out(dir).db.ch1.peak_params);
    peaknumberdb = tempdb(1);
    for nr = 1:peaknumberdb 
        if out(dir).db.ch1.peak_params(nr, 1) <= 4
            delta_amp_ch1 = delta_amp_ch1 + out(dir).db.ch1.peak_params(nr,2);
            delta_sd_ch1 = delta_sd_ch1 + out(dir).db.ch1.peak_params(nr,3);
            delta_freq_ch1 = delta_freq_ch1 + out(dir).db.ch1.peak_params(nr,1);
            count_delta_ch1 = count_delta_ch1 + 1;
        else
            beta_amp_ch1 = beta_amp_ch1 + out(dir).db.ch1.peak_params(nr,2);
            beta_sd_ch1 = beta_sd_ch1 + out(dir).db.ch1.peak_params(nr,3);
            beta_freq_ch1 = beta_freq_ch1 + out(dir).db.ch1.peak_params(nr,1);
            count_beta_ch1 = count_beta_ch1 + 1;
        end    
    end
    tempg = size(out(dir).g.ch1.peak_params);
    peaknumberg = tempg(1);
    for nr = 1:peaknumberg
        gamma_amp_ch1 = gamma_amp_ch1 + out(dir).g.ch1.peak_params(nr,2);
        gamma_sd_ch1 = gamma_sd_ch1 + out(dir).g.ch1.peak_params(nr,3);
        gamma_freq_ch1 = gamma_freq_ch1 + out(dir).g.ch1.peak_params(nr, 1);
        count_gamma_ch1 = count_gamma_ch1 + 1;
    end
end

% Then average it all and save in results struct
results(1).average_amp_delta = delta_amp_ch1/count_delta_ch1;
results(1).average_amp_beta = beta_amp_ch1/count_beta_ch1;
results(1).average_amp_gamma = gamma_amp_ch1/count_gamma_ch1;

results(1).average_sd_delta = delta_sd_ch1/count_delta_ch1;
results(1).average_sd_beta = beta_sd_ch1/count_beta_ch1;
results(1).average_sd_gamma = gamma_sd_ch1/count_gamma_ch1;

results(1).average_freq_delta = delta_freq_ch1/count_delta_ch1;
results(1).average_freq_beta = beta_freq_ch1/count_beta_ch1;
results(1).average_freq_gamma = gamma_freq_ch1/count_gamma_ch1;

% Now we do the same for Ch2
delta_amp_ch2 = 0;
beta_amp_ch2 = 0;
gamma_amp_ch2 = 0;

delta_freq_ch2 = 0;
beta_freq_ch2 = 0;
gamma_freq_ch2 = 0;

delta_sd_ch2 = 0;
beta_sd_ch2 = 0;
gamma_sd_ch2 = 0;

count_delta_ch2 = 0;
count_beta_ch2 = 0;
count_gamma_ch2 = 0;


% Loop through all directories, and then all the peaks in each directory,
% sum up all amplitudes, SDs, and frequencies and count them
for dir = 1:numberdirs
    tempdb = size(out(dir).db.ch2.peak_params);
    peaknumberdb = tempdb(1);
    for nr = 1:peaknumberdb 
        if out(dir).db.ch2.peak_params(nr, 1) <= 4
            delta_amp_ch2 = delta_amp_ch2 + out(dir).db.ch2.peak_params(nr,2);
            delta_sd_ch2 = delta_sd_ch2 + out(dir).db.ch2.peak_params(nr,3);
            delta_freq_ch2 = delta_freq_ch2 + out(dir).db.ch2.peak_params(nr,1);
            count_delta_ch2 = count_delta_ch2 + 1;
        else
            beta_amp_ch2 = beta_amp_ch2 + out(dir).db.ch2.peak_params(nr,2);
            beta_sd_ch2 = beta_sd_ch2 + out(dir).db.ch2.peak_params(nr,3);
            beta_freq_ch2 = beta_freq_ch2 + out(dir).db.ch2.peak_params(nr,1);
            count_beta_ch2 = count_beta_ch2 + 1;
        end    
    end
    tempg = size(out(dir).g.ch2.peak_params);
    peaknumberg = tempg(1);
    for nr = 1:peaknumberg
        gamma_amp_ch2 = gamma_amp_ch2 + out(dir).g.ch2.peak_params(nr,2);
        gamma_sd_ch2 = gamma_sd_ch2 + out(dir).g.ch2.peak_params(nr,3);
        gamma_freq_ch2 = gamma_freq_ch2 + out(dir).g.ch2.peak_params(nr, 1);
        count_gamma_ch2 = count_gamma_ch2 + 1;
    end
end

% Then average it all and save in results struct
results(2).average_amp_delta = delta_amp_ch2/count_delta_ch2;
results(2).average_amp_beta = beta_amp_ch2/count_beta_ch2;
results(2).average_amp_gamma = gamma_amp_ch2/count_gamma_ch2;

results(2).average_sd_delta = delta_sd_ch2/count_delta_ch2;
results(2).average_sd_beta = beta_sd_ch2/count_beta_ch2;
results(2).average_sd_gamma = gamma_sd_ch2/count_gamma_ch2;

results(2).average_freq_delta = delta_freq_ch2/count_delta_ch2;
results(2).average_freq_beta = beta_freq_ch2/count_beta_ch2;
results(2).average_freq_gamma = gamma_freq_ch2/count_gamma_ch2;

savefig('hist_freqs.fig');
movefile('hist_freqs.fig', foldertosave);
disp('Figure saved: do you want to continue?')
pause;
clc; close;

%% Process output: 3) plot the error and r squared

% initiation of the vector and finding the values needed

error_db_ch1 = zeros(numberdirs, 1);
error_db_ch2 = zeros(numberdirs, 1);
error_g_ch1 = zeros(numberdirs, 1);
error_g_ch2 = zeros(numberdirs, 1);
error_ei_ch1 = zeros(numberdirs, 1);
error_ei_ch2 = zeros(numberdirs, 1);

rsq_db_ch1 = zeros(numberdirs, 1);
rsq_db_ch2 = zeros(numberdirs, 1);
rsq_g_ch1 = zeros(numberdirs, 1);
rsq_g_ch2 = zeros(numberdirs, 1);
rsq_ei_ch1 = zeros(numberdirs, 1);
rsq_ei_ch2 = zeros(numberdirs, 1);

for c = 1:numberdirs
    error_db_ch1(c,1) = out(c).db.ch1.error;
    error_db_ch2(c,1) = out(c).db.ch2.error;
    error_g_ch1(c,1) = out(c).g.ch1.error;
    error_g_ch2(c,1) = out(c).g.ch2.error;
    error_ei_ch1(c,1) = out(c).ei.ch1.error;
    error_ei_ch2(c,1) = out(c).ei.ch2.error;
    
    rsq_db_ch1(c,1) = out(c).db.ch1.r_squared;
    rsq_db_ch2(c,1) = out(c).db.ch2.r_squared;
    rsq_g_ch1(c,1) = out(c).g.ch1.r_squared;
    rsq_g_ch2(c,1) = out(c).g.ch2.r_squared;
    rsq_ei_ch1(c,1) = out(c).ei.ch1.r_squared;
    rsq_ei_ch2(c,1) = out(c).ei.ch2.r_squared;
end

% Now also find averages and save in out struct

results(1).average_error_db = mean(error_db_ch1);
results(2).average_error_db = mean(error_db_ch2);
results(1).average_error_g = mean(error_g_ch1);
results(2).average_error_g = mean(error_g_ch2);
results(1).average_error_ei = mean(error_ei_ch1);
results(2).average_error_ei = mean(error_ei_ch2);

results(1).average_rsq_db = mean(rsq_db_ch1);
results(2).average_rsq_db = mean(rsq_db_ch2);
results(1).average_rsq_g = mean(rsq_g_ch1);
results(2).average_rsq_g = mean(rsq_g_ch2);
results(1).average_rsq_ei = mean(rsq_ei_ch1);
results(2).average_rsq_ei = mean(rsq_ei_ch2);

% Then plot

x1 = 0.5.*ones(1,numberdirs);
x2 = ones(1,numberdirs);
x3 = 1.5.*ones(1,numberdirs);
x4 = 2.*ones(1,numberdirs);

scattererror = figure('Name', 'Error');
plot(x1, error_db_ch1,'r*');
hold on;
plot(x2, error_db_ch2, 'k*');
hold on;
plot(x3, error_g_ch1, 'r*');
hold on;
plot(x4, error_g_ch2, 'k*');
xlim([0, 2.5]);

savefig('error_scatter.fig');
movefile('error_scatter.fig', foldertosave);

scatterrsq = figure('Name', 'R Squared');
plot(x1, rsq_db_ch1, 'r*');
hold on;   
plot(x2, rsq_db_ch2, 'k*');
hold on;
plot(x3, rsq_g_ch1, 'r*');
hold on;
plot(x4, rsq_g_ch2, 'k*');
xlim([0,2.5]);

savefig('r_squared_scatter.fig');
movefile('r_squared_scatter.fig', foldertosave);


disp('Figures saved: do you want to continue?')
pause;
clc; close all;


%% Process output: 4) plot the histogram of the slope values

% initiation of the vector and finding the values needed

offset_db_ch1 = zeros(numberdirs, 1);
offset_db_ch2 = zeros(numberdirs, 1);
offset_g_ch1 = zeros(numberdirs, 1);
offset_g_ch2 = zeros(numberdirs, 1);
offset_ei_ch1 = zeros(numberdirs, 1);
offset_ei_ch2 = zeros(numberdirs, 1);

exp_db_ch1 = zeros(numberdirs, 1);
exp_db_ch2 = zeros(numberdirs, 1);
exp_g_ch1 = zeros(numberdirs, 1);
exp_g_ch2 = zeros(numberdirs, 1);
exp_ei_ch1 = zeros(numberdirs, 1);
exp_ei_ch2 = zeros(numberdirs, 1);

for c = 1:numberdirs
    offset_db_ch1(c,1) = out(c).db.ch1.background_params(1);
    offset_db_ch2(c,1) = out(c).db.ch2.background_params(1);
    offset_g_ch1(c,1) = out(c).g.ch1.background_params(1);
    offset_g_ch2(c,1) = out(c).g.ch2.background_params(1);
    offset_ei_ch1(c,1) = out(c).ei.ch1.background_params(1);
    offset_ei_ch2(c,1) = out(c).ei.ch2.background_params(1);
    
    exp_db_ch1(c,1) = out(c).db.ch1.background_params(2);
    exp_db_ch2(c,1) = out(c).db.ch2.background_params(2);
    exp_g_ch1(c,1) = out(c).g.ch1.background_params(2);
    exp_g_ch2(c,1) = out(c).g.ch2.background_params(2);
    exp_ei_ch1(c,1) = out(c).ei.ch1.background_params(2);
    exp_ei_ch2(c,1) = out(c).ei.ch2.background_params(2);
end

% Now also find averages and save in out struct

results(1).average_offset_db = mean(offset_db_ch1);
results(2).average_offset_db = mean(offset_db_ch2);
results(1).average_offset_g = mean(offset_g_ch1);
results(2).average_offset_g = mean(offset_g_ch2);
results(1).average_offset_ei = mean(offset_ei_ch1);
results(2).average_offset_ei = mean(offset_ei_ch2);

results(1).average_exp_db = mean(exp_db_ch1);
results(2).average_exp_db = mean(exp_db_ch2);
results(1).average_exp_g = mean(exp_g_ch1);
results(2).average_exp_g = mean(exp_g_ch2);
results(1).average_exp_ei = mean(exp_ei_ch1);
results(2).average_exp_ei = mean(exp_ei_ch2);

% Then plot

scatteroffset = figure('Name', 'Offset');
plot(x1, offset_db_ch1,'r*');
hold on;
plot(x2, offset_db_ch2, 'k*');
hold on;
plot(x3, offset_g_ch1, 'r*');
hold on;
plot(x4, offset_g_ch2, 'k*');
xlim([0, 2.5]);

savefig('offset_scatter.fig');
movefile('offset_scatter.fig', foldertosave);

scatterexp = figure('Name', 'Exponent');
plot(x1, exp_db_ch1, 'r*');
hold on;   
plot(x2, exp_db_ch2, 'k*');
hold on;
plot(x3, exp_g_ch1, 'r*');
hold on;
plot(x4, exp_g_ch2, 'k*');
xlim([0,2.5]);

savefig('exp_scatter.fig');
movefile('exp_scatter.fig', foldertosave);


disp('Figures saved: do you want to continue?')
pause;
clc; close all;

%% Calculate the entropy

% First normalize the pws

% we do this for each directory, then we average

% initiate the variables that are going to hold the entropy for each
% directory
entropych1 = zeros(numberdirs, 1);
entropydeltach1 = zeros(numberdirs, 1);
entropybetach1 = zeros(numberdirs, 1);
entropygammach1 = zeros(numberdirs, 1);

entropych2 = zeros(numberdirs, 1);
entropydeltach2 = zeros(numberdirs, 1);
entropybetach2 = zeros(numberdirs, 1);
entropygammach2 = zeros(numberdirs, 1);


for dir = 1:numberdirs    
    % first find the 0.5-100 Hz range of the pws and freqs variables
    [temp, loc05] = min(abs(out(dir).freqs-0.5));
    [temp, loc100] = min(abs(out(dir).freqs-100));
    [temp, loc4] = min(abs(out(dir).freqs-4));
    [temp, loc9] = min(abs(out(dir).freqs-9));
    [temp, loc25] = min(abs(out(dir).freqs-25));
    [temp, loc40] = min(abs(out(dir).freqs-40));
    for bp = 1:4
        if bp == 1
            st = loc05;
            stop = loc100;
        elseif bp == 2
            st = loc05;
            stop = loc4;
        elseif bp == 3
            st = loc9;
            stop = loc25;
        elseif bp == 4
            st = loc40;
            stop = loc100;
        end
        pws_ch1 = out(dir).pws.ch1(st:stop);
        pws_ch2 = out(dir).pws.ch2(st:stop);
        freqsnew = out(dir).freqs(st:stop);
    
        % Then normalize: formula = (x-xmin)/(xmax-xmin)
        
        % Turns out this formula is not used for power spectrum: I need to
        % normalize in such a way as to visualize the signal as a
        % probability distribution. So code below is commented out      
    
%         xmin_ch1 = min(pws_ch1);
%         xmin_ch2 = min(pws_ch2);
%         xmax_ch1 = max(pws_ch1);
%         xmax_ch2 = max(pws_ch2);
%     
%         pwsnorm_ch1 = zeros(length(freqsnew),1);
%         pwsnorm_ch2 = zeros(length(freqsnew),1);
%     
%         for x = 1:length(freqsnew)
%             pwsnorm_ch1(x) = (pws_ch1(x)-xmin_ch1)/(xmax_ch1-xmin_ch1);
%             pwsnorm_ch2(x) = (pws_ch2(x)-xmin_ch2)/(xmax_ch2-xmin_ch2);
%         end

        % Follow new formula according to: https://www.mathworks.com/help/signal/ref/pentropy.html
        % and https://dsp.stackexchange.com/questions/23689/what-is-spectral-entropy
        
        % powerdensity_ch1 = abs(pws_ch1.^2)/length(freqsnew)
        powerdensity_ch1 = abs(pws_ch1).^2;
        pwsnorm_ch1 = powerdensity_ch1./sum(powerdensity_ch1);
        
        powerdensity_ch2 = abs(pws_ch2).^2;        
        pwsnorm_ch2 = powerdensity_ch2./sum(powerdensity_ch2);
    
        % Then continue to calclate the entropy according to formula             
    
        for f = 1:length(freqsnew)
            if bp == 1
                if pwsnorm_ch1(f)==0
                    ent_ch1(f) = NaN;
                else
                    ent_ch1(f) = pwsnorm_ch1(f)*log2(pwsnorm_ch1(f));
                end
                if pwsnorm_ch2(f)==0
                    ent_ch2(f) = NaN;
                else
                    ent_ch2(f) = pwsnorm_ch2(f)*log2(pwsnorm_ch2(f));
                end
            elseif bp == 2
                if pwsnorm_ch1(f)==0
                    ent_ch1_d(f) = NaN;
                else
                    ent_ch1_d(f) = pwsnorm_ch1(f)*log2(pwsnorm_ch1(f));
                end
                if pwsnorm_ch2(f)==0
                    ent_ch2_d(f) = NaN;
                else
                    ent_ch2_d(f) = pwsnorm_ch2(f)*log2(pwsnorm_ch2(f));
                end
            elseif bp == 3
                if pwsnorm_ch1(f)==0
                    ent_ch1_b(f) = NaN;
                else
                    ent_ch1_b(f) = pwsnorm_ch1(f)*log2(pwsnorm_ch1(f));
                end
                if pwsnorm_ch2(f)==0
                    ent_ch2_b(f) = NaN;
                else
                    ent_ch2_b(f) = pwsnorm_ch2(f)*log2(pwsnorm_ch2(f));
                end
            elseif bp == 4
                if pwsnorm_ch1(f)==0
                    ent_ch1_g(f) = NaN;
                else
                    ent_ch1_g(f) = pwsnorm_ch1(f)*log2(pwsnorm_ch1(f));
                end
                if pwsnorm_ch2(f)==0
                    ent_ch2_g(f) = NaN;
                else
                    ent_ch2_g(f) = pwsnorm_ch2(f)*log2(pwsnorm_ch2(f));
                end
            end            
        end        
    end
    entropych1(dir) = -nansum(ent_ch1);
    entropydeltach1(dir) = -nansum(ent_ch1_d);
    entropybetach1(dir) = -nansum(ent_ch1_b);
    entropygammach1(dir) = -nansum(ent_ch1_g);
    
    entropych2(dir) = -nansum(ent_ch2);
    entropydeltach2(dir) = -nansum(ent_ch2_d);
    entropybetach2(dir) = -nansum(ent_ch2_b);
    entropygammach2(dir) = -nansum(ent_ch2_g);
end

results(1).entropy = mean(entropych1);
results(1).entropy_d = mean(entropydeltach1);
results(1).entropy_b = mean(entropybetach1);
results(1).entropy_g = mean(entropygammach1);

results(2).entropy = mean(entropych2);
results(2).entropy_d = mean(entropydeltach2);
results(2).entropy_b = mean(entropybetach2);
results(2).entropy_g = mean(entropygammach2);

%% Finally print the output to save

% Row 1: frequency peak
str = ' ';

for ch = 1:2
    str = sprintf('%s\t', str, 'Oscillation Freq');
    for le = 1:length(results(ch).oscillation_freq)
        if le == length(results(ch).oscillation_freq)
            str = sprintf('%s%f', str, results(ch).oscillation_freq(le));
        else 
            str = sprintf('%s%f\t', str, results(ch).oscillation_freq(le));
        end
    end
    str = sprintf('%s\n', str);

    % Row 2: amplitude peak
    str = sprintf('%s\t', str, 'Oscillation Amplitude');
    for le = 1:length(results(ch).oscillation_freq)
            if le == length(results(ch).oscillation_freq)
                    str = sprintf('%s%f', str, results(ch).oscillation_amp_histcount(le));
            else
                    str = sprintf('%s%f\t', str, results(ch).oscillation_amp_histcount(le));
            end
    end
    str = sprintf('%s\n', str);

    % Row 3: width peak 
    str = sprintf('%s\t', str, 'Oscillation Width');
    for le = 1:length(results(ch).oscillation_freq)
            if le == length(results(ch).oscillation_freq)
                    str = sprintf('%s%f', str, results(ch).oscillation_width(le));
            else
                    str = sprintf('%s%f\t', str, results(ch).oscillation_width(le));
            end
    end
    str = sprintf('%s\n', str);

    % Row 4: delta oscillation freq, amp, sd
    str = sprintf('%s\t', str, 'Delta Oscillations');
    str = sprintf('%s%f\t', str, results(ch).average_freq_delta);
    str = sprintf('%s%f\t', str, results(ch).average_amp_delta);
    str = sprintf('%s%f', str, results(ch).average_sd_delta);
    str = sprintf('%s\n', str);

    % Row 5: beta oscillation freq, amp, sd
    str = sprintf('%s\t', str, 'Beta Oscillations');
    str = sprintf('%s%f\t', str, results(ch).average_freq_beta);
    str = sprintf('%s%f\t', str, results(ch).average_amp_beta);
    str = sprintf('%s%f', str, results(ch).average_sd_beta);
    str = sprintf('%s\n', str);

    % Row 6: gamma oscillation freq, amp, sd
    str = sprintf('%s\t', str, 'Gamma Oscillations');
    str = sprintf('%s%f\t', str, results(ch).average_freq_gamma);
    str = sprintf('%s%f\t', str, results(ch).average_amp_gamma);
    str = sprintf('%s%f', str, results(ch).average_sd_gamma);
    str = sprintf('%s\n', str);

    % Row 7: Error and r squared DeltaBeta
    str = sprintf('%s\t', str, 'Error and RSQ DeltaBeta');
    str = sprintf('%s%f\t', str, results(ch).average_error_db);
    str = sprintf('%s%f', str, results(ch).average_rsq_db);
    str = sprintf('%s\n', str);

    % Row 8: Error and r squared Gamma
    str = sprintf('%s\t', str, 'Error and RSQ Gamma');
    str = sprintf('%s%f\t', str, results(ch).average_error_g);
    str = sprintf('%s%f', str, results(ch).average_rsq_g);
    str = sprintf('%s\n', str);
    
    % Row 9: Error and r squared E/I ratio (30-50 Hz)
    str = sprintf('%s\t', str, 'Error and RSQ E/I');
    str = sprintf('%s%f\t', str, results(ch).average_error_ei);
    str = sprintf('%s%f', str, results(ch).average_rsq_ei);
    str = sprintf('%s\n', str);

    % Row 10: Offset and Exp DeltaBeta
    str = sprintf('%s\t', str, 'Offset and Exp DeltaBeta');
    str = sprintf('%s%f\t', str, results(ch).average_offset_db);
    str = sprintf('%s%f', str, results(ch).average_exp_db);
    str = sprintf('%s\n', str);

    % Row 11: Offset and Exp Gamma
    str = sprintf('%s\t', str, 'Offset and Exp Gamma');
    str = sprintf('%s%f\t', str, results(ch).average_offset_g);
    str = sprintf('%s%f', str, results(ch).average_exp_g);
    str = sprintf('%s\n', str);
    
    % Row 12: Offset and Exp E/I ratio (30-50 Hz)
    str = sprintf('%s\t', str, 'Offset and Exp E/I');
    str = sprintf('%s%f\t', str, results(ch).average_offset_ei);
    str = sprintf('%s%f', str, results(ch).average_exp_ei);
    str = sprintf('%s\n', str);
    
    % Row 13: Entropy
    str = sprintf('%s\t', str, 'Entropy');
    str = sprintf('%s%f\t', str, results(ch).entropy);
    str = sprintf('%s%f\t', str, results(ch).entropy_d);
    str = sprintf('%s%f\t', str, results(ch).entropy_b);
    str = sprintf('%s%f', str, results(ch).entropy_g);
    str = sprintf('%s\n', str);
    
    str = sprintf('%s\n', str);
    str = sprintf('%s\n', str);       
end

clipboard ('copy',str);
'finished!'