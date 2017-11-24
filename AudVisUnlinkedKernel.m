function AudVisUnlinkedKernel
%% Import data 
%cd into file that holds data
cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/Data/');
%load in data
behavior = load('Combined_Data_Brianna_171027_171030.mat');

%extract block led modes from matrix
behavior = behavior.all_data;
target = behavior(strcmp(behavior.LED, 'TARGET'), :);
masker = behavior(strcmp(behavior.LED, 'MASKER'), :);

%extract waveforms for masker and target on trials when led is with masker
i = find(strcmp(behavior.LED, 'MASKER'));
twaves = cell2mat(table2array(behavior(i, 8)));
mwaves = cell2mat(table2array(behavior(i, 9)));
correlations = [];
for j = 1:size(twaves,1)
    t = twaves(j, :);
    m = mwaves(j, :);
    bothon = find(t ~= 0 & m ~= 0);
    corr = length(bothon)/length(t);
    correlations = [correlations, corr];
end 

% %quartile split
% lq = prctile(correlations, 25) 
% uq = prctile(correlations, 75) 
% masker_hi = masker(correlations >= uq, :); 
% masker_low = masker(correlations <= lq, :);

% %median split
med = median(correlations);
masker_hi = masker(correlations >= med, :);
masker_low = masker(correlations < med, :);

%% From each category, extract 0% coherence trials only, separate by choice  
%hi linked target & masker with high choice
zhi_HI = masker_hi(masker_hi.cohLevel == 0.5 & masker_hi.choice == 2, :);
%lo linked target and masker with high choice
zlo_HI = masker_low(masker_low.cohLevel == 0.5 & masker_low.choice == 2, :);

% %hi linked target & masker with high choice
% zhi_HI = masker_hi(masker_hi.cohLevel == 1 & masker_hi.choice == 2, :);
% %lo linked target and masker with high choice
% zlo_HI = masker_low(masker_low.cohLevel == 1 & masker_low.choice == 2, :);

%hi linked low choice
zhi_LO = masker_hi(masker_hi.cohLevel(:) == 0.5 & masker_hi.choice(:) == 1, :);
%low linked low choice 
zlo_LO = masker_low(masker_low.cohLevel(:) == 0.5 & masker_low.choice(:) == 1, :);

% %hi linked low choice
% zhi_LO = masker_hi(masker_hi.cohLevel(:) == 0 & masker_hi.choice(:) == 1, :);
% %low linked low choice 
% zlo_LO = masker_low(masker_low.cohLevel(:) == 0 & masker_low.choice(:) == 1, :);

%target high choice 
zno_HI = target(target.cohLevel(:) == 0.5 & target.choice(:) == 2, :);
%target low choice 
zno_LO = target(target.cohLevel(:) == 0.5 & target.choice(:) == 1, :);

% %target high choice 
% zno_HI = target(target.cohLevel(:) == 1 & target.choice(:) == 2, :);
% %target low choice 
% zno_LO = target(target.cohLevel(:) == 0 & target.choice(:) == 1, :);

categories = {zhi_HI, zlo_HI, zhi_LO, zlo_LO, zno_HI, zno_LO};

%% Find the mean frequencies of each trial, subtract the mean of each trial 
%from the frequencies in each trial (up to RT)  
fs = 24414;

kernels = [];
for i = 1:length(categories)
    %extract category, frequencies
    c = categories{i};
    f = c.frequencies;
    %create matrix to store all the vectors for each category
    catstore = [];
    for j = 1:size(categories{i}, 1)
        %extract frequency component for trial, find overall mean
        f1 = f{j};
        nsamples = floor(c.RT(j) / 1000 * fs);
        f1 = f1(1:nsamples);
        z = (f1 == 0);
        f1(z) = nan;
        av = nanmean(f1);
        sd = nanstd(f1, 1, 1);
        if sd == 0 
            sd = 1;
        end 
        %subtract mean from each time point 
        fdiff = (f1 - av) ./ sd;
        %remove time between tones 
        %fdiff = fdiff(~isnan(fdiff));
        %fdiff = fdiff(1:50/1000 * fs:end);
        %pad with nan's to max vector length
        fdiff = [fdiff' nan(1,97656-length(fdiff))];
        catstore(end+1, :) = fdiff;
    end 
    %take the mean of each time point for the category
    adj_av = nanmean(catstore,1);
    if size(categories{i}, 1) ~= 0
        kernels(end+1, :) = adj_av;
    else 
        kernels(end+1, :) = nan(1, 97656);
    end
end

%% Plot kernels 
t = linspace(0, length(kernels(3,:))/fs * 1000, length(kernels(3,:)));
%t = linspace(0, length(kernels(3,:)), length(kernels(3,:)));

subplot(3, 1, 1);
h1 = plot(t, kernels(1,:), 'Linewidth', 2); hold on;
h2 = plot(t, kernels(3,:), 'Linewidth', 2);
ylabel('Freq. (Hz)'); xlabel('Time (ms)'); %xlim([-timeb4 0]);
title('High Masker Coherence w/ Target')
legend([h1, h2], 'High Choice', 'Low Choice', 'Location', 'NorthEastOutside')

subplot(3,1,2);
l1 = plot(t, kernels(2,:), 'Linewidth', 2); hold on;
l2 = plot(t, kernels(4,:), 'Linewidth', 2);
ylabel('Freq.(Hz)'); xlabel('Time (ms)'); %xlim([-timeb4 0]);
title('Low Masker Coherence w/ Target')
legend([l1, l2], 'High Choice', 'Low Choice', 'Location', 'NorthEastOutside')

subplot(3,1,3);
t1 = plot(t, kernels(5,:), 'Linewidth', 2); hold on;
t2 = plot(t, kernels(6,:), 'Linewidth', 2);
ylabel('Freq. (Hz)'); xlabel('Time (ms)'); ylim([-2 2]);
title('LED With Target')
legend([t1, t2], 'High Choice', 'Low Choice', 'Location', 'NorthEastOutside')

% %% Determine average response time for each coherence (intended) 
% %find unique coherence levels 
% cohs = unique(behavior.cohLevel);
% cohRT = [];
% ciRT1 = [];
% for i = 1:length(cohs) 
%     %get trials with those coherence levels 
%     rel_RT = behavior.RT(behavior.cohLevel == cohs(i));
%     cohRT(end+1) = mean(rel_RT);
%     sd = std(rel_RT);
%     conf = 1.96 * sd / sqrt(length(rel_RT));
%     ciRT1(end+1) = conf;
% end
% 
% figure();
% err = errorbar(cohs.*200-100, cohRT, ciRT1, '-ok', 'Linewidth', 1);
% xlim([-100 100]); ylim([0 2000])
% xlabel('Signed Coherence'); ylabel('Reaction Time (ms)'); hold on;
   
% %% Determine average response time for each coherence (played) 
% pcohs = unique(behavior.coh_played);
% pcohRT = [];
% ciRT2 = [];
% realpcohs = [];
% 
% bins = 13; 
% n = floor(length(pcohs)/bins);
% for i = 1:n:length(pcohs)
%     x = i + n-1;
%     if i + n-1 > length(pcohs)
%         x = length(pcohs);
%     end
%     %get trials with those coherence levels 
%     rel_RT = behavior.RT(behavior.coh_played >= pcohs(i) & behavior.coh_played <= pcohs(x));
%     pcohRT(end+1) = mean(rel_RT);
%     sd = std(rel_RT);
%     conf = 1.96 * sd / sqrt(length(rel_RT));
%     ciRT2(end+1) = conf;
%     realpcohs(end+1) = mean(pcohs(i:x));
% end
% 
% p = plot(realpcohs.*200-100, pcohRT, 'xr', 'MarkerSize', 6);
% legend([err p], 'Specified Coherence', 'Played Coherence')
end

