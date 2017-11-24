function AudVisUnlinkedKernel_v2
%% Import data 
%cd into file that holds data
cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/Data/');
%load in data
behavior = load('Combined_Data_BriannaUnlinked_171109_171110.mat');

%extract block led modes from matrix
behavior = behavior.all_data;
target = behavior(strcmp(behavior.LED, 'TARGET'), :);
masker = behavior(strcmp(behavior.LED, 'MASKER'), :);

%% From each category, extract 0% coherence trials only, separate by choice  
%hi linked target & masker with high choice
zm_HI = masker(masker.cohLevel == 0.25 & masker.choice == 2, :);
%lo linked target and masker with high choice
zm_LO = masker(masker.cohLevel == 0.25 & masker.choice == 1, :);

%target high choice 
zno_HI = target(target.cohLevel(:) == 0.25 & target.choice(:) == 2, :);
%target low choice 
zno_LO = target(target.cohLevel(:) == 0.25 & target.choice(:) == 1, :);

categories = {zm_HI, zm_LO, zno_HI, zno_LO};

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
        fdiff = fdiff(~isnan(fdiff));
        %convert to tone bursts
        fdiff = fdiff(1:50/1000 * fs:end);
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
%t = linspace(0, length(kernels(3,:))/fs * 1000, length(kernels(3,:)));
t = linspace(0, length(kernels(3,:)), length(kernels(3,:)));

subplot(2, 1, 1);
h1 = plot(t, kernels(1,:), 'Linewidth', 2); hold on;
h2 = plot(t, kernels(2,:), 'Linewidth', 2);
ylabel('Freq. (Hz)'); xlabel('# Tone Bursts'); %xlim([-timeb4 0]);
title('LED with Masker')
legend([h1, h2], 'High Choice', 'Low Choice', 'Location', 'NorthEastOutside')

subplot(2,1,2);
l1 = plot(t, kernels(3,:), 'Linewidth', 2); hold on;
l2 = plot(t, kernels(4,:), 'Linewidth', 2);
ylabel('Freq.(Hz)'); xlabel('# Tone Bursts'); %xlim([-timeb4 0]);
title('LED with Target')
legend([l1, l2], 'High Choice', 'Low Choice', 'Location', 'NorthEastOutside')

end

