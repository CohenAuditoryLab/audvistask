function AudVis_LogisticRegression
%% Fits a logistic function to task data 
    %   11/13/17 BMK
    %   Simultaneous Fits!! 
    %   Considering LED as one parameter.
%% File Selection 
% Select 'all' or 'combined' file
disp('Select task data (mac)');
[task_file, tpath] = uigetfile('*.mat');
% Load data 
load([tpath task_file]);
disp([tpath task_file]);
data = all_data;
%% Data manipulation 

% extract from table: 
%   Column 1: Played Coherence
%   Column 2: LED
%   Column 3: Choice
%   Column 4: RT 
coh = data.coh_played.*200-100; %converted to signed 
led = data.LED;
ch = data.choice - 1; %1 = high, 0 = low

%code LED data: 1 = target, -1 = masker
m = cellfun(@(x) strcmp(x,'MASKER'), led);
t = cellfun(@(y) strcmp(y, 'TARGET'), led);
led = zeros(length(m), 1);
led(m) = -1;
led(t) = 1;

%treat data as one matrix
%independent variables, then dependent variable (choice)
x = [led coh ch];

%get all target and all masker trials in separate data structures 
masker = data(m,:);
target = data(t,:);

mcoh = masker.coh_played.*200 - 100;
tcoh = target.coh_played.*200 - 100;
mled = -1 .* ones(length(mcoh), 1);
tled = ones(length(tcoh), 1);
mch = masker.choice - 1;
tch = target.choice - 1;

mx = [mled mcoh mch];
tx = [tled tcoh tch];

%bin coherences
cbins = [ ...
    -101  -99
    -99   -60
    -60   -40
    -40   -28
    -28   -12
    -12    -8
     -8     8
      8    12
     12    28
     28    40
     40    60
     60    99
     99   101];
%calculate number of coherence bins
nbins = size(cbins,1);
%calculate number of trials
ntrials = size(x, 1);
%calculate average value of each coherence bin - output vector
cax = mean(cbins,2);
pmf1   = NaN(nbins,1);
cmf1   = NaN(nbins,2);
%extract high choices
Hch = x(:,3)==1;
Hch1 = mx(:,3) == 1;
Hch2 = tx(:,3) == 1;
%sample coh vector
cfax  = -100:.1:100;

% make selection array, compute pmf cmf for target and masker data 

%masker (pmf1)
Lcoh  = false(length(mcoh), nbins);
for cc = 1:nbins
    Lcoh(:,cc) = mcoh>=cbins(cc,1) & mcoh<cbins(cc,2);
    cax(cc) = mean(mcoh(Lcoh(:,cc)));
    
    pmf1(cc) = sum(mx(Lcoh(:,cc)&Hch1,3))./sum(Lcoh(:,cc)).*100;
end

%target (pmf2)
Lcoh  = false(length(tcoh), nbins);
for cc = 1:nbins
    Lcoh(:,cc) = tcoh>=cbins(cc,1) & tcoh<cbins(cc,2);
    cax(cc) = mean(tcoh(Lcoh(:,cc)));
    
    pmf2(cc) = sum(tx(Lcoh(:,cc)&Hch2,3))./sum(Lcoh(:,cc)).*100;
end

pmf = {pmf1, pmf2};

%% Begin fitting process - Simultaneous
z = ones(size(x,1), 1);
[fits_,sems_,stats_,preds_,resids_] = logist_fit([z x], 'n');

nLED = 2;
LEDList = [-1; 1];
LEDName = {'Masker','Target'};
allTrialTable = array2table(x,'VariableNames',{'LED','coh','choice'});
choice_perGroup = grpstats(allTrialTable,{'LED','coh'},'mean','dataVars',{'choice'});
coh_space = linspace(-100,100)';
nCoh_sp = length(coh_space);
fake_data = [ones(nCoh_sp*nLED,1),...
    reshape(repmat(LEDList,1,nCoh_sp)',nCoh_sp*nLED,1),...
    repmat(coh_space,nLED,1)];
pred_line = feval(logist_setup([z x], 'n'), fits_ , fake_data);
pred_reshape = reshape(pred_line,nCoh_sp,nLED)';

% myColor = colormap(hsv(nLED)); close gcf;
color = {'r.', 'k.'};
colorl = {'r', 'k'};

fig_name = 'beh_audiDeci_allSubjects_choice_logisticFit_LED';
h = figure('Name',fig_name,'Position',get(0,'ScreenSize'));
p_all = [];
labels = [];
for pp = 1:1:nLED
    cur_cond = choice_perGroup.LED == LEDList(pp);
    plot(cax, pmf{pp}./100, color{pp}, 'MarkerSize', 16); hold on;
    %scatter(cax,pmf{pp}./100,[],myColor(pp,:),'filled')
    hold on;
    p = plot(coh_space,pred_reshape(pp,:),'Color',colorl{pp},'lineWidth',2);
    p_all = [p_all p];
    hold on;
    pse = find(pred_reshape(pp,:)>0.5,1,'first');
    if isempty(pse) || pse == 1
        pse = nan;
    else
        pse = coh_space(pse);
    end
    text(min(coh_space)*0.9,0.6+0.025*pp,['PSE (LED:' LEDName{pp} ' = ' num2str(pse)])
    labels = [labels;{['LED:' LEDName{pp}]}];
    xlim([-100 100])
    xlabel({'coherence';'Low(-) : High(+)'})
    ylabel('% high choice')
    title('Logistic Regression - Simultaneous Fit')
    legend(p_all,labels,'Location','northwest')
end

%% Begin fitting - independent 
zm = ones(size(mx,1), 1);
zt = ones(size(tx,1), 1);

[fits_m,sems_m,stats_m,preds_m,resids_m] = logist_fit([zm mx], 'n');
pred_line_m = feval(logist_setup([zm mx], 'n'), fits_m , fake_data);
pred_reshape_m = reshape(pred_line_m,nCoh_sp,nLED)';

figure();
p = plot(coh_space,pred_reshape_m(1,:),'r','lineWidth',2);
hold on; 
plot(cax,pmf1./100, 'r.', 'MarkerSize', 10);

[fits_t,sems_t,stats_t,preds_t,resids_t] = logist_fit([zt tx], 'n');
pred_line_t = feval(logist_setup([zt tx], 'n'), fits_t , fake_data);
pred_reshape_t = reshape(pred_line_t,nCoh_sp,nLED)';

p2 = plot(coh_space, pred_reshape_t(2,:), 'k', 'LineWidth', 2);
plot(cax, pmf2./100, 'k.', 'MarkerSize', 12);

title('Logistic Regression - Independent Fits'); 
xlabel('Coherence')
ylabel('% High Choice')
legend([p, p2], 'LED with Masker', 'LED with Target', 'Location', 'Southeast');
end