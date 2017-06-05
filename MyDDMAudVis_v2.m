function MyDDMAudVis_v2(block_size)

% Data File being used:
% 4th column: percentage of high tone (mapped this percentage to signed coherence like this: (proportion-50).*2)
% 8th column: monkey's choice (2: high-tone choice, 1: low-tone choice)
% 13th column: response time (sec)
% 9th column: success of choice (1: correct, 0: incorrect)

close all

%cd into file that holds data 
cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/Data/');
%load in data, starting with 3rd column (subtract 2 from all indices
%indicated above) 
PopBehavior = csvread('AudVisTask_v1_Beta_Diana_170602_1142.csv', 1, 2); %block size 75
Headings = load('AudVisTask_v1_Beta_Diana_170602_1142_table.mat');
h = Headings.data_table_stim(:, 2);

%extract block visual modes from matrix
block1 = h{1, 1};
block2 = h{block_size + 1,1};
block3 = h{2*block_size + 1,1};
block4 = h{3*block_size + 1,1};

%coherence bins
% cbins = [ ...
%     -100  -99
%     -99  -60
%     -60  -20
%     -20    20
%     20   60
%     60  100
%     100  101];

cbins = [ ...
    -100  -99
     -99  -50
     -50  -34
     -34  -20
     -20    0
       0   20
      20   34
      34   50
      50   99
      99  100];

%calculate number of coherence bins 
nbins = size(cbins,1);
%calculate average value of each coherence bin - output vector 
cax   = mean(cbins,2);
%create vector from -100 to 100
cfax  = -100:.1:100;
%initialize vectors to hold pmf and cmf values
pmf1   = NaN(nbins,1);
cmf1   = NaN(nbins,2);
pmf2   = NaN(nbins,1);
cmf2   = NaN(nbins,2);
pmf3   = NaN(nbins,1);
cmf3   = NaN(nbins,2);
pmf4   = NaN(nbins,1);
cmf4   = NaN(nbins,2);

%set data to variable
Behavior1 = PopBehavior(1:block_size, :);
Behavior2 = PopBehavior(block_size+1:2*block_size, :);
Behavior3 = PopBehavior(2*block_size+1:3*block_size, :);
Behavior4 = PopBehavior(3*block_size+1:4*block_size, :);
%enter file for fit functions
cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/DDM/');

%calculate number of trials 
ntrials = block_size;
%create vector of true = high false = low; ensures all values are 0 or 1 
Lch1     = Behavior1(:,6) == 2;
Lch2     = Behavior2(:,6) == 2;
Lch3     = Behavior3(:,6) == 2;
Lch4     = Behavior4(:,6) == 2;

% create a data matrix for each visual mode, columns are:
%   1. re-scale signal strength to [-100, 100]
%   2. choice (1: high-tone choice, 0: low-tone choice)
%   3. RT (ms)
%   4. correct (1) or error (0)
scoh1 = Behavior1(:,2).*200-100;
Lcor1 = Behavior1(:,7); %(Behavior(:,1)>= -20 & Behavior(:,1)<20) | (Behavior(:,1)>50&Lch) | (Behavior(:,1)<50&~Lch);
data1 = cat(2, scoh1, Behavior1(:,6) - 1, Behavior1(:,11), double(Lcor1));
scoh2 = Behavior2(:,2).*200-100;
Lcor2 = Behavior2(:,7); 
data2 = cat(2, scoh2, Behavior2(:,6) - 1, Behavior2(:,11), double(Lcor2));
scoh3 = Behavior3(:,2).*200-100;
Lcor3 = Behavior3(:,7); 
data3 = cat(2, scoh3, Behavior3(:,6) - 1, Behavior3(:,11), double(Lcor3));
scoh4 = Behavior4(:,2).*200-100;
Lcor4 = Behavior4(:,7); 
data4 = cat(2, scoh4, Behavior4(:,6) - 1, Behavior4(:,11), double(Lcor4));

% make selection array, compute pmf and cmf
% PMF - a function that gives the probability that a discrete random variable is exactly equal to some value
% CMF -  the probability that X will take a value less than or equal to x
Lcoh1  = false(ntrials, nbins);
for cc = 1:nbins
    Lcoh1(:,cc) = scoh1>=cbins(cc,1) & scoh1<cbins(cc,2);
    cax1(cc) = mean(scoh1(Lcoh1(:,cc)));
    
    pmf1(cc) = sum(Behavior1(Lcoh1(:,cc)&Lch1,2))./sum(Lcoh1(:,cc)).*100;
    
    cmf1(cc,1) = nanmean(data1(Lcoh1(:,cc)& Lch1,3));
    cmf1(cc,2) = nanmean(data1(Lcoh1(:,cc)&~Lch1,3));
end

Lcoh2  = false(ntrials, nbins);
for cc = 1:nbins
    Lcoh2(:,cc) = scoh2>=cbins(cc,1) & scoh2<cbins(cc,2);
    cax2(cc) = mean(scoh2(Lcoh2(:,cc)));
    
    pmf2(cc) = sum(Behavior2(Lcoh2(:,cc)&Lch2,2))./sum(Lcoh2(:,cc)).*100;
    
    cmf2(cc,1) = nanmean(data2(Lcoh2(:,cc)& Lch2,3));
    cmf2(cc,2) = nanmean(data2(Lcoh2(:,cc)&~Lch2,3));
end

Lcoh3  = false(ntrials, nbins);
for cc = 1:nbins
    Lcoh3(:,cc) = scoh3>=cbins(cc,1) & scoh3<cbins(cc,2);
    cax3(cc) = mean(scoh3(Lcoh3(:,cc)));
    
    pmf3(cc) = sum(Behavior3(Lcoh3(:,cc)&Lch3,2))./sum(Lcoh3(:,cc)).*100;
    
    cmf3(cc,1) = nanmean(data3(Lcoh3(:,cc)& Lch3,3));
    cmf3(cc,2) = nanmean(data3(Lcoh3(:,cc)&~Lch3,3));
end

Lcoh4  = false(ntrials, nbins);
for cc = 1:nbins
    Lcoh4(:,cc) = scoh4>=cbins(cc,1) & scoh4<cbins(cc,2);
    cax4(cc) = mean(scoh4(Lcoh4(:,cc)));
    
    pmf4(cc) = sum(Behavior4(Lcoh4(:,cc)&Lch4,2))./sum(Lcoh4(:,cc)).*100;
    
    cmf4(cc,1) = nanmean(data4(Lcoh4(:,cc)& Lch4,3));
    cmf4(cc,2) = nanmean(data4(Lcoh4(:,cc)&~Lch4,3));
end

X0  = [200 200 200 200 200];
Xlb = [0.01 0.01 0.01 0.01 0.01];
Xub = [50000 50000 50000 2000 2000];

Llapse1 = abs(scoh1)>=90;
lapse1 = 1-sum(Llapse1&Lcor1)./sum(Llapse1);
Llapse2 = abs(scoh2)>=90;
lapse2 = 1-sum(Llapse2&Lcor2)./sum(Llapse2);
Llapse3 = abs(scoh3)>=90;
lapse3 = 1-sum(Llapse3&Lcor3)./sum(Llapse3);
Llapse4 = abs(scoh4)>=90;
lapse4 = 1-sum(Llapse4&Lcor4)./sum(Llapse4);

% get err from initial fit
err0_1 = fitJT_err(X0, data1, lapse1);
err0_2 = fitJT_err(X0, data2, lapse2);
err0_3 = fitJT_err(X0, data3, lapse3);
err0_4 = fitJT_err(X0, data4, lapse4);

% fit it using pattern search
[fits1,err1] = patternsearch(@(x)fitJT_err(x, data1, lapse1), ...
    X0, [], [], [], [], Xlb, Xub, [], ...
    psoptimset('MaxIter', 5000, 'MaxFunEvals', 5000));
[fits2,err2] = patternsearch(@(x)fitJT_err(x, data2, lapse2), ...
    X0, [], [], [], [], Xlb, Xub, [], ...
    psoptimset('MaxIter', 5000, 'MaxFunEvals', 5000));
[fits3,err3] = patternsearch(@(x)fitJT_err(x, data3, lapse3), ...
    X0, [], [], [], [], Xlb, Xub, [], ...
    psoptimset('MaxIter', 5000, 'MaxFunEvals', 5000));
[fits4,err4] = patternsearch(@(x)fitJT_err(x, data4, lapse4), ...
    X0, [], [], [], [], Xlb, Xub, [], ...
    psoptimset('MaxIter', 5000, 'MaxFunEvals', 5000));

% possibly seed as init values
if err1 < err0_1
    X0g1  = fits1;
    err0_1 = err1;
else
    X0g1 = X0;
end

if err2 < err0_2
    X0g2  = fits2;
    err0_2 = err2;
else
    X0g2 = X0;
end

if err3 < err0_3
    X0g3  = fits3;
    err0_3 = err3;
else
    X0g3 = X0;
end

if err4 < err0_4
    X0g4  = fits4;
    err0_4 = err4;
else
    X0g4 = X0;
end

% now gradient descent
[fitsg1,errg1] = fmincon(@(x)fitJT_err(x, data1, lapse1), ...
    X0g1, [], [], [], [], Xlb, Xub, [], ...
    optimset('Algorithm', 'active-set', ...
    'MaxIter', 30000, 'MaxFunEvals', 30000));%, 'Display', 'iter'));
[fitsg2,errg2] = fmincon(@(x)fitJT_err(x, data2, lapse2), ...
    X0g2, [], [], [], [], Xlb, Xub, [], ...
    optimset('Algorithm', 'active-set', ...
    'MaxIter', 30000, 'MaxFunEvals', 30000));%, 'Display', 'iter'));
[fitsg3,errg3] = fmincon(@(x)fitJT_err(x, data3, lapse3), ...
    X0g3, [], [], [], [], Xlb, Xub, [], ...
    optimset('Algorithm', 'active-set', ...
    'MaxIter', 30000, 'MaxFunEvals', 30000));%, 'Display', 'iter'));
[fitsg4,errg4] = fmincon(@(x)fitJT_err(x, data4, lapse4), ...
    X0g4, [], [], [], [], Xlb, Xub, [], ...
    optimset('Algorithm', 'active-set', ...
    'MaxIter', 30000, 'MaxFunEvals', 30000));%, 'Display', 'iter'));

% save the best between the two methods
if errg1 < err0_1
    fits1 = fitsg1;
end
if errg2 < err0_2
    fits2 = fitsg2;
end
if errg3 < err0_3
    fits3 = fitsg3;
end
if errg4 < err0_4
    fits4 = fitsg4;
end

%%% PLOTZ
[ps1,rts1] = fitJT_val_simple5L(cfax, fits1, lapse1);
[ps2,rts2] = fitJT_val_simple5L(cfax, fits2, lapse2);
[ps3,rts3] = fitJT_val_simple5L(cfax, fits3, lapse3);
[ps4,rts4] = fitJT_val_simple5L(cfax, fits4, lapse4);

subplot(3,1,1); cla reset; hold on;
%%%block 1
plot(cax1, pmf1, 'k.', 'MarkerSize', 8);
p1 = plot(cfax, ps1.*100, 'k-');
plot([-100 0], [lapse1*100 lapse1*100], 'k--');
plot([0 100], [100-lapse1*100 100-lapse1*100], 'k--');
%%%block 2
plot(cax2, pmf2, 'r.', 'MarkerSize', 8);
p2 = plot(cfax, ps2.*100, 'r-');
plot([-100 0], [lapse2*100 lapse2*100], 'r--');
plot([0 100], [100-lapse2*100 100-lapse2*100], 'r--');
%%%block 3
plot(cax3, pmf3, 'b.', 'MarkerSize', 8); hold on;
p3 = plot(cfax, ps3.*100, 'b-'); hold on;
plot([-100 0], [lapse3*100 lapse3*100], 'b--'); hold on;
plot([0 100], [100-lapse3*100 100-lapse3*100], 'b--');
%%%block 4
plot(cax4, pmf4, 'g.', 'MarkerSize', 8);
p4 = plot(cfax, ps4.*100, 'g-');
plot([-100 0], [lapse4*100 lapse4*100], 'g--');
plot([0 100], [100-lapse4*100 100-lapse4*100], 'g--');
xlabel('Coherence (%): +100 means all high tones')
ylabel('high-tone choice (%)')
legend([p1, p2, p3, p4], [block1, block2, block3, block4])

%Horizontal shifts of these lines imply changes in the mean rate-of-rise 
%swivels about a fixed point at infinite RT imply changes in the bound height
subplot(3,1,2); cla reset; hold on;
%%%block 1
plot(cax1(cax1>=0), cmf1(cax1>=0,1), 'k.', 'MarkerSize', 8);
plot(cax1(cax1<=0), cmf1(cax1<=0,2), 'k.', 'MarkerSize', 8);
g1 = plot(cfax, rts1, 'k-');
plot([0 100], fits1([4 4]), 'k--');
plot([-100 0], fits1([5 5]), 'k--');
%%%block 2
plot(cax2(cax2>=0), cmf2(cax2>=0,1), 'r.', 'MarkerSize', 8);
plot(cax2(cax2<=0), cmf2(cax2<=0,2), 'r.', 'MarkerSize', 8);
g2 = plot(cfax, rts2, 'r-');
plot([0 100], fits2([4 4]), 'r--');
plot([-100 0], fits2([5 5]), 'r--');
%%%block 3
plot(cax3(cax3>=0), cmf3(cax3>=0,1), 'b.', 'MarkerSize', 8);
plot(cax3(cax3<=0), cmf3(cax3<=0,2), 'b.', 'MarkerSize', 8);
g3 = plot(cfax, rts3, 'b-');
plot([0 100], fits3([4 4]), 'b--');
plot([-100 0], fits3([5 5]), 'b--');
%%%block 4
plot(cax4(cax4>=0), cmf4(cax4>=0,1), 'g.', 'MarkerSize', 8);
plot(cax4(cax4<=0), cmf4(cax4<=0,2), 'g.', 'MarkerSize', 8);
g4 = plot(cfax, rts4, 'g-');
plot([0 100], fits4([4 4]), 'g--');
plot([-100 0], fits4([5 5]), 'g--');
xlabel('Coherence (%): +100 means all high tones')
ylabel('Response time (ms)')
legend([g1, g2, g3, g4], [block1, block2, block3, block4])

subplot(3,1,3); cla reset; hold on;
%%%block 1
plot(cfax(cfax<0), rts1(cfax<0)-fits1(5), 'k-');hold on;
v1 = plot(cfax(cfax>0), rts1(cfax>0)-fits1(4), 'k-');
%%%block 2
plot(cfax(cfax<0), rts2(cfax<0)-fits2(5), 'r-');hold on;
v2 = plot(cfax(cfax>0), rts2(cfax>0)-fits2(4), 'r-');
%%%block 3
plot(cfax(cfax<0), rts3(cfax<0)-fits3(5), 'b-');hold on;
v3 = plot(cfax(cfax>0), rts3(cfax>0)-fits3(4), 'b-');
%%%block 4
plot(cfax(cfax<0), rts4(cfax<0)-fits4(5), 'g-');hold on;
v4 = plot(cfax(cfax>0), rts4(cfax>0)-fits4(4), 'g-');
xlabel('Coherence (%): +100 means all high tones')
ylabel('Decision time (ms): RT-nonDT')
legend([v1 v2 v3 v4],[block1 block2 block3 block4])
end