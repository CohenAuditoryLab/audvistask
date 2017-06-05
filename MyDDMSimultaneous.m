function MyDDMSimultaneous(block_size)
%% Set up

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
block5 = h{4*block_size + 1, 1};

%coherence bins
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

%set data to variable
Behavior1 = PopBehavior(1:block_size, :);
Behavior2 = PopBehavior(block_size+1:2*block_size, :);
Behavior3 = PopBehavior(2*block_size+1:3*block_size, :);
Behavior4 = PopBehavior(3*block_size+1:4*block_size, :);
Behavior5 = PopBehavior(4*block_size+1:5*block_size, :);

%enter file for fit functions
cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/');

%calculate number of trials 
ntrials = block_size;

%% Data formatting

%create vector of true = high false = low; ensures all values are 0 or 1 
Lch1     = Behavior1(:,6) == 2;
Lch2     = Behavior2(:,6) == 2;
Lch3     = Behavior3(:,6) == 2;
Lch4     = Behavior4(:,6) == 2;
Lch5     = Behavior5(:,6) == 2;

% create a data matrix for each visual mode, columns are:
%   1. re-scale signal strength to [-100, 100]
%   2. choice (1: high-tone choice, 0: low-tone choice)
%   3. RT (ms)
%   4. correct (1) or error (0)

%these should act as parameter input into fit function
scoh1 = Behavior1(:,2).*200-100;
Lcor1 = Behavior1(:,7);
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
scoh5 = Behavior5(:,2).*200-100;
Lcor5 = Behavior5(:,7);
data5 = cat(2, scoh5, Behavior5(:,6) - 1, Behavior5(:,11), double(Lcor5));

%define x0, lower and upper bounds
X0  = [200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200, ...
    200 200 200 200];
Xlb = [0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01...
    0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
Xub = [50000 50000 50000 2000 2000 50000 50000 50000 2000 2000 50000, ...
    50000 50000 2000 2000 50000 50000 50000 2000 2000 50000];

%calculate lapse
Llapse1 = abs(scoh1)>=90;
lapse1 = 1-sum(Llapse1&Lcor1)./sum(Llapse1);
Llapse2 = abs(scoh2)>=90;
lapse2 = 1-sum(Llapse2&Lcor2)./sum(Llapse2);
Llapse3 = abs(scoh3)>=90;
lapse3 = 1-sum(Llapse3&Lcor3)./sum(Llapse3);
Llapse4 = abs(scoh4)>=90;
lapse4 = 1-sum(Llapse4&Lcor4)./sum(Llapse4);
Llapse5 = abs(scoh5)>=90;
lapse5 = 1-sum(Llapse5&Lcor5)./sum(Llapse5);

%% Independent Curve Fits

%fit the curves independently (full model)

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

Lcoh5 = false(ntrials, nbins);
for cc = 1:nbins
    Lcoh5(:,cc) = scoh5>=cbins(cc,1) & scoh5<cbins(cc,2);
    cax5(cc) = mean(scoh5(Lcoh5(:,cc)));
    
    pmf5(cc) = sum(Behavior5(Lcoh5(:,cc)&Lch5,2))./sum(Lcoh5(:,cc)).*100;
    
    cmf5(cc,1) = nanmean(data5(Lcoh5(:,cc)& Lch5,3));
    cmf5(cc,2) = nanmean(data5(Lcoh5(:,cc)&~Lch5,3));
end

% get err from initial fit
err0_1 = fitJT_err(X0, data1, lapse1);
err0_2 = fitJT_err(X0, data2, lapse2);
err0_3 = fitJT_err(X0, data3, lapse3);
err0_4 = fitJT_err(X0, data4, lapse4);
err0_5 = fitJT_err(X0, data5, lapse5);

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
[fits5,err5] = patternsearch(@(x)fitJT_err(x, data5, lapse5), ...
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

if err5 < err0_5
    X0g5  = fits5;
    err0_5 = err5;
else
    X0g5 = X0;
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
[fitsg5,errg5] = fmincon(@(x)fitJT_err(x, data5, lapse5), ...
    X0g5, [], [], [], [], Xlb, Xub, [], ...
    optimset('Algorithm', 'active-set', ...
    'MaxIter', 30000, 'MaxFunEvals', 30000));%, 'Display', 'iter'));

% save the best fit between the two methods
if errg1 < err0_1
    fits1 = fitsg1;
    err0_1 = errg1;
end
if errg2 < err0_2
    fits2 = fitsg2;
    err0_2 = errg2;
end
if errg3 < err0_3
    fits3 = fitsg3;
    err0_3 = errg3;
end
if errg4 < err0_4
    fits4 = fitsg4;
    err0_4 = errg4;
end
if errg5 < err0_5
    fits5 = fitsg5;
    err0_5 = errg5;
end

%total error 
err_indep = -(err0_1 + err0_2 + err0_3 + err0_4 + err0_5);

%% Simultaneous curve fit, holding accumulation rate parameter constant

data_all = cat(2, data1, data2, data3, data4, data5);

%fit holding mu constant
[fits_mu,err_mu] = fmincon(@(fits)fitBK_err_constDrift(fits, ...
    data_all, [lapse1, lapse2, lapse3, lapse4, lapse5]), ...
    X0g1, [], [], [], [], Xlb, Xub, [], optimset('Algorithm', 'active-set', ...
    'MaxIter', 30000, 'MaxFunEvals', 30000));

%% Simultaneous curve fit, holding bounds parameters constant

%fit holding both A and B bounds constant
[fits_AB,err_AB] = fmincon(@(fits)fitBK_err_constDrift(fits, ...
    data_all, [lapse1, lapse2, lapse3, lapse4, lapse5]), ...
    X0g1, [], [], [], [], Xlb, Xub, [], optimset('Algorithm', 'active-set', ...
    'MaxIter', 30000, 'MaxFunEvals', 30000));

%% Determine the best of the models 
%use BIC or AIC to determine the best fit model

errors = [err_indep, err_mu, err_AB];
aic = aicbic(errors, [5, 4, 3]);
[~, index] = min(aic);

if index == 1
    [ps1,rts1] = fitJT_val_simple5L(cfax, fits1, lapse1);
    [ps2,rts2] = fitJT_val_simple5L(cfax, fits2, lapse2);
    [ps3,rts3] = fitJT_val_simple5L(cfax, fits3, lapse3);
    [ps4,rts4] = fitJT_val_simple5L(cfax, fits4, lapse4);
    [ps5,rts5] = fitJT_val_simple5L(cfax, fits5, lapse5);
elseif index == 2 
    cohs  = linspace(-100, 100, length(data1(:, 1)));
    [ps, rts] = fitBK_val_constDrift4L(cohs, fits_mu, [lapse1, lapse2,...
        lapse3, lapse4, lapse5]);
elseif index == 3
    
end
    
    

end 