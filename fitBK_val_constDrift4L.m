function [ps_, rts_] = fitBK_val_constDrift4L(cohs, params, lapse)

% cohs are 0 ... 1. 
%   Assumes values are signed: 
%       + for stim corresponding to correct "A" choices
%       - for stim corresponding to correct "B" choices
%
% 5 parameters: (1 constant + 5(5-1) = 21 total)
%   1   ... k    = drift rate -> CONSTANT 
%   2   ... A    = A bound
%   3   ... B    = B bound
%   4   ... Andt = non-decision time for A choices in msec
%   5   ... Bndt = non-decision time for B choices in msec
%
% lapse is optional

params1 = params(1:5);
params2 = params([1 6 7 8 9]);
params3 = params([1, 10, 11, 12, 13]);
params4 = params([1, 14, 15, 16, 17]);
params5 = params([1, 18, 19, 20, 21]);

mu1 = (params(1)/100000) .* cohs(:,1);
mu2 = (params(1)/100000) .* cohs(:,2);
mu3 = (params(1)/100000) .* cohs(:,3);
if size(cohs, 2) > 3
    mu4 = (params(1)/100000) .* cohs(:,4);
else 
    mu4 = 0;
end 
if size(cohs,2) > 4
    mu5 = (params(1)/100000) .* cohs(:,5);
else 
    mu5 = 0;
end

% scale bounds
A1 = params1(2)/10;
B1 = params1(3)/10;
A2 = params2(2)/10;
B2 = params2(3)/10;
A3 = params3(2)/10;
B3 = params3(3)/10;
A4 = params4(2)/10;
B4 = params4(3)/10;
A5 = params5(2)/10;
B5 = params5(3)/10;

% scale NTD
Andt1 = params1(4);%./1000;
Bndt1 = params1(5);%./1000;
Andt2 = params2(4);%./1000;
Bndt2 = params2(5);%./1000;
Andt3 = params3(4);%./1000;
Bndt3 = params3(5);%./1000;
Andt4 = params4(4);%./1000;
Bndt4 = params4(5);%./1000;
Andt5 = params5(4);%./1000;
Bndt5 = params5(5);%./1000;

% PMF
eA1 = exp(2 .* mu1.* A1);
eB1 = exp(2 .* mu1.* B1);
ps_1 = (eB1 .* eA1 - eA1) ./ (eB1 .* eA1 - 1);
ps_1(abs(cohs(:,1))<=eps) = B1 ./ (A1 + B1);

eA2 = exp(2 .* mu2.* A2);
eB2 = exp(2 .* mu2.* B2);
ps_2 = (eB2 .* eA2 - eA2) ./ (eB2 .* eA2 - 1);
ps_2(abs(cohs(:,2))<=eps) = B2 ./ (A2 + B2);

eA3 = exp(2 .* mu3.* A3);
eB3 = exp(2 .* mu3.* B3);
ps_3 = (eB3 .* eA3 - eA3) ./ (eB3 .* eA3 - 1);
ps_3(abs(cohs(:,3))<=eps) = B3 ./ (A3 + B3);

if size(cohs, 2) > 3

    eA4 = exp(2 .* mu4.* A4);
    eB4 = exp(2 .* mu4.* B4);
    ps_4 = (eB4 .* eA4 - eA4) ./ (eB4 .* eA4 - 1);
    ps_4(abs(cohs(:,4))<=eps) = B4 ./ (A4 + B4);
end

if size(cohs, 2) > 4

    eA5 = exp(2 .* mu5.* A5);
    eB5 = exp(2 .* mu5.* B5);
    ps_5 = (eB5 .* eA5 - eA5) ./ (eB5 .* eA5 - 1);
    ps_5(abs(cohs(:,5))<=eps) = B5 ./ (A5 + B5);
end 

% if lapse given
if nargin > 2
    ps_1 = lapse(1) + (1-2.*lapse(1)).*ps_1;
    ps_2 = lapse(2) + (1-2.*lapse(2)).*ps_2;
    ps_3 = lapse(3) + (1-2.*lapse(3)).*ps_3;
    
    if length(lapse) > 3
        ps_4 = lapse(4) + (1-2.*lapse(4)).*ps_4;
    end 
    if length(lapse) > 4
        ps_5 = lapse(5) + (1-2.*lapse(5)).*ps_5;
    end
end

% CMF
%   E[Ta] = (A+B)/u * coth( (A+B)u/d^2) - B/u * coth(Bu/d^2) for u~=0
%   lim_E[Ta] = (A^2+2AB)/3d^2
%   E[Tb] = (A+B)/u * coth( (A+B)u/d^2) -  A/u * coth(Au/d^2) for u~=0
%   lim_E[Tb] = (B^2+2AB)/3d^2
rts_1 = NaN(size(cohs,1),1);
rts_2 = NaN(size(cohs,1),1);
rts_3 = NaN(size(cohs,1),1);
rts_4 = NaN(size(cohs,1),1);
rts_5 = NaN(size(cohs,1),1);

% positive ivar, T1 choice
Lpt1 = cohs(:,1) > eps;
rts_1(Lpt1) = (A1 + B1) ./ mu1(Lpt1) ./ tanh((A1 + B1) .* mu1(Lpt1)) - ...
    B1 ./ mu1(Lpt1) ./ tanh(B1 .* mu1(Lpt1)) + Andt1;
Lpt1 = cohs(:,2) > eps;
rts_2(Lpt1) = (A2 + B2) ./ mu2(Lpt1) ./ tanh((A2 + B2) .* mu2(Lpt1)) - ...
    B2 ./ mu2(Lpt1) ./ tanh(B2 .* mu2(Lpt1)) + Andt2;
Lpt1 = cohs(:,3) > eps;
rts_3(Lpt1) = (A3 + B3) ./ mu3(Lpt1) ./ tanh((A3 + B3) .* mu3(Lpt1)) - ...
    B3 ./ mu3(Lpt1) ./ tanh(B3 .* mu3(Lpt1)) + Andt3;
if size(cohs, 2) > 3
    Lpt1 = cohs(:,4) > eps;
    rts_4(Lpt1) = (A4 + B4) ./ mu4(Lpt1) ./ tanh((A4 + B4) .* mu4(Lpt1)) - ...
        B4 ./ mu4(Lpt1) ./ tanh(B4 .* mu4(Lpt1)) + Andt4;
end 
if size(cohs, 2) > 4
    Lpt1 = cohs(:,5) > eps;
    rts_5(Lpt1) = (A5 + B5) ./ mu5(Lpt1) ./ tanh((A5 + B5) .* mu5(Lpt1)) - ...
        B5 ./ mu5(Lpt1) ./ tanh(B5 .* mu5(Lpt1)) + Andt5;
end

% zero ivar, T1 choice
L0t1 = cohs(:,1) >= 0 & cohs(:,1) <= eps; 
rts_1(L0t1) = (A1.^2 + 2 .* A1 .* B1) ./ 3 + Andt1;
L0t1 = cohs(:,2) >= 0 & cohs(:,2) <= eps; 
rts_2(L0t1) = (A2.^2 + 2 .* A2 .* B2) ./ 3 + Andt2;
L0t1 = cohs(:,3) >= 0 & cohs(:,3) <= eps; 
rts_3(L0t1) = (A3.^2 + 2 .* A3 .* B3) ./ 3 + Andt3;
if size(cohs, 2) > 3
    L0t1 = cohs(:,4) >= 0 & cohs(:,4) <= eps; 
    rts_4(L0t1) = (A4.^2 + 2 .* A4 .* B4) ./ 3 + Andt4;
end 
if size(cohs, 2) > 4
    L0t1 = cohs(:,5) >= 0 & cohs(:,5) <= eps; 
    rts_5(L0t1) = (A5.^2 + 2 .* A5 .* B5) ./ 3 + Andt5;
end 

% negative ivar, T2 choice
Lnt2 = cohs(:,1) < -eps;
rts_1(Lnt2) = (A1 + B1) ./ mu1(Lnt2) ./ tanh((A1 + B1) .* mu1(Lnt2)) - ...
    A1 ./ mu1(Lnt2) ./ tanh(A1 .* mu1(Lnt2)) + Bndt1;
Lnt2 = cohs(:,2) < -eps;
rts_2(Lnt2) = (A2 + B2) ./ mu2(Lnt2) ./ tanh((A2 + B2) .* mu2(Lnt2)) - ...
    A2 ./ mu2(Lnt2) ./ tanh(A2 .* mu2(Lnt2)) + Bndt2;
Lnt2 = cohs(:,3) < -eps;
rts_3(Lnt2) = (A3 + B3) ./ mu3(Lnt2) ./ tanh((A3 + B3) .* mu3(Lnt2)) - ...
    A3 ./ mu3(Lnt2) ./ tanh(A3 .* mu3(Lnt2)) + Bndt3;
if size(cohs,2) > 3
    Lnt2 = cohs(:,4) < -eps;
    rts_4(Lnt2) = (A4 + B4) ./ mu4(Lnt2) ./ tanh((A4 + B4) .* mu4(Lnt2)) - ...
        A4 ./ mu4(Lnt2) ./ tanh(A4 .* mu4(Lnt2)) + Bndt4;
end 
if size(cohs,2) > 4
    Lnt2 = cohs(:,5) < -eps;
    rts_5(Lnt2) = (A5 + B5) ./ mu5(Lnt2) ./ tanh((A5 + B5) .* mu5(Lnt2)) - ...
        A5 ./ mu5(Lnt2) ./ tanh(A5 .* mu5(Lnt2)) + Bndt5;
end 

% zero ivar, T2 choice
L0t2 = cohs(:,1) <= 0 & cohs(:,1) >= -eps; 
rts_1(L0t2) = (B1.^2 + 2 .* A1 .* B1) ./ 3 + Bndt1;
L0t2 = cohs(:,2) <= 0 & cohs(:,2) >= -eps; 
rts_2(L0t2) = (B2.^2 + 2 .* A2 .* B2) ./ 3 + Bndt2;
L0t2 = cohs(:,3) <= 0 & cohs(:,3) >= -eps; 
rts_3(L0t2) = (B3.^2 + 2 .* A3 .* B3) ./ 3 + Bndt3;
if size(cohs, 2) > 3
    L0t2 = cohs(:,4) <= 0 & cohs(:,4) >= -eps; 
    rts_4(L0t2) = (B4.^2 + 2 .* A4 .* B4) ./ 3 + Bndt4;
end 
if size(cohs, 2) > 4
    L0t2 = cohs(:,5) <= 0 & cohs(:,5) >= -eps; 
    rts_5(L0t2) = (B5.^2 + 2 .* A5 .* B5) ./ 3 + Bndt5;
end 

ps_ = [ps_1, ps_2, ps_3]; 
rts_ = [rts_1, rts_2, rts_3];
if size(cohs, 2) > 3
    ps_ = [ps_, ps_4];
    rts_ = [rts_, rts_4];
end 
if size(cohs, 2) > 4
    ps_ = [ps_, ps_5];
    rts_ = [rts_, rts_5];
end 