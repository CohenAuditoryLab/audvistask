function [ps_, rts_] = fitJT_val_simple5L(cohs, params, lapse)
%
% L for Lapse (optional)
%
% cohs are 0 ... 1. 
%   Assumes values are signed: 
%       + for stim corresponding to correct "A" choices
%       - for stim corresponding to correct "B" choices
%
% 5 parameters:
%   1   ... k    = drift rate (scaled by coh)
%   2   ... A    = A bound
%   3   ... B    = B bound
%   4   ... Andt = non-decision time for A choices in msec
%   5   ... Bndt = non-decision time for B choices in msec
%
% lapse is optional

% Drift rate is linear function of ivar, scaled
mu = (params(1)/100000) .* cohs;

% scale bounds
A = params(2)/10;
B = params(3)/10;

% scale NTD
Andt = params(4);%./1000;
Bndt = params(5);%./1000;

% PMF
%   eA = exp(2uA) 
%   eB = exp(2uB)
%   pA = (eB-1)/(eB-1/eA)) for u~=0
%   lim_pA = B/(A+B) for u=0
eA = exp(2 .* mu.* A);
eB = exp(2 .* mu.* B);
ps_ = (eB .* eA - eA) ./ (eB .* eA - 1);
ps_(abs(cohs)<=eps) = B ./ (A + B);

% if lapse given
if nargin > 2
    ps_ = lapse + (1-2.*lapse).*ps_;
end

% CMF
%   E[Ta] = (A+B)/u * coth( (A+B)u/d^2) - B/u * coth(Bu/d^2) for u~=0
%   lim_E[Ta] = (A^2+2AB)/3d^2
%   E[Tb] = (A+B)/u * coth( (A+B)u/d^2) -  A/u * coth(Au/d^2) for u~=0
%   lim_E[Tb] = (B^2+2AB)/3d^2
rts_ = NaN(size(cohs,1),1);

% positive ivar, T1 choice
Lpt1 = cohs > eps;
rts_(Lpt1) = (A + B) ./ mu(Lpt1) ./ tanh((A + B) .* mu(Lpt1)) - ...
    B ./ mu(Lpt1) ./ tanh(B .* mu(Lpt1)) + Andt;

% zero ivar, T1 choice
L0t1 = cohs >= 0 & cohs <= eps; 
rts_(L0t1) = (A.^2 + 2 .* A .* B) ./ 3 + Andt;

% negative ivar, T2 choice
Lnt2 = cohs < -eps;
rts_(Lnt2) = (A + B) ./ mu(Lnt2) ./ tanh((A + B) .* mu(Lnt2)) - ...
    A ./ mu(Lnt2) ./ tanh(A .* mu(Lnt2)) + Bndt;

% zero ivar, T2 choice
L0t2 = cohs <= 0 & cohs >= -eps; 
rts_(L0t2) = (B.^2 + 2 .* A .* B) ./ 3 + Bndt;

% disp(fits)