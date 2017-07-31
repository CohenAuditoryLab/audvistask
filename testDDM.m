[fits1, fits2, fits3, fits4, fits5] = MyDDMAudVis_v2(80);

mu1 = fits1(1)/100000 .* cohs;
A1 = fits1(2);
B1 = fits1(3);
Andt1 = fits1(4);
Bndt1 = fits1(5);

cohs = -100:.1:100;
scoh = 2.*((cohs.*100)-50);
Lcor1 = (cohs>= -.20 & cohs<.20) | (cohs>.50&true) | (cohs<.50&~true);

eA1 = exp(2 .* mu1.* A1);
eB1 = exp(2 .* mu1.* B1);
ps_1 = (eB1 .* eA1 - eA1) ./ (eB1 .* eA1 - 1);
ps_1(abs(cohs)<=eps) = B1 ./ (A1 + B1);

data = cat(2, cohs, ps_1, ones(1, 2001), double(Lcor1));

Llapse = abs(cohs)>=90;
lapse = 1-sum(Llapse&Lcor1)./sum(Llapse);

 cbins = [ ...
    -100  -99
     -99  -80
     -80  -65
     -65  -50
     -50  -34
     -34  -20
     -20  -10
     -10    0
       0   10
      10   20
      20   34
      34   50
      50   65
      65   80
      80   99
      99  100];

%calculate number of coherence bins
nbins = size(cbins,1);

Lcoh1  = false(2001, nbins);
for cc = 1:nbins
    Lcoh1(:,cc) = scoh>=cbins(cc,1) & scoh<cbins(cc,2);
    cax1(cc) = mean(scoh(Lcoh1(:,cc)));
    
    pmf1(cc) = sum(cohs(Lcoh1(:,cc)&true))./sum(Lcoh1(:,cc)).*100;
    
    cmf1(cc,1) = nanmean(data(Lcoh1(:,cc)& true));
    cmf1(cc,2) = nanmean(data(Lcoh1(:,cc)&~ true));
end

X0  = [200 200 200 200 200];
Xlb = [0.001 0.001 0.001 0.001 0.001];
Xub = [50000 50000 50000 200 200];

err0 = fitJT_err(X0, data, lapse);

[fits,err1] = patternsearch(@(x)fitJT_err(x, data, lapse), ...
    X0, [], [], [], [], Xlb, Xub, [], ...
    psoptimset('MaxIter', 50000, 'MaxFunEvals', 50000));

if err1 < err0
    X0g1  = fits;
    err0 = err1;
else
    X0g1 = X0;
end

[fitsg,errg1] = fmincon(@(x)fitJT_err(x, data, lapse), ...
    X0g1, [], [], [], [], Xlb, Xub, [], ...
    optimset('Algorithm', 'active-set', ...
    'MaxIter', 40000, 'MaxFunEvals', 50000));

if errg1 < err0
    fits = fitsg;
end

[ps,rts] = fitJT_val_simple5L(cohs, fits1, lapse);

figure();
plot(cohs, ps_1, 'r.', 'MarkerSize', 12); hold on;
plot(cohs, ps.*100, 'k-', 'LineWidth', .75);

