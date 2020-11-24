clear;
close all;

ll = (300:2:800)*1e-9;

ncore = refractionIndexPMMA(ll);
sigmaabs = sigmaabs_C6(ll);
beta = (ncore - 1)./(2*ncore);
sigmaemi = sigmaemi_C6(ll);
wnsp = sigmaemi/sum(sigmaemi);
alfaPMMA = valuesalfaPMMA(ll);

lpump = 450e-9;
lsigma = 10e-9;
pumpDistribution = exp(-((ll-lpump)/lsigma).^2);
pumpDistribution = pumpDistribution/sum(pumpDistribution);
pumpPower = 2; % W
Ppump = pumpPower*pumpDistribution;

lightL = 0.01; % m
% darkL = 0; % m
N = 3e22; % m^-3
Nppm = N*1200/1.18e6*1e6/6.023e23; % ppm

M = length(ll);

C = beta.*wnsp./ll;

G1 = zeros(M);
for i = 1:M
    for j = 1:M
        if(i == j)
            G1(i, i) = -sigmaabs(i)*N*(1-beta(i)*wnsp(i))-alfaPMMA(i)+realmin^(1/3);
        else
            G1(j, i) = sigmaabs(i)*N*beta(j)*wnsp(j)*ll(i)/ll(j);
        end
    end
end

G2 = zeros(M);
for i = 1:M
    G2(i, i) = -sigmaabs(i)*N-alfaPMMA(i)+sqrt(realmin);
end

% result = expm(G*darkL)*(expm(G*lightL)-eye(M))*(G\C');
result1 = (expm(G1*lightL)-eye(M))*(G1\C');
result1 = result1/sum(result1);

result2 = (expm(G2*lightL)-eye(M))*(G2\C');
result2 = result2/sum(result2);

result3 = expm(G1*lightL)*Ppump';
result3 = result3/sum(result3);

result4 = expm(G2*lightL)*Ppump';
result4 = result4/sum(result4);

figure(1);
plot(ll, result1);
hold on;
plot(ll, result2);

figure(2);
plot(ll, result3);
hold on;
plot(ll, result4);