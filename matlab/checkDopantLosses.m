function checkDopantLosses(dopant, T)

c = 3e8; % m/s
h = 6.63e-34; % J*s
kB = 1.3806e-23; % J/K

thermalE = kB*T; % J

ll = (240:1:740)*1e-9;

[~, sigmaabsFun, sigmaemiFun] = getDyeDopantAttributes(dopant);

sigmaabs = sigmaabsFun(ll);
sigmaemi = sigmaemiFun(ll);

lossAverage = h*c*sum((1./ll-1./ll').*sigmaabs.*sigmaemi', 'all')/sum(sigmaabs)/sum(sigmaemi); % J
lossMoment2 = h^2*c^2*sum((1./ll-1./ll').^2.*sigmaabs.*sigmaemi', 'all')/sum(sigmaabs)/sum(sigmaemi); % J^2
lossDeviation = sqrt(lossMoment2-lossAverage^2); % J

fprintf('Loss average = %.4f thermal energies\n', 6*lossAverage/thermalE);
fprintf('Loss deviation = %.4f thermal energies\n', 6*lossDeviation/thermalE);

end