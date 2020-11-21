L = [0.1 0.2 0.4 0.5 0.7 1 1.5 2]; % m

lightPout = zeros(1, length(L));
electricPout = zeros(1, length(L));

dopant = "LumogenRed";
N = 1.776e23; % m^-3 (0.03 mol%)
diameter = 1e-3; % m
darkL = .03; % m

for i = 1:length(L)
    [lightPout(i), electricPout(i)] = oneDopantIterative(dopant, N, diameter, L(i), darkL);
end

figure(1);
plot(L, lightPout*1e6);

% figure(2);
% plot(L, electricPout*1e6);