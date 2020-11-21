N = (0.25:0.25:5)*1e22;
L = 0.1;
dopant = "LumogenRed";
D = 1e-3;

lightPoutI = zeros(1, length(N));
lightPoutT = zeros(1, length(N));

for i=1:length(N)
    lightPoutI(i) = oneDopantIterative(dopant, N(i), D, L, 0);
    lightPoutT(i) = oneDopantTruncated(dopant, N(i), D, L, 0);
end

figure(3);
plot(N, lightPoutI);
hold on;
plot(N, lightPoutT);