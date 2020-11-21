function [tau, sigmaabs, sigmaemi, quenchingA, quenchingB] = getDyeDopantAttributes(dopant)

if strcmp(dopant, 'RhB')
    tau = 4.8e-9;
    sigmaabs = @sigmaabs_RhB;
    sigmaemi = @sigmaemi_RhB;
    quenchingA = 0; % Not known -> no quenching happens
    quenchingB = -3;
elseif strcmp(dopant, 'C1')
    tau = 3.4e-9;
    sigmaabs = @sigmaabs_C1;
    sigmaemi = @sigmaemi_C1;
    quenchingA = 0; % Not known -> no quenching happens
    quenchingB = -3;
elseif strcmp(dopant, 'C6')
    tau = 3.4e-9; % Need to check again
    sigmaabs = @sigmaabs_C6;
    sigmaemi = @sigmaemi_C6;
    quenchingA = 0; % Not known -> no quenching happens
    quenchingB = -3;
elseif strcmp(dopant, 'LumogenRed')
    tau = 6e-9;
    sigmaabs = @sigmaabs_LumogenRed;
    sigmaemi = @sigmaemi_LumogenRed;
    quenchingA = 7.3e9;
    quenchingB = -3.92;
else % Default is Rh6G
    tau = 4.8e-9;
    sigmaabs = @sigmaabs_Rh6G;
    sigmaemi = @sigmaemi_Rh6G;
    quenchingA = 7e9;
    quenchingB = -2.4;
end

end