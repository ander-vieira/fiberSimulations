function absorption = fiberAbsorptionTwoInterfaces(ncore, nclad, diameter, q, alfaDopant, alfaCore, alfaClad)
%FIBERABSORPTIONTWOINTERFACES Absorbed fraction model without reflections
%   Calculates the fraction of incident sunlight absorbed by a fiber
%   with a two interface model: air-cladding and cladding-core
%   Exact solution via integrating over u

    % Sunlight outside this value doesn't enter the core
    % If integrating from 0 to 1, dclad would have complex values
    uCut = min(1, q*nclad);

    function result = integrand(u)
        % Distance traveled inside core and cladding, respectively
        dcore = @(u) diameter*sqrt(q^2 - u.^2/ncore^2);
        dclad = @(u) diameter*(sqrt(1 - u.^2/nclad^2) - sqrt(q^2 - u.^2/nclad^2));

        cosAir = @(u) sqrt(1 - u.^2);
        cosClad1 = @(u) sqrt(1 - u.^2/nclad^2);
        cosClad2 = @(u) sqrt(1 - u.^2/(q*nclad)^2);
        cosCore = @(u) sqrt(1 - u.^2/(q*ncore)^2);

        % Fresnel reflection coefficients in each interface
        R_F1 = @(u) (((cosAir(u) - nclad*cosClad1(u))./(cosAir(u) + nclad*cosClad1(u)))^2 + ((nclad*cosAir(u) - cosClad1(u))./(nclad*cosAir(u) + cosClad1(u)))^2)/2;
        R_F2 = @(u) (((nclad*cosClad2(u) - ncore*cosCore(u))./(nclad*cosClad2(u) + ncore*cosCore(u)))^2 + ((ncore*cosClad2(u) - nclad*cosCore(u))./(ncore*cosClad2(u) + nclad*cosCore(u)))^2)/2;

        expCore = @(u) exp(-alfaCore*dcore(u));
        expClad = @(u) exp(-alfaClad*dclad(u));

        % Result is obtained analytically by solving a system of linear
        % equations
        % The corresponding matrix and vector are defined as functions
        v = @(u) [0 ; 0 ; (1 - expCore(u))*alfaDopant/alfaCore];
        A = @(u) [1 -R_F2(u).*expClad(u) -(1-R_F2(u)).*expClad(u) ; -R_F1(u).*expClad(u) 1 0 ; 0 -(1-R_F2(u)).*expCore(u) 1-R_F2(u).*expCore(u)];

        result = zeros(1, length(u));

        % To use matrices and vectors, you need to iterate over the values
        % of u
        % This is also why you need a named function like this
        for i = 1:length(u)
            % Solve the system, get a vector
            solvedVector = A(u(i))\v(u(i));
            
            % Solution is the vector's first element, plus transmission
            % coefficient in the air-cladding interface
            result(i) = solvedVector(1)*(1-R_F1(u(i)));
        end
    end

    % Integrate numerically
    absorption = quadgk(@integrand, 0, uCut);

end