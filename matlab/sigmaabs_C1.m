function sigmaabs = sigmaabs_C1(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  sigmaabs_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [500e-9 750e-9] nm.')
    ll = linspace(250e-9,700e-9,1000);
    %figure
    plot(ll,sigmaabs_C1(ll))
    xlabel('\lambda (m)')
    ylabel('\sigma^e (m^2)')
    %axis([350e-9 750e-9 0e-20 1.6e-20])
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    sigmaabs = zeros(size(lambda));
    for i = 1:length(lambda)
        sigmaabs(i) = sigmaabs_C1(lambda(i));
    end
    return
end

% SI LAMBDA ES ESCALAR:

% Comprobación de rango:
if lambda>750e-9, sigmaabs = 0; return; end
% if lambda<249e-9, error('Seguramente te has equivocado de lambda. Demasiado pequeña.'); end
if lambda<249e-9, sigmaabs = 0; return; end

% Datos brutos de Igor (recortados en lambda según rango razonable):
lambda = lambda*1e9;    % trabajar en nm para mejorar el condicionamiento
datos = [               % primera columna lambda (nm), segunda columna sigma emisión sin normalizar (a.u.)    
    230        0
    235        0
    237.18  0
    249.38  0
219.49            0
       282.53            0
       285.19   3.4839e-22
       291.27    1.012e-21
       299.62   2.1235e-21
       309.11   3.5503e-21
        315.3   4.6784e-21
       320.51   5.9725e-21
       325.82   7.6978e-21
       328.86   8.7264e-21
       331.52   9.2905e-21
        333.8   9.6887e-21
       336.46   9.9873e-21
       337.97   1.0087e-20
       341.01    1.017e-20
       344.43    1.002e-20
       347.09   9.7716e-21
        351.5   9.2905e-21
        358.1   8.4775e-21
       364.94   7.3494e-21
       370.63   5.9558e-21
       374.81   4.5457e-21
       381.27   2.5382e-21
       385.44   1.5097e-21
       389.24   8.6268e-22
        393.8   3.6498e-22
        396.5   1.4931e-22
       399.49            0
          600            0
    752               0
    753               0
    754               0
    755               0
    ];
%Poner los datos en valores absolutos:
%sigmapico = 1.5e-20;      % Valor al que hay que normalizar
datos(:,2) = datos(:,2)*0.8207;

% Recojo los datos "reales":
lambdas = datos(:,1); % en nm
sigmas  = datos(:,2); % normalizadas (a.u.)

% Empieza la aproximación.
% Encontramos el índice del dato bruto inmediatamente a la izquierda de la lambda de entrada:
menores = find(lambdas<lambda);
mayores = find(lambdas>=lambda);
iizq  = menores(end);
ider  = mayores(1);
if ider-iizq ~= 1
    error('Revisa el código de sigmaabs_Rh6G')
end

sigmaabs=0;
number=2;
if iizq>number
    % Para suavizar cogemos number puntos a cada lado sin tocar, pero el punto nº 5 se coge interpolado:
    % si p.ej. lambda está a un 30% de distancia del dato bruto más cercano a su izquierda y a un 70% de distancia
    % del dato bruto más cercano a su derecha, se pondera el quinto punto entre el cuarto y el quinto a cada lado:
    % Fracción de lejanía al punto izquierdo fli:
    fli = (lambda-lambdas(iizq))/(lambdas(ider)-lambdas(iizq));
    fld = 1-fli;    % Fracción de lejanía a la derecha.
    % Punto ficticio izquierdo:
    x1 = fld*lambdas(iizq-(number-1)) + fli*lambdas(iizq-(number-2));
    y1 = fld*sigmas(iizq-(number-1))  + fli*sigmas(iizq-(number-2));
    % Punto ficticio derecho:
    xend = fld*lambdas(ider+(number-2)) + fli*lambdas(ider+(number-1));
    yend = fld*sigmas(ider+(number-2))  + fli*sigmas(ider+(number-1));
    % Montar xk,yk y aprovechar el código de la parábola mínimo-cuadrática:
    xk = [x1; lambdas(iizq-(number-2):iizq+(number-1)); xend]; % abscisas
    yk = [y1; sigmas(iizq-(number-2):iizq+(number-1)); yend]; % ordenadas
    % Aproximamos mediante parábola de grado 2 mínimo-cuadrática:
    A  = [xk ones(size(xk))];
    b  = yk;
    abc = A\b;
    a = abc(1);
    b = abc(2);
    %c = abc(3);
    sigmaabs = a*lambda + b;
end
% if sigmaemi <0 || iizq<=number
%     sigmaemi = (sigmas(iizq)+sigmas(ider))/2;
% end


% Esto es para dibujar en rojo la parábola de interpolación utilizada junto con los puntos dato utilizados y el
% valor devuelto:
% figure
% hold on
% plot(lambdas,sigmas,'.')
% plot(lambda,sigmaemi,'or')
% xx = linspace(xk(1),xk(end));
% yy = a*xx.^2 + b*xx + c;
% plot(xx,yy,'r')

return
