function sigmaemi = sigmaemi_LumogenRed(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  sigmaabs_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [500e-9 750e-9] nm.')
    ll = linspace(500e-9,750e-9,1000);
    %figure
    plot(ll,sigmaemi_LumogenRed(ll))
    xlabel('\lambda (m)')
    ylabel('\sigma^e (m^2)')
    %axis([350e-9 750e-9 0e-20 1.6e-20])
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    sigmaemi = zeros(size(lambda));
    for i = 1:length(lambda)
        sigmaemi(i) = sigmaemi_LumogenRed(lambda(i));
    end
    return
end

% SI LAMBDA ES ESCALAR:

% Comprobación de rango:
if lambda>750e-9, sigmaemi = 0; return; end
% if lambda<301e-9, error('Seguramente te has equivocado de lambda. Demasiado pequeña.'); end
if lambda<345e-9, sigmaemi = 0; return; end

% Datos brutos de Igor (recortados en lambda según rango razonable):
lambda = lambda*1e9;    % trabajar en nm para mejorar el condicionamiento
datos = [               % primera columna lambda (nm), segunda columna sigma emisión sin normalizar (a.u.)
    345               0
    346               0
    347               0
    348               0
    349               0
    350               0
    500          0.0005
    515          0.0005
    516          0.0005
    520           0.013
    524.84     0.019417
    532.69     0.025242
    538.76      0.03495
    544.83     0.054367
    552.68      0.10291
    557.67      0.15339
    561.6       0.21359
    564.09      0.27184
    566.95      0.32038
    572.66      0.48348
    577.65      0.68153
    582.65      0.90483
    585.86       1.0757
    589.43       1.2446
    592.64       1.3689
    595.5        1.4582
    597.99       1.5165
    600.85        1.532
    602.99       1.5165
    606.56       1.4563
    610.12       1.3514
    614.41       1.1961
    619.4        1.0194
    624.04      0.87959
    631.18      0.72425
    636.53      0.63688
    645.81      0.52232
    655.08      0.41746
    662.58       0.3262
    670.07      0.23106
    679.35      0.15922
    687.91      0.11456
    697.55     0.081551
    707.18     0.062134
    718.96     0.044659
    728.95      0.03495
    736.08     0.031067
    738.94       0.0233
    750          0.0233
    751               0
    752               0
    753               0
    754               0
    755               0
    ];
%Poner los datos en valores absolutos:
%sigmapico = 1.5e-20;      % Valor al que hay que normalizar
datos(:,2) = datos(:,2)*1e-20;

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

sigmaemi=0;
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
    sigmaemi = a*lambda + b;
end
% if sigmaemi <0 || iizq<=number
%     sigmaemi = (sigmas(iizq)+sigmas(ider))/2;
% end

sigmaemi = sigmaemi * 1.4215; % Fuchtbauer-Ladenburg correction

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
