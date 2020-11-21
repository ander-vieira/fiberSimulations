function sigmaabs = sigmaabs_LumogenRed(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  sigmaabs_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [350e-9 750e-9] nm.')
    ll = linspace(350e-9,750e-9,1000);
    %figure
    plot(ll,sigmaabs_LumogenRed(ll),'m')
    xlabel('\lambda (m)')
    ylabel('\sigma^a (m^2)')
    %axis([350e-9 750e-9 0e-20 1.6e-20])
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    sigmaabs = zeros(size(lambda));
    for i = 1:length(lambda)
        sigmaabs(i) = sigmaabs_LumogenRed(lambda(i));
    end
    return
end

% SI LAMBDA ES ESCALAR:

% Comprobación de rango:
if lambda>750e-9, sigmaabs = 0; return; end
% if lambda<301e-9, error('Seguramente te has equivocado de lambda. Demasiado pequeña.'); end
if lambda<345e-9, sigmaabs = 0; return; end

% Datos brutos de Igor (recortados en lambda según rango razonable):
lambda = lambda*1e9;    % trabajar en nm para mejorar el condicionamiento
datos = [               % primera columna lambda (nm), segunda columna sigma emisión sin normalizar (a.u.)
    345               0
    346               0
    347               0
    348               0
    349               0
    350               0
    365.7             0
    367.48     0.055971
    368.91      0.10615
    369.63      0.11001
    370.7       0.10229
    372.48     0.075271
    374.62     0.055971
    375.33     0.055971
    378.9      0.079131
    384.97       0.1158
    387.47      0.13124
    392.82      0.16212
    398.53      0.20265
    403.17      0.23353
    407.09      0.27214
    411.02      0.30881
    414.23      0.34548
    417.8       0.38022
    421.01      0.41496
    424.93      0.45163
    429.57      0.49216
    435.64      0.52304
    441.7       0.54234
    444.92      0.54234
    448.84      0.52111
    452.77      0.48251
    456.33      0.44198
    458.12      0.40724
    460.62      0.37443
    463.11      0.33583
    465.61      0.30495
    469.89      0.26828
    472.75      0.25476
    474.53      0.25669
    477.39      0.27021
    483.1       0.29916
    487.38       0.3339
    491.3       0.36864
    495.58      0.40531
    498.8       0.44005
    502.72      0.47865
    506.65      0.53076
    510.93      0.60024
    514.14      0.67551
    517.71      0.74885
    520.56      0.82027
    523.77      0.88975
    526.27      0.94572
    528.77      0.98239
    532.34        1.021
    535.55       1.0306
    539.83       1.0326
    543.4        1.0364
    546.25       1.0499
    550.18       1.0828
    551.96       1.1194
    556.24       1.2082
    559.1        1.2835
    562.31       1.3568
    565.17       1.4398
    566.95       1.4784
    568.73       1.5151
    570.87       1.5363
    573.37       1.5382
    574.8        1.5247
    576.23       1.5093
    578.37       1.4707
    579.44        1.434
    581.22         1.38
    583.36       1.3086
    586.22       1.1677
    589.07        1.021
    591.57      0.87238
    593.71      0.73148
    596.57      0.58866
    599.78      0.44391
    603.35      0.30109
    608.34      0.16019
    611.91      0.10229
    615.12     0.067551
    620.12     0.030881
    627.61    0.0057901
    629.39    0
    750       0
    751       0
    752       0
    753       0
    754       0
    755       0    ];

%Poner los datos en valores absolutos:
%sigmapico = 1.5e-20; %Para Rh6G. Valor al que hay que normalizar
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

% if sigmaabs <0 || iizq<=4
%     sigmaabs = (sigmas(iizq)+sigmas(ider))/2;
% end

% Esto es para dibujar en rojo la parábola de interpolación utilizada junto con los puntos dato utilizados y el
% valor devuelto:
% figure
% hold on
% plot(lambdas,sigmas,'.')
% plot(lambda,sigmaabs,'or')
% xx = linspace(xk(1),xk(end));
% yy = a*xx.^2 + b*xx + c;
% plot(xx,yy,'r')

return
