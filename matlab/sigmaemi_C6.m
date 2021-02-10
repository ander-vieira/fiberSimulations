function sigmaemi = sigmaemi_C6(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  sigmaabs_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [500e-9 750e-9] nm.')
    ll = linspace(400e-9,700e-9,1000);
    %figure
    plot(ll,sigmaemi_C6(ll))
    xlabel('\lambda (m)')
    ylabel('\sigma^e (m^2)')
    %axis([350e-9 750e-9 0e-20 1.6e-20])
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    sigmaemi = zeros(size(lambda));
    for i = 1:length(lambda)
        sigmaemi(i) = sigmaemi_C6(lambda(i));
    end
    return
end

% SI LAMBDA ES ESCALAR:

% Comprobación de rango:
if lambda>750e-9, sigmaemi = 0; return; end
% if lambda<301e-9, error('Seguramente te has equivocado de lambda. Demasiado pequeña.'); end
if lambda<301e-9, sigmaemi = 0; return; end

% Datos brutos de Igor (recortados en lambda según rango razonable):
lambda = lambda*1e9;    % trabajar en nm para mejorar el condicionamiento
datos = [               % primera columna lambda (nm), segunda columna sigma emisión sin normalizar (a.u.)
    300        0
    350        0
    375        0
    400        0
    464   3.9298e-23
       468.72   6.6009e-23
       473.45   1.8698e-22
       478.17    6.018e-22
        482.9   2.2119e-21
       487.62   5.8522e-21
       492.35   1.0588e-20
       497.07   1.4141e-20
        501.8   1.4717e-20
       506.52   1.3893e-20
       511.25   1.2745e-20
       515.97   1.1911e-20
        520.7   1.1242e-20
       525.42   1.0577e-20
       530.15    9.621e-21
       534.87   8.5065e-21
        539.6   7.2891e-21
       544.32   6.1409e-21
       549.05   5.1048e-21
       553.77   4.2764e-21
       558.49   3.6473e-21
       563.22   3.1214e-21
       567.94   2.6517e-21
       572.67     2.28e-21
       577.39   2.0143e-21
       582.12    1.838e-21
       586.84   1.7054e-21
       591.57   1.5746e-21
       596.29   1.2855e-21
       601.02   8.7548e-22
       605.74   6.1496e-22
       610.47   5.3208e-22
       615.19   3.7295e-22
       619.92   3.3256e-22
       624.64   2.9503e-22
       629.37   2.6027e-22
       634.09    2.282e-22
       638.81    1.987e-22
       643.54   1.7168e-22
       648.26   1.4703e-22
       652.99   1.2462e-22
       657.71   1.0435e-22
       662.44   8.6097e-23
       667.16   6.9753e-23
       671.89   5.5206e-23
       676.61   4.2344e-23
       681.34   3.1061e-23
       686.06   2.1255e-23
       690.79   1.2827e-23
       695.51   5.6822e-24
        700         0
    751               0
    760               0
    765               0
    770              0
    775               0
    ];
%Poner los datos en valores absolutos:
%sigmapico = 1.5e-20;      % Valor al que hay que normalizar
datos(:,2) = datos(:,2)*0.996793;

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

sigmaemi = sigmaemi * 1.5624; % Fuchtbauer-Ladenburg correction

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
