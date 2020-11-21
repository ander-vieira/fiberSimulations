function sigmaemi = sigmaemi_C1(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  sigmaabs_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [500e-9 750e-9] nm.')
    ll = linspace(350e-9,700e-9,1000);
    %figure
    plot(ll,sigmaemi_C1(ll))
    xlabel('\lambda (m)')
    ylabel('\sigma (m^2)')
    %axis([350e-9 750e-9 0e-20 1.6e-20])
    legend('a','b','c','d')
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    sigmaemi = zeros(size(lambda));
    for i = 1:length(lambda)
        sigmaemi(i) = sigmaemi_C1(lambda(i));
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
    340        0
    355        0
    356        0
       360   4.9145e-23
       366.41   2.6489e-22
       372.81   1.0917e-21
       379.22   3.1423e-21
       385.63   6.1143e-21
       392.03   8.6605e-21
       398.44   9.9482e-21
       404.84   1.0166e-20
       411.25   9.8653e-21
       417.66   9.0304e-21
       424.06   7.8652e-21
       430.47   6.6329e-21
       436.88   5.4467e-21
       443.28   4.3509e-21
       449.69   3.3446e-21
        456.1   2.4781e-21
        462.5   1.8202e-21
       468.91   1.3567e-21
       475.32   9.9307e-22
       481.72   7.1811e-22
       488.13   5.2067e-22
       494.53   3.8692e-22
       500.94   2.9519e-22
       507.35   2.2926e-22
       513.75   1.8869e-22
       520.16   1.6171e-22
       526.57   1.4069e-22
       532.97   1.0012e-22
       539.38   9.1501e-23
       545.79    8.333e-23
       552.19   7.5599e-23
        558.6     6.83e-23
       565.01   6.1426e-23
       571.41    5.497e-23
       577.82   4.8922e-23
       584.22   4.3275e-23
       590.63   3.8019e-23
       597.04   3.3145e-23
       603.44   2.8644e-23
       609.85   2.4505e-23
       616.26   2.0718e-23
       622.66   1.7274e-23
       629.07   1.4163e-23
       635.48   1.1373e-23
       641.88    8.894e-24
       648.29   6.7159e-24
       654.69   4.8279e-24
        661.1   3.2192e-24
       667.51   1.8794e-24
       673.91   7.9787e-25  
    700         0
    751               0
    752               0
    753               0
    754               0
    755               0
    ];
%Poner los datos en valores absolutos:
%sigmapico = 1.5e-20;      % Valor al que hay que normalizar
datos(:,2) = datos(:,2)*0.823464;

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

% sigmaemi = sigmaemi * 8.1564; % Fuchtbauer-Ladenburg correction

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
