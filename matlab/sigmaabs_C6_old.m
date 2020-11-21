function sigmaemi = sigmaabs_C6(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  sigmaabs_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [500e-9 750e-9] nm.')
    ll = linspace(250e-9,700e-9,1000);
    %figure
    plot(ll,sigmaabs_C6(ll))
    xlabel('\lambda (m)')
    ylabel('\sigma (m^2)')
    %axis([350e-9 750e-9 0e-20 1.6e-20])
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    sigmaemi = zeros(size(lambda));
    for i = 1:length(lambda)
        sigmaemi(i) = sigmaabs_C6(lambda(i));
    end
    return
end

% SI LAMBDA ES ESCALAR:

% Comprobación de rango:
if lambda>750e-9, sigmaemi = 0; return; end
% if lambda<249e-9, error('Seguramente te has equivocado de lambda. Demasiado pequeña.'); end
if lambda<249e-9, sigmaemi = 0; return; end

% Datos brutos de Igor (recortados en lambda según rango razonable):
lambda = lambda*1e9;    % trabajar en nm para mejorar el condicionamiento
datos = [               % primera columna lambda (nm), segunda columna sigma emisión sin normalizar (a.u.)    
    230        0
    235        0
    237.18  0
    249.38  0
  221.67            0
       380.23            0
          384   1.3381e-22
       387.04   3.5184e-22
       390.07   7.5779e-22
       396.88   3.0041e-21
       401.04   4.8175e-21
       406.72   6.6849e-21
       411.26   8.2547e-21
          418   1.0068e-20
          427   1.2044e-20
       435.86   1.3722e-20
       439.26   1.4263e-20
       441.53   1.4452e-20
       444.94   1.4561e-20
        449.1   1.4588e-20
       452.13   1.4506e-20
       456.67   1.4452e-20
        459.7   1.4506e-20
       462.72   1.4669e-20
          465   1.4777e-20
       467.27   1.4777e-20
       469.16   1.4723e-20
       471.81   1.4425e-20
       475.21   1.3316e-20
       479.75   1.0744e-20
       484.67   7.3073e-21
       488.46   4.4385e-21
       491.86   2.1651e-21
       495.65   1.0014e-21
        498.3   5.9542e-22
       499.81   4.0597e-22
       502.08   2.1651e-22
          506   1.0035e-22
       516.46            0
       699.62            0
    752               0
    753               0
    754               0
    755               0
    ];
%Poner los datos en valores absolutos:
%sigmapico = 1.5e-20;      % Valor al que hay que normalizar
datos(:,2) = datos(:,2)*0.96113;

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
