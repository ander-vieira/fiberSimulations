function sigmaabs = sigmaabs_Rh6G(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  sigmaabs_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [410e-9 750e-9] nm.')
    ll = linspace(418e-9,750e-9,1000);
    %figure
    plot(ll,sigmaabs_Rh6G(ll))
    xlabel('\lambda (m)')
    ylabel('\sigma (m^2)')
    %axis([400e-9 700e-9 0e-20 5e-20])
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    sigmaabs = zeros(size(lambda));
    for i = 1:length(lambda)
        sigmaabs(i) = sigmaabs_Rh6G(lambda(i));
    end
    return
end

% SI LAMBDA ES ESCALAR:

% Comprobación de rango:
if lambda>700e-9, sigmaabs = 0; return; end
% if lambda<301e-9, error('Seguramente te has equivocado de lambda. Demasiado pequeña.'); end
if lambda<301e-9, sigmaabs = 0; return; end

% Datos brutos de Igor (recortados en lambda según rango razonable):
lambda = lambda*1e9;    % trabajar en nm para mejorar el condicionamiento
datos = [               % primera columna lambda (nm), segunda columna sigma emisión sin normalizar (a.u.)
    300                                      0
    414                                      0
    415                                      0
    416                                      0
    417                                      0
    418.348623853211                         0
    428.440366972477        0.0127298444130127
    437.920489296636        0.0325318246110325
    449.54128440367         0.0707213578500708
    457.798165137615         0.101838755304102
    466.97247706422          0.145685997171146
    477.064220183486         0.209335219236209
    483.180428134557         0.271570014144272
    488.990825688073         0.339462517680339
    493.88379204893          0.391796322489392
    498.470948012232          0.43988684582744
    504.281345565749         0.499292786421499
    507.951070336392         0.548797736916549
    512.844036697248         0.640735502121641
    516.51376146789          0.731258840169731
    519.877675840979         0.823196605374823
    523.547400611621          0.90947666195191
    525.382262996942         0.951909476661952
    527.217125382263         0.975954738330976
    527.82874617737          0.987270155586987
    529.663608562691         0.998585572842999
    530.581039755352                         1
    532.415902140673         0.995756718528996
    533.639143730887          0.98019801980198
    535.474006116208         0.951909476661952
    537.920489296636         0.892503536067893
    541.590214067278         0.770862800565771
    545.259938837921         0.636492220650637
    548.929663608563         0.497878359264498
    555.351681957187         0.311173974540311
    561.162079510704         0.203677510608204
    564.220183486239         0.164073550212164
    567.889908256881          0.13012729844413
    572.782874617737         0.104667609618105
    578.593272171254        0.0834512022630836
    588.073394495413        0.0608203677510609
    603.363914373089        0.0410183875530411
    620.489296636086        0.0254596888260255
    646.788990825688        0.0113154172560113
    670.642201834862       0.00282885431400299
    700                                      0
    701                                      0
    702                                      0
    703                                      0
    704                                      0];
       
%Poner los datos en valores absolutos:
sigmapico = 4.3298e-20; %Para Rh6G. Valor al que hay que normalizar
datos(:,2) = datos(:,2)*sigmapico;


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
if iizq>4
% Para suavizar cogemos los 4 puntos a cada lado sin tocar, pero el punto nº 5 se coge interpolado:
% si p.ej. lambda está a un 30% de distancia del dato bruto más cercano a su izquierda y a un 70% de distancia
% del dato bruto más cercano a su derecha, se pondera el quinto punto entre el cuarto y el quinto a cada lado:
% Fracción de lejanía al punto izquierdo fli:
fli = (lambda-lambdas(iizq))/(lambdas(ider)-lambdas(iizq));
fld = 1-fli;    % Fracción de lejanía a la derecha.
% Punto ficticio izquierdo:
x1 = fld*lambdas(iizq-3) + fli*lambdas(iizq-2);
y1 = fld*sigmas(iizq-3)  + fli*sigmas(iizq-2);
% Punto ficticio derecho:
xend = fld*lambdas(ider+2) + fli*lambdas(ider+3);
yend = fld*sigmas(ider+2)  + fli*sigmas(ider+3);
% Montar xk,yk y aprovechar el código de la parábola mínimo-cuadrática:
xk = [x1; lambdas(iizq-1:ider+1); xend]; % abscisas
yk = [y1; sigmas(iizq-1:ider+1); yend]; % ordenadas
% Aproximamos mediante parábola de grado 2 mínimo-cuadrática:
A  = [xk ones(size(xk))];
b  = yk;
abc = A\b;
a = abc(1);
b = abc(2);
%c = abc(3);
sigmaabs = a*lambda + b;
end

if sigmaabs <0 || iizq<=4
    sigmaabs = (sigmas(iizq)+sigmas(ider))/2;
end

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
