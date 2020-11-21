function sigmaemi = sigmaemi_Rh6G(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  sigmaabs_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [500e-9 750e-9] nm.')
    ll = linspace(500e-9,750e-9,1000);
    %figure
    plot(ll,sigmaemi_Rh6G(ll))
    xlabel('\lambda (m)')
    ylabel('\sigma (m^2)')
    %axis([400e-9 700e-9 0e-20 5e-20])
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    sigmaemi = zeros(size(lambda));
    for i = 1:length(lambda)
        sigmaemi(i) = sigmaemi_Rh6G(lambda(i));
    end
    return
end

% SI LAMBDA ES ESCALAR:

% Comprobación de rango:
if lambda>700e-9, sigmaemi = 0; return; end
% if lambda<301e-9, error('Seguramente te has equivocado de lambda. Demasiado pequeña.'); end
if lambda<301e-9, sigmaemi = 0; return; end

% Datos brutos de Igor (recortados en lambda según rango razonable):
lambda = lambda*1e9;    % trabajar en nm para mejorar el condicionamiento
datos = [               % primera columna lambda (nm), segunda columna sigma emisión sin normalizar (a.u.)
    300                                      0
    498                                      0
    499                                      0
    500                                      0
    501                                      0
    502.568644818423        0.0195121951219512
    504.694419840567        0.0463414634146341
    508.945969884854         0.123170731707317
    513.197519929141         0.228048780487805
    518.5119574845           0.360975609756097
    523.56067316209          0.510975609756097
    528.875110717449         0.670731707317073
    532.5952170062           0.792682926829268
    535.252435783879         0.870731707317073
    538.972542072631         0.953658536585365
    541.62976085031          0.986585365853658
    542.161204605846                     0.993
    544.021257750221                         1
    546.412754650133                     0.993
    548.272807794508         0.969512195121951
    551.727192205492         0.908536585365853
    557.04162976085          0.785365853658536
    562.621789193977         0.628048780487805
    568.467670504872         0.485365853658536
    574.579273693534         0.373170731707317
    586.536758193091         0.230487804878049
    597.962798937112         0.142682926829268
    608.325952170062        0.0853658536585365
    622.674933569531        0.0390243902439024
    637.555358724535        0.0158536585365853
    655.358724534987       0.00365853658536577
    678.210806023029       0.00121951219512187
    700                                      0
    701                                      0
    702                                      0
    703                                      0
    704                                      0];
    
%Poner los datos en valores absolutos:
sigmapico = 1.6553e-20;      % Valor al que hay que normalizar
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

sigmaemi=0;
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
xk = [x1; lambdas(iizq-2:iizq+3); xend]; % abscisas
yk = [y1; sigmas(iizq-2:iizq+3); yend]; % ordenadas
% Aproximamos mediante parábola de grado 2 mínimo-cuadrática:
A  = [xk ones(size(xk))];
b  = yk;
abc = A\b;
a = abc(1);
b = abc(2);
%c = abc(3);
sigmaemi = a*lambda + b;
end
if sigmaemi <0 || iizq<=4
    sigmaemi = (sigmas(iizq)+sigmas(ider))/2;
end

sigmaemi = sigmaemi*1.230876; % Fuchtbauer-Ladenburg correction

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
