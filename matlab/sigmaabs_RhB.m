function sigmaabs = sigmaabs_RhB(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  sigmaabs_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [350e-9 750e-9] nm.')
    ll = linspace(350e-9,750e-9,1000);
    %figure
    plot(ll,sigmaabs_RhB(ll))
    xlabel('\lambda (m)')
    ylabel('\sigma^a (m^2)')
    %axis([350e-9 750e-9 0e-20 1.6e-20])
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    sigmaabs = zeros(size(lambda));
    for i = 1:length(lambda)
        sigmaabs(i) = sigmaabs_RhB(lambda(i));
    end
    return
end

% SI LAMBDA ES ESCALAR:

% Comprobación de rango:
if lambda>700e-9, sigmaabs = 0; return; end
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
   400                     0.0427135678391962
    412.477876106195        0.0389447236180906
    427.345132743363        0.0326633165829147
    439.823008849558        0.0238693467336685
    446.725663716814        0.0226130653266334
    453.628318584071        0.0226130653266334
    462.12389380531         0.0288944723618092
    470.353982300885        0.0402010050251258
    482.300884955752        0.0690954773869349
    490.265486725664        0.0942211055276384
    495.575221238938         0.115577889447236
    503.008849557522         0.165829145728643
    508.053097345133         0.218592964824121
    514.159292035398         0.285175879396985
    518.672566371681         0.324120603015076
    524.24778761062          0.360552763819096
    529.026548672566         0.385678391959799
    532.743362831858         0.417085427135679
    536.991150442478         0.473618090452261
    543.628318584071         0.634422110552764
    549.203539823009         0.809045226130653
    552.920353982301         0.918341708542714
    555.575221238938         0.971105527638191
    557.16814159292          0.993718592964824
    559.29203539823                          1
    561.681415929204         0.992462311557789
    563.53982300885          0.956030150753769
    567.787610619469          0.82286432160804
    570.973451327434         0.685929648241206
    575.221238938053         0.502512562814071
    579.734513274336         0.340452261306533
    586.902654867257         0.160804020100503
    592.743362831858        0.0841708542713569
    597.522123893805        0.0515075376884424
    603.893805309735        0.0263819095477388
    612.920353982301        0.0113065326633168
    624.601769911505       0.00502512562814079
    641.061946902655       0.00251256281407049
    750               0
    751               0
    752               0
    753               0
    754               0
    755               0
    ];

%Poner los datos en valores absolutos:
%sigmapico = 1.5e-20; %Para Rh6G. Valor al que hay que normalizar
datos(:,2) = datos(:,2)*3.37e-20;


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
