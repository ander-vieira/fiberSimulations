function alfa = valuesalfaPMMA(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  alfa_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [350e-9 750e-9] nm.')
    ll = linspace(400e-9,800e-9,1000);
    %figure
    plot(ll,valuesalfaPMMA(ll),'m')
    xlabel('\lambda (m)')
    ylabel('\alpha (m^{-1})')
    %axis([350e-9 750e-9 0e-20 1.6e-20])
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    alfa = zeros(size(lambda));
    for i = 1:length(lambda)
        alfa(i) = valuesalfaPMMA(lambda(i));
    end
    return
end

% SI LAMBDA ES ESCALAR:

% Comprobación de rango:
if lambda>800e-9, alfa = 0; return; end
if lambda<400e-9, alfa = 0; return; end

% Datos brutos de Igor (recortados en lambda según rango razonable):
lambda = lambda*1e9;    % trabajar en nm para mejorar el condicionamiento
datos = [               % primera columna lambda (nm), segunda columna alfa. % Nota:
% Las abscisas van en nm lineales, pero las ordenadas van en escala logarítmica.
% El cero en la segunda columna de `datos` corresponde a 0.01 dB/m,
% y el 10 corresponde a 10 dB/m.
%
% Para pasar ordenadas (2ª columna de datosbrutos) a valores absolutos de alfa:
%
% alfa = 0.01 * 10^(0.3 * h))
    396          4.434
    397          4.434 
    398          4.344
    399          4.344
    400.00       4.3430
    414.52       3.9696
    427.68       3.6791
    440.38       3.4302
    455.35       3.1674
    475.32       2.9322
    488.93       2.8216
    497.55       2.8077
    507.08       2.8216
    518.42       2.9184
    526.13       3.0014
    532.49       3.1120
    539.29       3.2780
    543.83       3.3748
    547.46       3.3887
    552.00       3.3748
    557.44       3.2227
    560.62       3.1120
    565.15       2.9876
    570.60       2.8492
    576.50       2.7801
    579.67       2.7524
    584.21       2.8077
    588.29       2.9046
    592.38       3.1259
    597.82       3.5408
    603.72       4.2462
    609.62       4.9101
    613.70       5.2144
    618.24       5.4495
    621.42       5.5048
    623.23       5.5187
    627.31       5.4772
    630.04       5.3804
    633.21       5.2144
    637.30       4.8824
    642.29       4.4675
    645.46       4.2185
    648.19       4.0664
    650.00       4.0111
    653.18       3.9696
    656.35       4.0111
    661.34       4.1909
    666.33       4.4675
    669.51       4.6473
    670.87       4.7026
    674.50       4.7718
    679.04       4.8271
    683.58       4.8963
    687.21       5.0069
    692.20       5.2835
    697.19       5.7538
    702.18       6.5145
    708.98       7.3306
    713.52       7.7178
    718.97       7.9806
    721.69       8.0498
    723.96       8.0636
    727.13       8.0360
    730.31       7.9530
    735.30       7.6349
    740.74       7.2199
    746.64       6.8050
    751.63       6.5145
    756.17       6.3485
    758.44       6.2794
    759.80       6.2517
    762.07       6.2241
    765.70       6.2932
    769.78       6.4730
    776.13       6.9018
    781.13       7.2199
    786.12       7.4274
    792.92       7.5795
    798.37       7.6487
    809.26       7.6349
    814.25       7.7040
    819.24       7.8008
    826.50       8.0775
    838.75       8.6030
    849.64       9.0456
    860.07       9.4191
    869.60       9.6819
    876.41       9.8617
    884.12      10.0000];

%Poner los datos en valores absolutos:
datos(:,2) = 0.01 * 10.^(0.3 * datos(:,2)); %puesto en dB/m
datos(:,2) = datos(:,2)*0.5; %fibra con un mínimo de 0.08 dB/m el el verde
%quito datos(:,2) = datos(:,2)*0.25; %fibra con un mínimo de 0.04 dB/m el el verde
%quito datos(:,2) = datos(:,2)*1e-6; %fibra sin atenuación

% Recojo los datos "reales":
lambdas = datos(:,1); % en nm
alfas  = datos(:,2)*log(10); % pasado a Np/m

% Empieza la aproximación.
% Encontramos el índice del dato bruto inmediatamente a la izquierda de la lambda de entrada:
menores = find(lambdas<lambda);
mayores = find(lambdas>=lambda);
iizq  = menores(end);
ider  = mayores(1);
if ider-iizq ~= 1
    error('Revisa el código de alfa_Rh6G')
end

alfa=0;
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
    y1 = fld*alfas(iizq-(number-1))  + fli*alfas(iizq-(number-2));
    % Punto ficticio derecho:
    xend = fld*lambdas(ider+(number-2)) + fli*lambdas(ider+(number-1));
    yend = fld*alfas(ider+(number-2))  + fli*alfas(ider+(number-1));
    % Montar xk,yk y aprovechar el código de la parábola mínimo-cuadrática:
    xk = [x1; lambdas(iizq-(number-2):iizq+(number-1)); xend]; % abscisas
    yk = [y1; alfas(iizq-(number-2):iizq+(number-1)); yend]; % ordenadas
    % Aproximamos mediante parábola de grado 2 mínimo-cuadrática:
    A  = [xk ones(size(xk))];
    b  = yk;
    abc = A\b;
    a = abc(1);
    b = abc(2);
    %c = abc(3);
    alfa = a*lambda + b;
end

% if alfa <0 || iizq<=4
%     alfa = (alfas(iizq)+alfas(ider))/2;
% end

% Esto es para dibujar en rojo la parábola de interpolación utilizada junto con los puntos dato utilizados y el
% valor devuelto:
% figure
% hold on
% plot(lambdas,alfas,'.')
% plot(lambda,alfa,'or')
% xx = linspace(xk(1),xk(end));
% yy = a*xx.^2 + b*xx + c;
% plot(xx,yy,'r')

return
