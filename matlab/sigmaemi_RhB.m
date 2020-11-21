function sigmaemi = sigmaemi_RhB(lambda)
% Sigma de absorción.
% Para ver la curva a golpe de vista, ejecutar  sigmaabs_Rh6G  sin parámetros.
% Medidas nuestras del laboratorio, febrero de 2013.

% Jon Arrue, Felipe Jiménez, Igor Ayesta.


% Dominio de definición:
if nargin==0
    %disp('Rango de lambdas en que esta función funciona = [500e-9 750e-9] nm.')
    ll = linspace(500e-9,750e-9,1000);
    %figure
    plot(ll,sigmaemi_RhB(ll))
    xlabel('\lambda (m)')
    ylabel('\sigma^e (m^2)')
    %axis([350e-9 750e-9 0e-20 1.6e-20])
    return
end

% Si lambda es vector, autollamarse:
if length(lambda)>1
    sigmaemi = zeros(size(lambda));
    for i = 1:length(lambda)
        sigmaemi(i) = sigmaemi_RhB(lambda(i));
    end
    return
end

% SI LAMBDA ES ESCALAR:

% Comprobación de rango:
if lambda>750e-9, sigmaemi = 0; return; end
% if lambda<301e-9, error('Seguramente te has equivocado de lambda. Demasiado pequeña.'); end
if lambda<350e-9, sigmaemi = 0; return; end

% Datos brutos de Igor (recortados en lambda según rango razonable):
lambda = lambda*1e9;    % trabajar en nm para mejorar el condicionamiento
datos = [               % primera columna lambda (nm), segunda columna sigma emisión sin normalizar (a.u.)
    350               0
     500.53050397878         0.0100502512562819
    508.488063660477        0.0201005025125631
    514.854111405836        0.0364321608040205
    518.832891246684        0.0502512562814074
    524.137931034483        0.0753768844221109
    531.299734748011         0.131909547738694
    539.522546419098         0.253768844221106
    546.153846153846         0.409547738693468
    551.724137931034         0.576633165829146
    555.702917771883         0.703517587939699
    558.885941644562         0.802763819095478
    562.068965517241         0.880653266331659
    563.129973474801         0.903266331658292
    566.57824933687          0.964824120603015
    568.435013262599         0.987437185929649
    569.761273209549          0.99748743718593
    571.087533156499                         1
    572.413793103448          0.99497487437186
    574.005305039788         0.986180904522614
    575.596816976127         0.968592964824121
    578.779840848806         0.915829145728644
    583.023872679045          0.81532663316583
    585.411140583554         0.742462311557789
    588.063660477454         0.658291457286433
    591.246684350133         0.569095477386935
    594.429708222812         0.492462311557789
    599.20424403183          0.384422110552764
    603.978779840849         0.300251256281407
    610.079575596817         0.224874371859297
    618.302387267905         0.155778894472362
    627.055702917772          0.10427135678392
    636.604774535809        0.0640703517587943
    646.684350132626        0.0364321608040205
    659.151193633952        0.0163316582914576
    674.005305039788       0.00502512562814102
    751               0
    752               0
    753               0
    754               0
    755               0
    ];
%Poner los datos en valores absolutos:
%sigmapico = 1.5e-20;      % Valor al que hay que normalizar
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

sigmaemi = sigmaemi * 0.741042; % Fuchtbauer-Ladenburg correction

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
