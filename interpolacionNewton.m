disp('Vamos a dibujar una función y su polinomio de interpolación de Lagrange en el intervalo (a,b).');
a = input('Introduzca a: ');
b = input('Introduzca b: ');
x = linspace(a, b);

fstr = input('Introduce una función de x con notación vectorial:', 's');
%Obtenemos valores de la función introducida
y = eval(fstr);
%La guardamos como una función
eval(['f = @(x) ', fstr]);

n = input('Introduce el número de puntos de interpolación equiespaciados en [a,b]: ');
h = (b - a) / n;
n = n + 1;
%En realidad tendremos un polinomio de grado <= n
%Interpolando n+1 puntos equidistantes, entre los que están a y b

p = a;
T = zeros(n, 2);
%Rellenamos la tabla con los valores de f en los puntos de
%interpolación
for i = 1:n
    T(i, 1) = p;
    T(i, 2) = f(p);
    p = p + h;
end    

hold on
[P, Pi, T] = newton(T);
plot(x, y);
hold off

b = 1;
while (b == 1)
    b = input('¿Deseas añadir otro punto de interpolación (1/0)?:');
    if (b == 1)
        x = input('Coordenada x del punto:');

        hold on
        [P, Pi, T] = nuevoPunto(x, f(x), P, Pi, T);
        hold off
    else
        break;
    end
end

function [P, Pi, T] = newton(T)
    % T es una matriz (n+1)x2
    % En la primera columna están x0, x1, ..., xn
    % En la segunda columna están f(x0), f(x1), ..., f(xn)
    
    n = size(T, 1) - 1; % número de puntos menos 1.
    
    Pi = 1;
    P = [];
    for k = 0:n-1      
        %Calculamos el polinomio de grado k
        %P_k = P_k-1 + Pi_k-1*DifDiv(k)
        P = [0 P] + Pi * T(1, 2);
        %Calculamos Pi_k
        Pi = [Pi 0] - [0, Pi .* T(k+1, 1)];
        
        %Calculamos en la segunda columna de T 
        %las diferencias divididas de orden k + 1
        %Para calcular el polinomio solo nos será util T(1,2)
        %El resto se necesitan para calcular las de ordenes mayores
        for i = 1:n-k
            T(i, 2) = (T(i, 2) - T(i+1, 2)) / (T(i, 1) - T(i+k+1, 1));
        end
    end
    
    P = [0 P] + Pi * T(1, 2);
    %Dejamos Pi calculado para poder añadir un nuevo punto
    Pi = [Pi 0] - [0, Pi .* T(n+1, 1)]; 
    r = linspace(min(T(:,1)), max(T(:,1)));
    %Pintamos el polinomio de interpolación obtenido
    plot(r, polyval(P, r));
end

function [P, Pi, T] = nuevoPunto(x, y, P, Pi, T)
    %Número de puntos menos 1 (ya que vamos a añadir uno nuevo)
    n = size(T, 1); 
    %Añadimos el nuevo punto a la tabla
    T = [T; [x, y]];
    
    %Calculamos diferencias divididas
    for i = n:-1:1
        T(i,2) = (T(i,2) - T(i+1,2))/(T(i,1) - T(n+1,1));
    end
    
    %Calculamos el nuevo polinomio
    P = [0, P] + Pi * T(1,2);
    %Dejamos Pi calculado para un posible nuevo punto
    Pi = [Pi, 0] - [0, Pi .* T(n+1,1)];
   
    %Pintamos el nuevo polinomio de interpolación
    r = linspace(min(T(:,1)), max(T(:,1)));
    plot(r, polyval(P, r));
end