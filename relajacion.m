disp('Resolveremos un sistema lineal con el método de relajación por puntos');
n = input('Dimensión de la matriz A:');

A = zeros(n, n);
for i = 1:n
    m = ['Introduzca la fila ', num2str(i), ':'];
    A(i, 1:n) = input(m);
end            
b = zeros(1,n);
b = input('Introduzca b:');

w = input('Introduzca el parámetro de relajación: ');
iter = input('Número máximo de iteraciones:');
precision = input('Precisión buscada(epsilon):');
[sol, preciso] = relajacionPunt(A, b, iter, precision,w);

disp('La solución es: ')
disp(sol')
if preciso == 1
    disp('Se ha alcanzado la precisión buscada');
else 
    disp('No se ha alcanzado la precisión buscada');
end


function [sol, preciso] = relajacionPunt(A, b, iter, prec, w)
    preciso = 0;
    n = size(A,1);
    normb = norm(b);
    %Vectores con notación del libro
    r = zeros(n, 1);
    d = zeros(n, 1);
    u_k = zeros(n, 1);
    u_k1 = zeros(n, 1);
    for k = 1:iter
        %Inicializamos los primeros
        %Hacemos esta iteracion a parte para no tener problemas con
        %u_k1(1:-1)
        r(1) = b(1) - A(1, 1:n)*u_k(1:n);
        d(1) = w*(r(1)/A(1, 1));
        u_k1(1) = u_k(1) + d(1);
        for i = 2: n
            r(i) = b(i) - A(i,1:i-1)*u_k1(1:i-1) - A(i, i:n)*u_k(i:n);
            d(i) = w*(r(i)/A(i, i));
            u_k1(i) = u_k(i) + d(i);
        end
        if prec*normb > norm(r)
            preciso = 1;
            break;
        end
        u_k = u_k1;
    end
    sol = u_k1;
end
