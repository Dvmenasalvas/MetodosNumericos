disp('Resolveremos un sistema lineal con el método de Jacobi por puntos');
n = input('Dimensión de la matriz A:');

A = zeros(n, n);
for i = 1:n
    m = ['Introduzca la fila ', num2str(i), ':'];
    A(i, 1:n) = input(m);
end            
b = zeros(1,n);
b = input('Introduzca b:');

iter = input('Número máximo de iteraciones:');
precision = input('Precisión buscada(epsilon):');       
[sol, preciso] = jacobiPunt(A, b, iter, precision);

disp('La solución es: ')
disp(sol')
if preciso == 1
    disp('Se ha alcanzado la precisión buscada');
else 
    disp('No se ha alcanzado la precisión buscada');
end

function [sol, preciso] = jacobiPunt(A, b, iter, precision)
    preciso = 0;
    n = size(A,1);
    %Norma 2 de b
    normb = norm(b);
    %Vectores con la notación del libro
    r = zeros(n, 1);
    d = zeros(n, 1);
    u_k = zeros(n, 1);
    for k = 1:iter
        r = b' - A*u_k;
        for i = 1:n
            d(i) = r(i)/A(i,i);
        end
        u_k = u_k + d;
        if precision*normb > norm(r)
            %Si alcanzamos la precision buscada
            %lo indicamos y paramos las iteraciones
            preciso = 1;
            break;
        end
    end
    sol = u_k;
end
