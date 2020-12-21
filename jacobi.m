disp('Resolveremos un sistema lineal con el m�todo de Jacobi por puntos');
n = input('Dimensi�n de la matriz A:');

A = zeros(n, n);
for i = 1:n
    m = ['Introduzca la fila ', num2str(i), ':'];
    A(i, 1:n) = input(m);
end            
b = zeros(1,n);
b = input('Introduzca b:');

iter = input('N�mero m�ximo de iteraciones:');
precision = input('Precisi�n buscada(epsilon):');       
[sol, preciso] = jacobiPunt(A, b, iter, precision);

disp('La soluci�n es: ')
disp(sol')
if preciso == 1
    disp('Se ha alcanzado la precisi�n buscada');
else 
    disp('No se ha alcanzado la precisi�n buscada');
end

function [sol, preciso] = jacobiPunt(A, b, iter, precision)
    preciso = 0;
    n = size(A,1);
    %Norma 2 de b
    normb = norm(b);
    %Vectores con la notaci�n del libro
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
