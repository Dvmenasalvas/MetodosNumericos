disp('Resolveremos varios sistemas lineales que compartan la misma matriz A invertible usando el método de factorización LU');
n = input('Dimensión de la matriz A:');

A = zeros(n, n);
for i = 1:n
    m = ['Introduzca la fila ', num2str(i), ':'];
    A(i, 1:n) = input(m);
end            

b = leerVector(size(A,1));

[u, A] = A_LU(A, b);

disp('La solución es: ')
disp(u')

while (true)
    b = input('¿Deseas resolver más sistemas con otro vector b (1/0)?');
    if b == 0
        break;
    end
    b = leerVector(size(A, 1));
    u = resolverLU(A, b);
    disp('La solución es: ')
    disp(u')
end

function v = leerVector(n)
    v = zeros(1,n);
    v = input('Introduzca b:');
end

%Resuelve sistemas por el metodo de factorizacion LU
%Calcula la factorización A=LU
function [sol, Amod] = A_LU(A, b) 
    n = size(A);
    n = n(1);
    ok = true;
    
    %Guardamos LU en A (omitimos la diagonal de L)
    for i = 1:n
        A(i,i) = A(i,i) - (A(i, 1:i-1)*A(1:i-1, i));
        %U(i,i) = A(i,i);
        if A(i,i) == 0
            %Si U(i,i) = 0, A no se puede factorizar LU
            disp('No se puede hacer LU')
            ok = false;
            break;
        end
        
        for j = i+1:n
            A(i,j) = A(i,j) - (A(i, 1:i-1)*A(1:i-1, j));
            %U(i,j) = A(i,j);
        end
        
        for j = i+1:n
            A(j,i) = (A(j,i) - (A(j, 1:i-1)*A(1:i-1, i)))/A(i,i);
            %L(j,i) = A(j,i);
        end
    end
    if ok
        %Una vez obtenido LU, resolvemos Lw = b, Uu=w
        %por el método del remonte
        sol = resolverLU(A, b);
        Amod = A;
    end
end

%Método del remonte para resolver Lw = b, Uu=w
%Permite resolver varios sistemas con la misma matriz
function u = resolverLU(Amod, b) 
    n = size(Amod,1);
    w = zeros(n,1);
    
    %Resuelve Lw = b
    for i = 1:n
        w(i) = b(i) - Amod(i, 1:i-1) * w(1:i-1);
    end

    %Resuelve Uu = w
    u = zeros(n,1);
    for i = n:-1:1
        u(i) = (w(i) - Amod(i, i+1:n) * u(i+1:n))/Amod(i, i);
    end
end
