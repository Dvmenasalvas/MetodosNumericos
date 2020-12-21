disp('Resolveremos varios sistemas lineales que compartan la misma matriz A simetrica usando el método de Cholesky');
n = input('Dimensión de la matriz A:');

A = zeros(n, n);
for i = 1:n
    m = ['Introduzca los elementos  A(',num2str(i),',',num2str(i),':', num2str(n),')',':'];
    A(i, i:n) = input(m);
    for j = i+1:n
        A(j,i) = A(i,j);
    end
end            

b = leerVector(size(A,1));

[u, A] = cholesky(A, b);

disp('La solución es: ')
disp(u')

while (true)
    b = input('¿Deseas resolver más sistemas con otro vector b (1/0)?');
    if b == 0
        break;
    end
    b = leerVector(size(A, 1));
    u = resolverCholesky(A, b);
    disp('La solución es: ')
    disp(u')
end

function v = leerVector(n)
    v = zeros(1,n);
    v = input('Introduzca b:');
end

%Resuelve sistemas por el método de Cholesky.
function [sol, cholmat] = cholesky(A, b) 
    n = size(A,1);
    ok = true;
    
    %Calculamos la matriz triangular superior de la factorización
    %(la denotaremos B)
    %Será la U de LU traspuesta y haciendo sqrt de su diagonal
    for i = 1:n
        A(i,i) = A(i,i) - dot(A(i, 1:i-1) , A(i, 1:i-1));
        if A(i,i) <= 0
            disp('No se puede hacer Cholesky')
            %Si no se puede hacer LU, no se puede hacer Cholesky
            ok = false;
            break;
        end
        A(i,i) = sqrt(A(i,i));

        for j = i+1:n
            A(j,i) = (A(i,j) - dot(A(i, 1:i-1) , A(j, 1:i-1)))/A(i,i);
        end
    end
    if ok
        %Una vez tenemos B, resolvemos
        %Bw = b ; B'u = w
        sol = resolverCholesky(A, b);
        cholmat = A;
    end
end

%Permite resolver varios sistemas con la misma matriz
%Resuelve B'w = b ; Bu = w con el método de pivotaje
function u = resolverCholesky(Amod, b) 
    n = size(Amod,1);

    w = zeros(n,1);
    
    %Bw = b
    for i = 1:n
        w(i) = (b(i) - Amod(i, 1:i-1) * w(1:i-1))/Amod(i, i);
    end

    %B'u = w
    u = zeros(n,1);
    for i = n:-1:1
        u(i) = (w(i) - dot(Amod(i+1:n, i) , u(i+1:n)))/Amod(i, i);
    end
end
