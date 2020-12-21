disp('Resolveremos varios sistemas lineales que compartan la misma matriz A invertible usando el m�todo de eliminaci�n gaussiana');
n = input('Dimensi�n de la matriz A:');

A = zeros(n, n);
for i = 1:n
    m = ['Introduzca la fila ', num2str(i), ':'];
    A(i, 1:n) = input(m);
end            

b = leerVector(size(A,1));

[u, A, punt] = PA_LU(A, b);

disp('La soluci�n es: ')
disp(u')

while (true)
    b = input('�Deseas resolver m�s sistemas con otro vector b (1/0)?');
    if b == 0
        break;
    end
    b = leerVector(size(A, 1));
    u = resolverGauss(A, b, punt);
    disp('La soluci�n es: ')
    disp(u')
end

function v = leerVector(n)
    v = zeros(1,n);
    v = input('Introduzca b:');
end

%A partir de la matriz A, busca la descomposicion PA = LU
function [sol, Amod, punt] = PA_LU(A, b) 
    n = size(A,1);
    punt = 1:n;

    for i = 1:n-1
        %Buscaremos PA = LU
        %Guardaremos LU en A
        
        %Pivotaremos con el pivote parcial
        [m pos] = max(abs(A(punt(i:n), i)));

        %A([i, pos+i-1], :) = A([pos+i-1, i], :);
        %No permutamos las filas de A, 
        %solo indicamos la permutaciones en punt
        punt([i, pos+i-1]) = punt([pos+i-1, i]);

        %Guardamos los opuestos de los multiplicadores en la parte de 
        %la columna i que pertenece a L
        A(punt(i+1:n), i) = (A(punt(i+1:n), i)/A(punt(i), i));
        
        %Actualizamos la fila de U 
        %U es la matriz que queda tras aplicar el m�todo de Gauss a Ax = b
        %A cada A(j, i+1:n) se le suma A(j,i+1:n) por el multiplicador
        for j = i+1:n
            A(punt(j), i+1:n) = A(punt(j), i+1:n) - (A(punt(i), i+1:n)*A(punt(j), i));
        end
    end
    
    %Una vez tenemos PA = LU, procedemos a resolver el sistema
    % Lw = Pb ; Uu = w . Por el m�todo del remonte
    sol = resolverGauss(A, b, punt);
    Amod = A;
end

function u = resolverGauss(Amod, b, punt)
    %Resolvemos Lw = Pb ; Uu = w . Por el m�todo del remonte
    %A tiene la informaci�n de LU y punt la informaci�n de P
    n = size(Amod);
    n = n(1);
    %Inicializamos w como vector de 0's
    w = zeros(n,1);

    w(1) = b(punt(1));
    %Resolvemos Lw = Pb 
    for i = 2:n
        w(i) = b(punt(i)) - Amod(punt(i), 1:i-1) * w(1:i-1);
    end

    %Inicializamos w como vector de 0's
    u = zeros(n,1);
    
    u(n) = w(n)/Amod(punt(n), n);
    %Resolvemos Uu = w 
    for i = n-1:-1:1
        u(i) = (w(i)- Amod(punt(i), i+1:n) * u(i+1:n))/Amod(punt(i), i);
    end
end
