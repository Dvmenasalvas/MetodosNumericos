disp('Resolveremos un sistema lineal tridiagonal');
n = input('Dimensión de la matriz tridiagonal:');
    
a = zeros(1,n-1);
a = input('Introduzca la diagonal inferior:');
b = zeros(1,n);
b = input('Introduzca la diagonal principal:');
c = zeros(1,n-1);
c = input('Introduzca la diagonal superior:');
d = zeros(1,n);
d = input('Introduzca el vector de términos independientes:');

u = tridiag(a,b,c,d);

disp('La solución es: ')
disp(u)



function sol = tridiag(a,b,c,d) %Resuelve sistemas tridiagonales.
    n = length(d);
    
    m = zeros(1,n);
    g = zeros(1,n);

    m(1) = b(1);
    g(1) = d(1)/m(1);

    %Calculamos m y g
    for k = 2:n
        m(k) = b(k) - (c(k-1)*a(k-1))/m(k-1);
         g(k) = (d(k)-(g(k-1)*a(k-1)))/m(k);
    end

    sol = zeros(1,n);
    sol(n) = g(n);
    %Calculamos la solucion
    for k = n-1:-1:1
        sol(k) = g(k) - (sol(k+1)*c(k))/m(k);
    end
end

