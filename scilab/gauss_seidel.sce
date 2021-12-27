//Methode de Gauss-Seidel
//Jean-Francois Brandon

A = rand(n,n);


for i = 1:n
    for j = i:n
        A(j,i) = A(i,j);
    end
    A(i,i) = 10*A(i,i);
end
disp("Matrice symetrique");
disp(A);
x_sol = rand(n,1);
b = A * x_sol;
disp("methode de Gauss-Seidel");
function [X] = gauss_seidel(A,b)
    n = size(A,1);
    X = rand(n,1);
    r = norm(b-A*X)/norm(b);
    eps = 1e-15;
    k= 0
    while r > eps
        for i =1:n
            s1 = 0.0;
            s2 = 0.0;
            for j  = 1:i
                if j != i-1
                    s1 = s1 + A(i,j)*X(j);
                end
                //disp(X(i));
            end
            for j  = i+1:n
                if j != i
                    s2 = s2 + A(i,j)*X(j);
                end
                //disp(X(i));
            end
            X(i) = X(i) + (b(i) - s1 -s2)/A(i,i);
        end 
    k = k +1;
    r = norm(b-A*X)/norm(b);
    end
    disp(k);
    disp(r)
endfunction

disp("x");
x = gauss_seidel(A,b);

disp(x)
