//Implementation de la methode Jacobi
//Jean-Francois brandon
n = 6;
/*
A = zeros(n,n);
for i = 1:n
    A(i,i)= 10*rand();
    if i <= n-1
        A(i,i+1) = rand();
       A(i+1,i) = A(i,i+1);
    end
end
*/
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
function [X] = Jacobi(A,b)
    n = size(A,1);
    X = rand(n,1);
    r = norm(b-A*X)/norm(b);
    eps = 1e-15;
    k= 0
    while r > eps
        for i =1:n
            s = 0.0;
            for j  = 1:n
                if j != i
                    s = s + A(i,j)*X(j);
                end
                //disp(X(i));
            end
            X(i) = X(i) + (b(i) - s)/A(i,i);
        end 
    k = k +1;
    r = norm(b-A*X)/norm(b);
    end
    disp("k",k);
    disp("r",r)
endfunction
disp("x");
x = Jacobi(A,b);

disp(x)
