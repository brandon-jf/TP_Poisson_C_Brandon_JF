n = 4
disp('matrice random A');
A =rand(n,n);
disp(A);
function [A] = mylu(A)
    for k = 1:n-1
        for i=k+1:n
            A(i,k) = A(i,k)/A(k,k);
        end 
        for i = k+1:n
            for j = k+1:n
                A(i,j) = A(i,j) - A(i,k)*A(k,j);
            end
        end
    end
    U = triu(A);
endfunction
A = mylu(A);
L = tril(A);

U = triu(A);
for i=1:n
    L(i,i)=1.0;
end
disp('Matrix unit lower L')
disp(L);
disp('Matrice upper U');
disp(U);
disp('Resultat de LU');
disp(L*U);
