//Jean-Francois Brandon
//Algorithme de factorisation A = LDL' en forme compacte 
n=4;
disp('A randonm');
A = rand(n,n);
disp(A);
for i = 1:n
    for j = i:n
        A(j,i) = A(i,j);
    end
    A(i,i) = 10*A(i,i);
end
disp('A symetrique');
disp(A);
function [A] =myldl(A)
    n = size(A,1);
    v = zeros(n,1);
    for j = 1:n do
        for i = 1:j-1 do
            v(i) =A(j,i)*A(i,i);
            //D(i,i) = v(i);
        end
        //A(j,j) = A(j,j) - A(j,1:j-1)*v(1:j-1);
        for k = 1:j-1
            A(j,j) = A(j,j) - A(j,k)*v(k);  
        end
        //A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n,1:j-1)*v(1:j-1))/v(j);
        for i = j+1:n
            for k = 1:j-1
                A(i,j) = (A(i,j) - A(i,k)*v(k));
            end
            A(i,j) = A(i,j)/A(j,j);
        end
    end
endfunction
A = myldl(A);
L1=tril(A);
D1 = zeros(n,n);
for i = 1:n
    L1(i,i) =1.0;
    D1(i,i) =A(i,i);
end

disp("Matrix unit lower L");
disp(L1);
disp("Matrice diagonal D");
disp(D1);
disp("resultat L*D*L =");
disp(L1*D1*(L1'));
