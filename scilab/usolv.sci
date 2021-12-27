function [x]= usolv(U,b)
    n = size(U,1);
    x = zeros(n,1);
    x(n) = b(n)/U(n,n);
    i=1;
    j=1;
    for i = 1:n-1
        x(n-i) = b(n-i);
        for j = n-i+1:n;
            x(n-i) = x(n-i)- x(j)*U(n-i,j);
        end
        x(n-i) = x(n-i)/U(n-i,n-i);
    end
endfunction
