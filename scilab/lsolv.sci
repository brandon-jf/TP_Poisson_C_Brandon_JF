function [x]= lsolv(L,b)
    n = size(L,1);
    x = zeros(n,1);
    x(1) = b(1)/L(1,1);
    i=1;
    j=1;
    for i = 2:n
        x(i) = b(i);
        for j = 1:i-1;
            x(i) = x(i)- x(j)*L(i,j);
        end
        x(i) = x(i)/L(i,i);
    end
endfunction
