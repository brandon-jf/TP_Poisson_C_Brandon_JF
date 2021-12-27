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

temps=zeros(10,1);
k=1;
for n = [10,20,30,40,50,60,70,80,90,100]
    for l = 1:100
        //Matrice quelquonque 
        A = rand(n,n);
        tic();
        mylu(A);
        temps(k)=temps(k) + toc()/100; 
    end
    k=k+1;
    disp(k);
end   
disp(temps);
vect = [10,20,30,40,50,60,70,80,90,100];
subplot(221);
plot2d('ll',vect,temps,3);
xtitle("duree d execution en fonction de la taille en echelle logarithmique ");
conda = zeros(99,1);
f_err = zeros(n,1);
b_err = zeros(n,1);
exec('lsolv.sci',-1);
exec('usolv.sci',-1);
for n = 2:100
    for k = 1:100
    //creation de matrice symetrique
        A = rand(n,n);
        S = A //sauvegarde de A
        conda(n) = cond(A);
    //Choix d'un vecteur solution
        xs= rand(n,1);
        b=A*xs;
        A = mylu(A);
        L = tril(A);
        U = triu(A);
        for i=1:n
            L(i,i)=1.0;
        end
        y = usolv(U,b);
        xc = lsolv(L,y);
        f_err(n) = f_err(n) + 0.01*norm(xs - xc)/ norm(xs);
        b_err(n) = b_err(n) + 0.01*norm(b - S*xc)/(norm(S)*norm(xc));
    end
end

subplot(223);
plot(1:100,f_err,3);
xtitle("forward error of usolv");
subplot(224);
plot(1:100,b_err,3);
xtitle("backward error of usolv");
