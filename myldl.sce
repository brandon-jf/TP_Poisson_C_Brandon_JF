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
temps=zeros(10,1);
k=1;
for n = [10,20,30,40,50,60,70,80,90,100]
    for l = 1:100
        A = rand(n,n);
        for i = 1:n
            for j = i:n
                A(j,i) = A(i,j);
            end
            A(i,i) = 10*A(i,i); 
        end
        tic();
        myldl(A);
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
        for i = 1:n
            for j = i:n
                A(j,i) = A(i,j);
            end
            A(i,i) = 10*A(i,i);
        end
    S = A //sauvegarde de A
        conda(n) = cond(A);
    //Choix d'un vecteur solution
        xs= rand(n,1);
        b=A*xs;
        A = myldl(A);
        L1=tril(A);
        D1 = zeros(n,n);
        for i = 1:n
            L1(i,i) =1.0;
            D1(i,i) =A(i,i);
        end
        y = usolv(L1',b);
        z = zeros(n,1);
        for i = 1:n
            z(i) = y(i)/D1(i,i);
        end
        xc = lsolv(L1,z);
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
