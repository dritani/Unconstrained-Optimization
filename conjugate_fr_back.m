THEX = 0;
THEY = 0;
CONV = 0;
x = 1;
y = 1;
% These are simply here for quick reference
f = (x - 1)^2 + 100*(- x^2 + y)^2;
grad = [2*x - 400*x*(- x^2 + y) - 2 ; - 200*x^2 + 200*y];
p_k_formula =  [400*x*(- x^2 + y) - 2*x + 2;  200*x^2 - 200*y];
hess = [ 1200*x^2 - 400*y + 2, -400*x; -400*x,    200];

k=0;
limit = 100000;
tol = 0.01;
alpha = 0.0001;
theC = 10^(-4);
rho = 0.5;

x_old = 0;
y_old = 0;
x = 3;
y = 4;
p_k  = [400*x*(- x^2 + y) - 2*x + 2;  200*x^2 - 200*y];

while(1)
    epsilon = norm([2*x - 400*x*(- x^2 + y) - 2 ; - 200*x^2 + 200*y]);
    
    k=k+1;

    CONV(k) = epsilon;
    THEX(k) = x;
    THEY(k) = y;
    
    if (epsilon < tol)
        break;
    end  
    
    if (k==limit)
        break;
    end 
    
    % Backtracking Start
    f_new = (1-x)^2+100*((y-(x^2))^2);
    f_old = (1-x_old)^2+100*((y_old-(x_old^2))^2);
    grad_old = [2*x_old - 400*x_old*(- x_old^2 + y_old) - 2 ; - 200*x_old^2 + 200*y_old];
    p_k_old = [400*x_old*(- x_old^2 + y_old) - 2*x_old + 2;  200*x_old^2 - 200*y_old];
    RHS = f_old + theC * alpha * (grad_old' * p_k_old);
    
    tempx = x;
    tempy = y;
    while (f_new > RHS)
        alpha = alpha * rho;
        tempx = x_old + alpha * p_k_old(1);
        tempy = y_old + alpha * p_k_old(2);
        f_new = (1-tempx)^2+100*((tempy-(tempx^2))^2);
    end
    % Backtracking End
    
    x_old = x;
    y_old = y;

    x = x_old + alpha * p_k(1);
    y = y_old + alpha * p_k(2);
    r_k = [2*x_old - 400*x_old*(- x_old^2 + y_old) - 2 ; - 200*x_old^2 + 200*y_old];
    r_k1 = [2*x - 400*x*(- x^2 + y) - 2 ; - 200*x^2 + 200*y];
    
    % Fletcher-Reeves definition
    sigma = (r_k1'*r_k1)/(r_k'*r_k);
    
    p_k_old = p_k;
    p_k = -r_k1 + sigma*p_k_old;

end

x
y
k
epsilon
alpha

semilogy(CONV)
xlabel('Number of iterations')
ylabel('Norm of \nablaf(x,y)')
title('Conjugate Gradient FR - Convergence of Algorithm')
grid on
size = 8;
a = linspace(-size,size);
b = linspace(-size,size);
[A,B] = meshgrid(a,b);
C = (1-A).^2+100*((B-(A.^2)).^2);
levels = 100:100:100;
figure
hold on
contour(A,B,C,200)
plot(THEX,THEY)
plot(3,4,'g*')
plot(x,y,'r*')
legend('\nablaf(x,y)','Path','Start','End')
xlabel('X')
ylabel('Y')
title('Conjugate Gradient FR - Contour Plot and Algorithm Path')

dlmwrite('cg_CONV_back.txt',CONV);