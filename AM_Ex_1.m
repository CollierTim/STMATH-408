%Alternate Minimization Algorithim

%A -> 100x100 diagonal matrix of eigenvalues from 0.1 to 100
%x -> 100x1 vector of zeroes
%b -> 100x1 vector of ones
%g -> search direction
%alpha -> minimizes x multiplied by g to find x

%Set up A, x0, b
%A = QDQ'

%creates identity of desired size
I = eye(5000, 5000);

%w1
w1 = rand(5000,5000);
w1_norm = norm(w1);
%unit_vec
w1 = w1/w1_norm;

%w2
w2 = rand(5000,5000);
w2_norm = norm(w2);
%unit_vec
w2 = w2/w2_norm;

%w3
w3 = rand(5000,5000);
w3_norm = norm(w3);
%unit_vec
w3 = w3/w3_norm;

Q = (I - 2*(w3*w3'))*(I - 2*(w2*w2'))*(I - 2*(w1*w1'));

D = ones(5000,5000);
m = 2;
%Converts A to be diagonal with values 0.1 -> 100
for i = 1:5000
    for j = 1:5000
        %not equals in matlab is ~=
        if i ~= j
            D(i,j) = 0;
        end
        if i == 5000 && j == 5000
            D(i,j) = 5000;
        end
        if (i == j) && i ~= 1
            %finds random number between 1 and 100
            a = 1;
            b = 100;
            D(i,j) = a + (b-a).*rand(1,1);
        end
    end
end

A = Q*D*Q';
x = zeros(5000,1);

a = 10;
c = -10;
%generates a random b between -10 and 10
b = a + (c-a).*rand(5000,1);
%set g0 and calculate alpha initials
g0 = A*x - b;

alphaSD = (g0'*g0)/(g0'*A*g0);

%first iteration done outside to be used in while condition
alpha = alphaSD;

%finds x1 and g1
x = x - alpha*g0;
g = A*x - b;

%sets iteration to 1 for x1 and finds norm
iter = 1;
normofg = norm(g);

%prints the norm and iteration number
fprintf('iter = %2d  norm = %.6f\n', iter, normofg)

%algorithim begins
%continues iterating till we minimize our direction to a tol of 10^-6
while norm(g0)*10^-3 <= norm(g)
    iter = iter+1;
    %calculates alphas using new g
    alphaSD = (g'*g)/(g'*A*g);
    alphaMG = (g'*A*g)/(g'*A^2*g);

    if mod(iter, 2) == 0
        alpha = alphaMG;
    else 
        alpha = alphaSD;
    end

    x = x - alpha*g;
    g = A*x - b;
    %prints the norm and iteration number
    normofg = norm(g);
    fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
end