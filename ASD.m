%Adaptive Steepest Descent Algorithim

%A -> 100x100 diagonal matrix of eigenvalues from 0.1 to 100
%x -> 100x1 vector of zeroes
%b -> 100x1 vector of ones
%delta -> constant to shorten steepest descent alpha
%kappa -> constant used as condition to alternate alphas
%g -> search direction
%alpha -> minimizes x multiplied by g to find x

%Set up A, x0, b
A = ones(100,100);
m = 2;
%Converts A to be diagonal with values 0.1 -> 100
for i = 1:100
    for j = 1:100
        if i == 1 && j == 1
            A(i,j) = 0.1;
        end
        %not equals in matlab is ~=
        if i ~= j
            A(i,j) = 0;
        end
        if (i == j) && i ~= 1
            A(i,j) = m;
            m =  m + 1;
        end
    end
end
x = zeros(100,1);
b = ones(100,1);

%set delta and kappa (constants)
delta = 0.5;
kappa = 0.5;

%set g0 and calculate alpha initials
g0 = A*x - b;

alphaSD = (g0'*g0)/(g0'*A*g0);
alphaMG = (g0'*A*g0)/(g0'*A*A*g0);

%first iteration done outside to be used in while condition
if alphaMG/alphaSD > kappa
    alpha = alphaMG;
else 
    alpha = alphaSD - delta*alphaMG;
end

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
while norm(g0)*10^-6 <= norm(g)
    iter = iter+1;
    %calculates alphas using new g
    alphaSD = (g'*g)/(g'*A*g);
    alphaMG = (g'*A*g)/(g'*A*A*g);

    if alphaMG/alphaSD > kappa
        alpha = alphaMG;
    else 
        alpha = alphaSD - delta*alphaMG;
    end

    x = x - alpha*g;
    g = A*x - b;
    %prints the norm and iteration number
    normofg = norm(g);
    fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
end