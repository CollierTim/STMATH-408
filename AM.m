%Alternate Minimization Algorithim

%A -> 100x100 diagonal matrix of eigenvalues from 0.1 to 100
%x -> 100x1 vector of zeroes
%b -> 100x1 vector of ones
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
while norm(g0)*10^-6 <= norm(g)
    iter = iter+1;
    %calculates alphas using new g
    alphaSD = (g'*g)/(g'*A*g);
    alphaMG = (g'*A*g)/(g'*A*A*g);

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