%Alternate Step Algorithim

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

%For alpha0 we have, alpha0 = alphaSD
alpha0 = (g0'*g0)/(g0'*A*g0);

%finds x1 and g1 for which g1 is ued to find the two alphas
x = x - alpha0*g0;
g = A*x - b;

%sets iteration to 1 for x1 and finds norm
iter = 1;
normofg = norm(g);

%prints the norm and iteration number
fprintf('iter = %2d  norm = %.6f\n', iter, normofg)

%Sets all of the new variable names that will be used in the loop
g_old = g0;
g_new = g;
alpha = alpha0;
%algorithim begins
%continues iterating till we minimize our direction to a tol of 10^-6
while norm(g0)*10^-6 <= norm(g_new)
    iter = iter+1;

    %find alphaBB1 and alphaSD
    alphaBB1 = alpha*((g_old'*g_old)/(g_old'*g_old - g_old'*g_new));
    alphaSD = (g_new'*g_new)/(g_new'*A*g_new);

    %alternates between each every iteration
    if mod(iter, 2) == 0
        alpha = alphaBB1;
    else
        alpha = alphaSD;
    end
    %calculate x
    x = x - alpha*g_new;
    
    %set and find g_old and g_new
    g_old = g_new;
    g_new = A*x - b;

    %prints the norm and iteration number
    normofg = norm(g_new);
    fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
end