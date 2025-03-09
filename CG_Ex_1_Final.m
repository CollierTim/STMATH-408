%Conjugate Gradient, Ex_1

%A -> 100x100 diagonal matrix of eigenvalues from 0.1 to 100
%x -> 100x1 vector of zeroes
%b -> 100x1 vector of ones
%kappa -> constant used as condition to alternate alphas
%g -> search direction
%alpha -> minimizes x multiplied by g to find x

%Set cond matrix and tols matrix for iterating between different conditions
%and tolerances in the problem
conditions = [10^2, 10^3, 10^4, 10^5, 10^6];
tols = [10^-1, 10^-2, 10^-3, 10^-4, 10^-5];
%Set up A, x0, b
%A = QDQ'

%creates identity of desired size
I = eye(3000, 3000);

for it = 1:5
        %w1
        w1 = rand(3000,3000);
        w1_norm = norm(w1);
        %unit_vec
        w1 = w1/w1_norm;
        
        %w2
        w2 = rand(3000,3000);
        w2_norm = norm(w2);
        %unit_vec
        w2 = w2/w2_norm;
        
        %w3
        w3 = rand(3000,3000);
        w3_norm = norm(w3);
        %unit_vec
        w3 = w3/w3_norm;
        
        Q = (I - 2*(w3*w3'))*(I - 2*(w2*w2'))*(I - 2*(w1*w1'));
        D = ones(3000,3000);
        m = 2;
        %Converts A to be diagonal with values 0.1 -> 100
        for i = 1:3000
            for j = 1:3000
                %not equals in matlab is ~=
                if i ~= j
                    D(i,j) = 0;
                end

                if (i == j) && i ~= 1
                    %finds random number between 1 and 100
                    a = 1;
                    b = conditions(it);
                    D(i,j) = a + (b-a).*rand(1,1);
                end
                if i == 3000 && j == 3000
                    D(i,j) = conditions(it); %sets final element to be condition
                end
            end
        end
        
        A = Q*D*Q';
        x = zeros(3000,1);
        
        a = -10;
        c = 10;
        %generates a random b between -10 and 10
        b = a + (c-a).*rand(3000,1);
    for its = 1:5
        x = zeros(3000,1);
        r0 = b - A*x;
        g0 = r0;
        normg0 = norm(g0);
        k = 0;
        
        alpha0 = (r0'*r0)/(g0'*A*g0);
        x = x +alpha0*g0;
        r = r0 - alpha0*A*g0;
        
        beta = (r'*r)/(r0'*r0);
        g = r + beta*g0;
        k = k+1;
        
        normofg = norm(g);
        %prints the norm and iteration number
        %fprintf('iter = %2d  norm = %.6f\n', k, normofg)
        
        %sets variables for loop
        r_old = r;
        g_old = g;
        
        while normg0*tols(its) < normofg && k < 10000
           alpha = (r_old'*r_old)/(g_old'*A*g_old); 
           x = x + alpha*g_old;
           r = r_old - alpha*A*g_old;
        
           beta = (r'*r)/(r_old'*r_old);
           g = r + beta*g_old;
           k = k+1;
        
           normofg = norm(g);
           %fprintf('iter = %2d  norm = %.6f\n', k, normofg)
        
           %sets old variables
            r_old = r;
            g_old = g;
        end
        fprintf('cond = %d, tol = %d, iterations = %d\n',conditions(it), tols(its), k)
    end
end