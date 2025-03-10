%EXAMPLE 2 FOR PRESENTATION

%Adaptive Barzilai Borwein Algorithim, Ex_2

%A -> 100x100 diagonal matrix of eigenvalues from 0.1 to 100
%x -> 100x1 vector of zeroes
%b -> 100x1 vector of ones
%delta -> constant to shorten steepest descent alpha
%kappa -> constant used as condition to alternate alphas
%g -> search direction
%alpha -> minimizes x multiplied by g to find x


%Creates a matrix of the kappa values
kappa = linspace(0.01, 0.80, 80);
%Set I outside loops (identity 5000x5000)
I = eye(200, 200);
vals = zeros(1,80);
num = 0;
for kap = 1:80
    total = 0;
    ave = 0;
    for prob = 1:200 
    %Set up A, x0, b
    %A = QDQ'
        %w1
        w1 = rand(200,200);
        w1_norm = norm(w1);
        %unit_vec
        w1 = w1/w1_norm;
        
        %w2
        w2 = rand(200,200);
        w2_norm = norm(w2);
        %unit_vec
        w2 = w2/w2_norm;
        
        %w3
        w3 = rand(200,200);
        w3_norm = norm(w3);
        %unit_vec
        w3 = w3/w3_norm;
        
        Q = (I - 2*(w3*w3'))*(I - 2*(w2*w2'))*(I - 2*(w1*w1'));
        D = ones(200,200);
        %Converts A to be diagonal with values 0.1 -> 100
        for i = 1:200
            for j = 1:200
                %not equals in matlab is ~=
                if i ~= j
                    D(i,j) = 0;
                end
                if (i == j) && i ~= 1
                    %finds random number between 1 and 100
                    a = 1;
                    b = 1000;
                    D(i,j) = a + (b-a).*rand(1,1);
                end
                if i == 200 && j == 200
                    D(i,j) = 1000; %sets final element to be condition
                end
            end
        end
        
        A = Q*D*Q';
        x = zeros(200,1);
        a = -10;
        c =  10;
        %generates a random b between -10 and 10
        b = a + (c-a).*rand(200,1);
        %resets x as old x value is kept in next iteration
        %set delta and kappa (constants)
        delta = 0.5;
        
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
        %norms = [normofg];
        %prints the norm and iteration number
        %fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
        
        %Sets all of the new variable names that will be used in the loop
        g_old = g0;
        g_new = g;
        alpha = alpha0;
        %algorithim begins
        %continues iterating till we minimize our direction to a tol of 10^-6
        while norm(g0)*10^-5 <= norm(g_new) && iter < 10000
            iter = iter+1;
        
            %find alphaBB1 and alphaBB2
            alphaBB1 = alpha*((g_old'*g_old)/(g_old'*g_old - g_old'*g_new));
            alphaBB2 = alpha*((g_old'*g_old - g_old'*g_new)/(g_old'*g_old - 2*g_old'*g_new + g_new'*g_new));
        
            %condition for alternating alpha
            if alphaBB2/alphaBB1 < kappa(kap)
                alpha = alphaBB2;
            else 
                alpha = alphaBB1;
            end
        
            %calculate x
            x = x - alpha*g_new;
            
            %set and find g_old and g_new
            g_old = g_new;
            g_new = A*x - b;
        
            %prints the norm and iteration number
            normofg = norm(g_new);
            %norms = [norms normofg];
            %fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
        end
        total = total + iter;
        %fprintf('cond = %d, tol = %d, iterations = %d\n',conditions(it), tols(its), iter)
    end
    ave = total/200;
    vals(num+1) = ave;
    num = num+1;
    fprintf('ABB: kappa = %2d, ave iter = %8d\n', kappa(kap), ave);
end

%plots ABB

figure(1);
plot(kappa, vals, 'g-', DisplayName='ABB')
hold on

%Adaptive Steepest Descent Algorithim, Ex_2

%A -> 100x100 diagonal matrix of eigenvalues from 0.1 to 100
%x -> 100x1 vector of zeroes
%b -> 100x1 vector of ones
%delta -> constant to shorten steepest descent alpha
%kappa -> constant used as condition to alternate alphas
%g -> search direction
%alpha -> minimizes x multiplied by g to find x


%Creates a matrix of the kappa values
kappa = linspace(0.01, 0.80, 80);
%Set I outside loops (identity 5000x5000)
I = eye(200, 200);
vals = zeros(1,80);
num = 0;
for kap = 1:80
    total = 0;
    ave = 0;
    for prob = 1:200 
    %Set up A, x0, b
    %A = QDQ'
        %w1
        w1 = rand(200,200);
        w1_norm = norm(w1);
        %unit_vec
        w1 = w1/w1_norm;
        
        %w2
        w2 = rand(200,200);
        w2_norm = norm(w2);
        %unit_vec
        w2 = w2/w2_norm;
        
        %w3
        w3 = rand(200,200);
        w3_norm = norm(w3);
        %unit_vec
        w3 = w3/w3_norm;
        
        Q = (I - 2*(w3*w3'))*(I - 2*(w2*w2'))*(I - 2*(w1*w1'));
        D = ones(200,200);
        %Converts A to be diagonal with values 0.1 -> 100
        for i = 1:200
            for j = 1:200
                %not equals in matlab is ~=
                if i ~= j
                    D(i,j) = 0;
                end
                if (i == j) && i ~= 1
                    %finds random number between 1 and 100
                    a = 1;
                    b = 1000;
                    D(i,j) = a + (b-a).*rand(1,1);
                end
                if i == 200 && j == 200
                    D(i,j) = 1000; %sets final element to be condition
                end
            end
        end
        
        A = Q*D*Q';
        x = zeros(200,1);
        a = -10;
        c =  10;
        %generates a random b between -10 and 10
        b = a + (c-a).*rand(200,1);
        %resets x as old x value is kept in next iteration
        %set delta and kappa (constants)
        delta = 0.5;
        
        %set g0 and calculate alpha initials
        g0 = A*x - b;
        
        alphaSD = (g0'*g0)/(g0'*A*g0);
        alphaMG = (g0'*A*g0)/(g0'*A^2*g0);
        
        %first iteration done outside to be used in while condition
        if alphaMG/alphaSD > kappa(kap)
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
        
        %algorithim begins
        %continues iterating till we minimize our direction to a tol of 10^-6
        while norm(g0)*10^-5 <= norm(g) && iter <10000
            iter = iter+1;
            %calculates alphas using new g
            alphaSD = (g'*g)/(g'*A*g);
            alphaMG = (g'*A*g)/(g'*A^2*g);
        
            if alphaMG/alphaSD > kappa(kap)
                alpha = alphaMG;
            else 
                alpha = alphaSD - delta*alphaMG;
            end
        
            x = x - alpha*g;
            g = A*x - b;
            %prints the norm and iteration number
            normofg = norm(g);
        end
        total = total + iter;
        %fprintf('cond = %d, tol = %d, iterations = %d\n',conditions(it), tols(its), iter)
    end
    ave = total/200;
    vals(num+1) = ave;
    num = num+1;
    fprintf('ASD: kappa = %2d, ave iter = %8d\n', kappa(kap), ave);
end

%plots asd
plot(kappa, vals, 'b-', DisplayName='ASD')
hold on



%Barzilai Borwein, Ex_2

%A -> 100x100 diagonal matrix of eigenvalues from 0.1 to 100
%x -> 100x1 vector of zeroes
%b -> 100x1 vector of ones
%delta -> constant to shorten steepest descent alpha
%kappa -> constant used as condition to alternate alphas
%g -> search direction
%alpha -> minimizes x multiplied by g to find x


%Creates a matrix of the kappa values
kappa = linspace(0.01, 0.80, 80);
%Set I outside loops (identity 5000x5000)
I = eye(200, 200);
vals = ones(1,80);
    %Set up A, x0, b
    %A = QDQ'
        %w1
        w1 = rand(200,200);
        w1_norm = norm(w1);
        %unit_vec
        w1 = w1/w1_norm;
        
        %w2
        w2 = rand(200,200);
        w2_norm = norm(w2);
        %unit_vec
        w2 = w2/w2_norm;
        
        %w3
        w3 = rand(200,200);
        w3_norm = norm(w3);
        %unit_vec
        w3 = w3/w3_norm;
        
        Q = (I - 2*(w3*w3'))*(I - 2*(w2*w2'))*(I - 2*(w1*w1'));
        D = ones(200,200);
        %Converts A to be diagonal with values 0.1 -> 100
        for i = 1:200
            for j = 1:200
                %not equals in matlab is ~=
                if i ~= j
                    D(i,j) = 0;
                end
                if (i == j) && i ~= 1
                    %finds random number between 1 and 100
                    a = 1;
                    b = 1000;
                    D(i,j) = a + (b-a).*rand(1,1);
                end
                if i == 200 && j == 200
                    D(i,j) = 1000; %sets final element to be condition
                end
            end
        end
        
        A = Q*D*Q';
        x = zeros(200,1);
        a = -10;
        c =  10;
        %generates a random b between -10 and 10
        b = a + (c-a).*rand(200,1);
        %resets x as old x value is kept in next iteration
        %set delta and kappa (constants)
        delta = 0.5;
        
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
        %fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
        
        %Sets all of the new variable names that will be used in the loop
        g_old = g0;
        g_new = g;
        alpha = alpha0;
        %algorithim begins
        %continues iterating till we minimize our direction to a tol of 10^-6
        while norm(g0)*10^-5 <= norm(g_new) && iter < 10000
            iter = iter+1;
        
            %find alphaBB1 and alphaBB2
            alphaBB1 = alpha*((g_old'*g_old)/(g_old'*g_old - g_old'*g_new));
            alphaBB2 = alpha*((g_old'*g_old - g_old'*g_new)/(g_old'*g_old - 2*g_old'*g_new + g_new'*g_new));
        
            %alternates between each every iteration
            if mod(iter, 2) == 0
                alpha = alphaBB1;
            else
                alpha = alphaBB2;
            end
            %calculate x
            x = x - alpha*g_new;
            
            %set and find g_old and g_new
            g_old = g_new;
            g_new = A*x - b;
        
            %prints the norm and iteration number
            normofg = norm(g_new);
            %fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
        end
        total = total + iter;
        %fprintf('cond = %d, tol = %d, iterations = %d\n',conditions(it), tols(its), iter)
    fprintf('BB: kappa = %2d, ave iter = %8d\n', kappa(kap), ave);
    vals = vals*iter;
%plots bb
plot(kappa, vals, 'r-', DisplayName='BB');


%Alternating Step, Ex_2

%A -> 100x100 diagonal matrix of eigenvalues from 0.1 to 100
%x -> 100x1 vector of zeroes
%b -> 100x1 vector of ones
%delta -> constant to shorten steepest descent alpha
%kappa -> constant used as condition to alternate alphas
%g -> search direction
%alpha -> minimizes x multiplied by g to find x


%Creates a matrix of the kappa values
kappa = linspace(0.01, 0.80, 80);
%Set I outside loops (identity 5000x5000)
I = eye(200, 200);
vals = ones(1,80);
    %Set up A, x0, b
    %A = QDQ'
        %w1
        w1 = rand(200,200);
        w1_norm = norm(w1);
        %unit_vec
        w1 = w1/w1_norm;
        
        %w2
        w2 = rand(200,200);
        w2_norm = norm(w2);
        %unit_vec
        w2 = w2/w2_norm;
        
        %w3
        w3 = rand(200,200);
        w3_norm = norm(w3);
        %unit_vec
        w3 = w3/w3_norm;
        
        Q = (I - 2*(w3*w3'))*(I - 2*(w2*w2'))*(I - 2*(w1*w1'));
        D = ones(200,200);
        %Converts A to be diagonal with values 0.1 -> 100
        for i = 1:200
            for j = 1:200
                %not equals in matlab is ~=
                if i ~= j
                    D(i,j) = 0;
                end
                if (i == j) && i ~= 1
                    %finds random number between 1 and 100
                    a = 1;
                    b = 1000;
                    D(i,j) = a + (b-a).*rand(1,1);
                end
                if i == 200 && j == 200
                    D(i,j) = 1000; %sets final element to be condition
                end
            end
        end
        
        A = Q*D*Q';
        x = zeros(200,1);
        a = -10;
        c =  10;
        %generates a random b between -10 and 10
        b = a + (c-a).*rand(200,1);
        %resets x as old x value is kept in next iteration
        %set delta and kappa (constants)
        delta = 0.5;
        
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
        %fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
        
        %Sets all of the new variable names that will be used in the loop
        g_old = g0;
        g_new = g;
        alpha = alpha0;
        %algorithim begins
        %continues iterating till we minimize our direction to a tol of 10^-6
        while norm(g0)*10^-5 <= norm(g_new) && iter < 10000
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
           % fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
        end
        total = total + iter;
        %fprintf('cond = %d, tol = %d, iterations = %d\n',conditions(it), tols(its), iter)
    fprintf('AS: kappa = %2d, ave iter = %8d\n', kappa(kap), iter);
    vals = vals*iter;

%plots AS
plot(kappa, vals, 'c-', DisplayName='AS');

%Alternating Minimal, Ex_2

%A -> 100x100 diagonal matrix of eigenvalues from 0.1 to 100
%x -> 100x1 vector of zeroes
%b -> 100x1 vector of ones
%delta -> constant to shorten steepest descent alpha
%kappa -> constant used as condition to alternate alphas
%g -> search direction
%alpha -> minimizes x multiplied by g to find x


%Creates a matrix of the kappa values
kappa = linspace(0.01, 0.80, 80);
%Set I outside loops (identity 5000x5000)
I = eye(200, 200);
vals = ones(1,80);
    %Set up A, x0, b
    %A = QDQ'
        %w1
        w1 = rand(200,200);
        w1_norm = norm(w1);
        %unit_vec
        w1 = w1/w1_norm;
        
        %w2
        w2 = rand(200,200);
        w2_norm = norm(w2);
        %unit_vec
        w2 = w2/w2_norm;
        
        %w3
        w3 = rand(200,200);
        w3_norm = norm(w3);
        %unit_vec
        w3 = w3/w3_norm;
        
        Q = (I - 2*(w3*w3'))*(I - 2*(w2*w2'))*(I - 2*(w1*w1'));
        D = ones(200,200);
        %Converts A to be diagonal with values 0.1 -> 100
        for i = 1:200
            for j = 1:200
                %not equals in matlab is ~=
                if i ~= j
                    D(i,j) = 0;
                end
                if (i == j) && i ~= 1
                    %finds random number between 1 and 100
                    a = 1;
                    b = 1000;
                    D(i,j) = a + (b-a).*rand(1,1);
                end
                if i == 200 && j == 200
                    D(i,j) = 1000; %sets final element to be condition
                end
            end
        end
        
        A = Q*D*Q';
        x = zeros(200,1);
        a = -10;
        c =  10;
        %generates a random b between -10 and 10
        b = a + (c-a).*rand(200,1);
        %resets x as old x value is kept in next iteration
        %set delta and kappa (constants)
        delta = 0.5;
        
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
        %fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
        
        %algorithim begins
        %continues iterating till we minimize our direction to a tol of 10^-6
        while norm(g0)*10^-5 <= norm(g) && iter < 10000
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
          %  fprintf('iter = %2d  norm = %.6f\n', iter, normofg)
        end
        total = total + iter;
        %fprintf('cond = %d, tol = %d, iterations = %d\n',conditions(it), tols(its), iter)
    fprintf('AM: kappa = %2d, ave iter = %8d\n', kappa(kap), iter);
    vals = vals*iter;
%plots AM
plot(kappa, vals, 'm-', DisplayName= 'AM')
xlabel('\kappa')
ylabel('Average Number of Iterations')
legend
hold off;