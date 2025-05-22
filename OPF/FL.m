function [x, f_vals, hi_vals, he_vals, KKT_gaps] = FL(Df, Dhi, Dhe, x0, eta, k, max_iters,Hess)
    % g is a (d+1)-dimensional function where g = {f, h}
    % x0 is the initial point
    % eta is the learning rate
    % k is a matrix for the constraint
    % epsilon is the finite difference step size
    % max_iters is the number of iterations
    assert(nargin >=7)

    x = x0;
    
    % Initialize history
    f_vals = zeros(max_iters, 1);
    hi_vals = zeros(max_iters, 1);
    he_vals = zeros(max_iters, 1);
    KKT_gaps = zeros(max_iters, 1);

    for iter = 1:max_iters
        % Estimate the Jacobian of g = {f, h}

        % Extract gradient of f and Jacobian of h
        [f_x, grad_f] = Df(x); % Transpose to get the column vector
        [hi_x, J_hi] = Dhi(x);
        [he_x, J_he] = Dhe(x);

        
        % Solve for lambda using the solve_lambda function
        


        % Update x using the gradient descent step
        if nargin == 7
            try
                [li,le] = solve_lambda(J_hi, J_he, grad_f, hi_x, he_x, k);
            catch
                a = a+1;
            end
            x = x - eta * (grad_f + J_hi' * li + J_he'* le); %TBM
        elseif nargin == 8
            Hinv = inv(Hess);
            try
                [li,le] = solve_lambda2(J_hi, J_he, grad_f, hi_x, he_x, Hinv);
            catch
                a = a+1;
            end
            delta = grad_f + J_hi' * li + J_he'* le;
            x = x-eta*Hinv*delta;
        end

        % Store history of f(x) and h(x)
        f_vals(iter) = f_x;
        hi_vio = max(0, max(hi_x));
        he_vio = max(abs(he_x));
        hi_vals(iter) = hi_vio;
        he_vals(iter) = he_vio;

        lag_grad = norm(grad_f + J_hi' * li + J_he'* le);
        KKT_gap = max([lag_grad, dot(hi_x, li), hi_vio, he_vio]);
        KKT_gaps(iter) = KKT_gap;
        disp(KKT_gap)

    end
end

function [li_opt,le_opt] = solve_lambda(J_hi, J_he, grad_f, hi_x, he_x, k)
    % Solve for lambda using Gurobi given the Jacobian of h (J_h),
    % the gradient of f (grad_f), the values of h(x) (h_x), and matrix K.
    
    % Input:
    % J_h - Jacobian of h(x), size (m, n)
    % grad_f - Gradient of f(x), size (n, 1)
    % h_x - Values of h(x), size (m, 1)
    % K - Hyperparameter matrix, size (m, m)
    
    mi = size(J_hi, 1); % Number of constraints
    me = size(J_he,1);
    h_x = [hi_x;he_x];

    J_h = [J_hi; J_he];

    % Objective function
    % 1/2 * lambda' * J_h * J_h' * lambda + lambda' * J_h * grad_f - lambda' * K * h_x
    H =  (J_h * J_h');
    H = 1e-5*sparse(eye(size(H,1))) + H;
    f =(J_h * grad_f-k * h_x);
   
    lb = [zeros(mi, 1); -inf(me,1)];
    
    options =  optimoptions('quadprog','Display','off');
    [lambda,fval,exitflag] = quadprog(H, f,[],[],[],[],lb, [],[],options);
    %disp(exitflag)
    
    % Extract the optimal lambda
    li_opt = lambda(1:mi);
    le_opt = lambda(mi+1:end);
end

function [li_opt,le_opt] = solve_lambda2(J_hi, J_he, grad_f, hi_x, he_x, Hinv)
    % Solve for lambda using Gurobi given the Jacobian of h (J_h),
    % the gradient of f (grad_f), the values of h(x) (h_x), and matrix K.
    
    % Input:
    % J_h - Jacobian of h(x), size (m, n)
    % grad_f - Gradient of f(x), size (n, 1)
    % h_x - Values of h(x), size (m, 1)
    % K - Hyperparameter matrix, size (m, m)
    
    mi = size(J_hi, 1); % Number of constraints
    me = size(J_he,1);
    h_x = [hi_x;he_x];

    J_h = [J_hi; J_he];

    % Objective function
    % 1/2 * lambda' * J_h * J_h' * lambda + lambda' * J_h * grad_f - lambda' * K * h_x
    H =  (J_h *Hinv* J_h');
    f = J_h * (Hinv*grad_f)- h_x;
   
    lb = [zeros(mi, 1); -inf(me,1)];
    
    options =  optimoptions('quadprog','Display','off');
    [lambda,fval,exitflag] = quadprog(H, f,[],[],[],[],lb, [],[],options);
    %disp(exitflag)
    
    % Extract the optimal lambda
    li_opt = lambda(1:mi);
    le_opt = lambda(mi+1:end);
end
