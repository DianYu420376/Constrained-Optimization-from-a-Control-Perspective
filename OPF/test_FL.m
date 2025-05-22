
n=100;
mi = 50;
me =20;
m= mi+me;
x0 = zeros(n,1); % Input vector [x1; x2]

eta = 0.1;
k = 1/eta * eye(m);
max_iters = 1000;
A = rand(mi, n); % Random matrix for constraints
b = 10 * ones(mi, 1);
C = rand(me, n); % Random matrix for equality constraints
d = zeros(me, 1); % Constraint targets (example: vector of 1s)



[x, f_vals, hi_vals, he_vals, KKT_gaps] = FL(@Df, @(x)Dhi(x,A,b), @(x)Dhe(x,C,d), x0, eta, k, max_iters);
plot(log(KKT_gaps))


function [f_val, grad_f] = Df(x)
    % Quadratic objective: f(x) = x'Qx + c'x
    n = 100; % Number of parameters
    Q = eye(n); % Quadratic coefficient matrix (identity for simplicity)
    c = (1:n)'; % Linear coefficient (example: [1, 2, ..., n]')
    f_val = 0.5 * x' * Q * x + c' * x; % Objective value
    grad_f = Q * x + c; % Gradient: df/dx
end

% Inequality Constraints hi(x)

function [hi_val, grad_hi] = Dhi(x,A,b)
    % 50 inequality constraints: hi(x) = Ax - b
     % Constraint bounds (example: vector of 10s)
    hi_val = A * x - b; % Constraint values
    grad_hi = A; % Jacobian: dhi/dx
end


% Equality Constraints he(x)
function [he_val, grad_he] = Dhe(x,C,d)
    % 20 equality constraints: he(x) = Cx - d
    he_val = C * x - d; % Constraint values
    grad_he = C; % Jacobian: dhe/dx
end

% Example Usage
