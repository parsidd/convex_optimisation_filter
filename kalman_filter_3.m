N = 10; % number of data points
n = 2;  % number of state variables
x_0 = [100 100]';   % initial point

A = rand(n, n);     % generate a random matrix for next state
A = 0.98*A/(eigs(A, 1, "largestabs"));  % scale so that maximum eigenvalue
                                        % is 0.98. This is done for
                                        % stability as we are using a
                                        % steady state filter
C = eye(2);         % use Identity as scale factor for measurement so that
                    % it is easy to see the measurements on graph
W = eye(2);         % covariance matrix for process noise
w = randn(n, N);    % process noise
V = eye(n);         % covariance matrix for measurement noise
v = randn(n, N);    % measurement noise 

P = 1.5 * eye(2);   % use a larger estimate for initial error covariance
                    % matrix
x = [x_0];
x_n = x_0;
y = [C*x_0+v(:,1)];

% generate data
for i = 2:N
    x = [x A*x_n + w(:, i)];
    x_n = x(:, length(x));
    outlier = [1 1]';
    if(rand()<0.5)          % randomly use outliers. assign measurement to 
                            % noise without signal for outliers
        outlier = [0 0]';
    end
    y = [y (C*x_n).*outlier + v(:, i)];
end

x_prev = A*x_0;
x_kf = [x_0];
t = 0.001;

%kalman filter for t = 0.001
for i = 2:N
    cvx_begin
        variables x_est(n) v(n) z(n);
        minimize sum(v'*V*v + (x_est-x_prev)'*inv(P)*(x_est-x_prev)) + t*norm(z,1);
        subject to 
            y(:, 1) ==C*x_est + v + z;
    cvx_end
    x_prev = A*x_est;
    x_kf = [x_kf x_est];
    P = A*P*A' + W - A*P*C'*inv(C*P*C' + V)*C*P*A';
end
hold on
plot(x_kf(2, 1:N),x_kf(1, 1:N), "go-",'DisplayName',"Recovered state with t = 0.001")

x_prev = A*x_0;
x_kf = [x_0];
t = 0.9;
P = 1.5*eye(2);

%kalman filter for t = 0.1
for i = 2:N
    cvx_begin
        variables x_est(n) v(n) z(n);
        minimize sum(v'*V*v + (x_est-x_prev)'*inv(P)*(x_est-x_prev)) + t*norm(z,1);
        subject to 
            y(:, 1) ==C*x_est + v + z;
    cvx_end
    x_prev = A*x_est;
    x_kf = [x_kf x_est];
    P = A*P*A' + W - A*P*C'*inv(C*P*C' + V)*C*P*A';
end

plot(x_kf(2, 1:N),x_kf(1, 1:N), "c*-",'DisplayName',"Recovered states with t = 0.9")
% plot(y(2,1:N), y(1, 1:N), "k*-")
plot(x(2,:),x(1,:),"ro-",'DisplayName',"Actual state")
legend

