function [L, P] = lqr_infinite_horizon_solution(Q, R)

%% find the infinite horizon L and P through running LQR back-ups
%%   until norm(L_new - L_current, 2) <= 1e-4  
dt = 0.1;
mc = 10; mp = 2.; l = 1.; g= 9.81;

% TODO write A,B matrices
a1 = mp*g / mc;
a2 = (mc+mp)*g / (l*mc);

df_ds = [0 0 1 0; 
        0 0 0 1; 
        0 a1 0 0; 
        0 a2 0 0];
A = eye(4) + dt*df_ds;

df_du = [0; 0; 1/mc; 1/(l*mc)];
B = dt* df_du;

% TODO implement Riccati recursion
P_cur = Q; % Initialize P_{k+1}
L_cur = zeros(1,4); % Initialize L_{k+1}, for first break comparison
%count = 0;
while 1
    %count = count + 1
    % Iterate  
    L_new = -inv(R + B'*P_cur*B) * (B'*P_cur*A);
    P_new = Q + L_new'*R*L_new + (A+B*L_new)'*P_cur*(A+B*L_new);
    % Check convergence
    if norm(L_new - L_cur, 2) <= 1e-4
        break
    end
    % Update
    P_cur = P_new;
    L_cur = L_new;
end

P = P_new;
L = L_new;

end