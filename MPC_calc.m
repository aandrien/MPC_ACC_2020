function [u, utilde, status, iA1,  Pout,rout, V] = MPC_calc(x, xref, mpc_params)
%#codegen
%% Retrieve parameters
Atil = mpc_params.Atil;
Btil = mpc_params.Btil;
Q = mpc_params.Q;
R = mpc_params.R;
N = mpc_params.N;
Deltatil = mpc_params.Deltatil;
Delta = mpc_params.Delta;
Delta_unscaled = mpc_params.Delta_unscaled;
c7 = mpc_params.c7;

%% Calculate P_epsilon (solution to DARE with scheduling epsilon)
r = 0.999;
% P = dare(Atil,B,r*Q,R);
P = 1000*eye(3);
% for k = 1:500
%     P = Atil'*P*Atil+r*Q-Atil'*P*B*inv(R+B'*P*B)*B'*P*Atil;
% end
while(x'*P*x*trace(P)>R*Delta^2/(2*trace(Btil*Btil')))
    r = 0.9*r;
    for k = 1:50
        P = Atil'*P*Atil+r*Q-Atil'*P*Btil*inv(R+Btil'*P*Btil)*Btil'*P*Atil;
    end
end
Pout = P;
rout = r;

F = -inv(Btil'*P*Btil+R)*Btil'*P*Atil;
D = R*F;
C = Btil'*P*Btil;

% xr0 = xref(1:3);
S = sign(D*x);

%% Inequality constraints
Aineq = [-eye(N); % U_{l|k} <= Delta
    eye(N); % -Delta <= U_{l|k}
    S*[1,zeros(1,N-1)]; % 0<= sign(R*F*x)[u_{0|k}-F*x]
    -S*C*[1,zeros(1,N-1)]]; % 0<= sign(R*F*x)(2*R*F-B'*P*B*[u_{0|k}-F*x])

% Aineq = [-eye(N); % U_{l|k} <= Delta
%     eye(N)]; % -Delta <= U_{l|k}

% Aineq = [];

bineq = [-Deltatil; % U_{l|k} <= Delta
    -Deltatil; % -Delta <= U_{l|k}
    S*F*x; % 0<= sign(RFx)[u_{0|k}-Fx]
    -S*2*D*x-sign(D*x)*C*F*x]; % 0<= sign(R*F*x)(2*R*F-B'*P*B*[u_{0|k}-F*x])

% bineq = [-Deltatil; % U_{l|k} <= Delta
%     -Deltatil]; % -Delta <= U_{l|k}

% bineq = zeros(0,1);

%% Equality constraints
Aeq = []; % no equality constraints
beq = zeros(0,1); % no equality constraints

[uMPC, status, iA1] = mpcqpsolver(mpc_params.Linv, mpc_params.f*(mpc_params.Ftil*x-xref), Aineq, bineq, Aeq, beq, false(size(bineq,1),1), mpc_params.opt);
uMPC = uMPC.*Delta_unscaled;
if(status<= 0)
    utilde = F*x;
else
    utilde = uMPC(1);
end
% utilde = F*x;
V = x'*Pout*x;

c7 = mpc_params.c7;

u = utilde-x(3)/c7;
