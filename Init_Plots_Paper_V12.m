%% Default code
clc; close all; clear;

scale = false;
base = true;

%% Simulation parameters
dt = 1/1000;
tend = 16;

%% Quadcopter parameters
m = 1;
quad_params.m = m;
quad_params.g = 9.81;
c = 0.5;
quad_params.c = c;


tau = 1/20; %discretization step MPC
dt_MPC = tau;

% Model
Ac = [0 1 0;
    0 -c/m 1/m;
    0 0 0];

Bc = [0;0;1];

sysc = ss(Ac,Bc,[],[]);
sysd = c2d(sysc,tau);
Ad = sysd.A; Bd = sysd.B;
nx = size(Ad,1);
nu = size(Bd,2);

c1 = Ad(1,2); c2 = Ad(1,3); c3 = Bd(1,1);
c4 = Ad(2,2); c5 = Ad(2,3); c6 = Bd(2,1);
c7 = Bd(3,1);

Atil = [1, c1, c2-c3/c7;
    0, c4, c5-c6/c7;
    0, 0, 0];

% Atil = Ad;
ctrbAtil = rank(ctrb(Atil,Bd));
if(ctrbAtil~=nx)
    disp('Error, system not controllable')
    return
end

L = 2.5;
quad_params.L = L;
Delta = L/c7;

if(scale)
    Btil = Bd.*Delta;
    Delta_unscaled = Delta;
    Delta =1;
else
    Btil = Bd.*1;
    Delta_unscaled = 1;
end

% MPC settings
N = 50; % MPC horizon
Q = 1.*diag([1,1,1]); % State Cost
% Q = 1*eye(3);
% R = 1; % Input Cost
R = 0.001;
nconstr = 2*N+2;

%% Reference
amp = -1;
freq = 1/6;
[p_ref,v_ref,a_ref] = create_reference_trajectories(tend,dt,amp,freq);



%% Saturated PD controller settings
% Kp = 6;
% Kv = 1.5;

Kp = 1.5;
Kv = 6;

%% Initial conditions
p0 = [-10;5;3];
v0 = [-3;-1;2];
a0 = -c*v0;

x0 = [p0(1);v0(1);a0(1)];
y0 = [p0(2);v0(2);a0(2)];
z0 = [p0(3);v0(3);a0(3)];

% p0 = p_ref(1,2:4);
% v0 = v_ref(1,2:4);
% a0 = a_ref(1,2:4);

r30     = [0;0;1];
r30     = r30./norm(r30,2);

%% MPC
%Settings
[x_ref_MPC,y_ref_MPC,z_ref_MPC] = create_reference_MPC(p_ref,v_ref,a_ref);


Ftil = eye(nx);
G = zeros(nx,nu);
for k = 1:N
    Ftil = [Ftil;Atil^k];
    G = [G;Atil^(k-1)*Btil];
end
Gtilde = [];
for k = 0:N-1
    Gtilde = [Gtilde,[zeros(nx*k,nu);G(1:end-k*nx,:)]];
end

% Cost
Qtilde = blkdiag(kron(eye(N+1),Q));
% Qtilde = zeros((N+1)*nx,(N+1)*nx);
% gamma = 1.1;
% for k = 1:N+1
%     Qtilde(1+(k-1)*nx:k*nx,1+(k-1)*nx:k*nx) = gamma^k*Q;
% end



% Input cost (necessary because else problem ill posed?)
Rtilde = kron(eye(N),R);

% Write in standard QP form
H = 2*(Gtilde'*Qtilde*Gtilde+Rtilde);
[Ltmp,p] = chol(H,'lower');
if(p==0)
    Linv = Ltmp\eye(size(H,1));
else
    disp('H not pos def!')
end
if (~(istril(Linv) && all(diag(Linv)>0)))
    disp('H not pos def!')
end
ftilde = 2*Gtilde'*Qtilde;


% Inequality constraints
Deltatil = repmat(Delta,N,1);

% MPC options
opt = mpcqpsolverOptions;
opt.FeasibilityTol = 1e-8;
opt.MaxIter = 500;

mpc_params.Linv = Linv;
mpc_params.f = ftilde;
mpc_params.Atil = Atil;
mpc_params.Ftil = Ftil;
mpc_params.Deltatil = Deltatil;
mpc_params.Delta = Delta;
mpc_params.Delta_unscaled = Delta_unscaled;
mpc_params.A = Ad;
mpc_params.L = L;
mpc_params.opt = opt;
mpc_params.H = H;
mpc_params.N = N;
mpc_params.dt = dt_MPC;
mpc_params.x_ref = x_ref_MPC;
mpc_params.y_ref = y_ref_MPC;
mpc_params.z_ref = z_ref_MPC;
mpc_params.tol = 0;
mpc_params.Qtilde = Qtilde;
mpc_params.R = R;
mpc_params.Q = Q;
mpc_params.B = Bd;
mpc_params.Btil = Btil;
mpc_params.nx = nx;
mpc_params.c7 = c7;
mpc_params.base = base;

% Initialize P
r = 0.9;
P = dare(Atil,Btil,r*Q,R);
while(x0'*P*x0*trace(P)>Delta/(2*trace(Btil*Btil')))
    r = 0.99*r;
    P = dare(Atil,Btil,r*Q,R);
end
P0 = P
r

F = -inv(Btil'*P*Btil+R)*Btil'*P*Atil;
S = sign(R*F*x0);

% Calculate intitial input
xk = x0;
yk = y0;
zk = z0;
ux = MPC_calc(xk, zeros((N+1)*3,1), mpc_params);
uy = MPC_calc(yk, zeros((N+1)*3,1), mpc_params);
uz = MPC_calc(zk, zeros((N+1)*3,1), mpc_params);

aref0 = a0+c*v0

uMPC0 = [ux;uy;uz];

%% Functions
function [p_ref,v_ref,a_ref] = create_reference_trajectories(tend,dt,amp,freq)
% Create sinusoidal reference trajectories
N = tend/dt+1;
tvec = [0:dt:tend]';

% Position
x = zeros(N,1);
% y = amp*sin(2*pi*freq.*tvec);
y = zeros(N,1);
z = zeros(N,1);

% Velocity
vx = zeros(N,1);
% vy = amp*2*pi*freq*cos(2*pi*freq.*tvec);
vy = zeros(N,1);
vz = zeros(N,1);

% Acceleration
ax = zeros(N,1);
% ay = -amp*(2*pi*freq)^2*sin(2*pi*freq.*tvec);
ay = zeros(N,1);
az = zeros(N,1);

% Combine
p_ref = [tvec,x,y,z];
v_ref = [tvec,vx,vy,vz];
a_ref = [tvec,ax,ay,az];
end

function [x_ref,y_ref,z_ref] = create_reference_MPC(p_ref,v_ref,a_ref)
n = size(p_ref,1);
x_ref = zeros(n*3,1);
y_ref = zeros(n*3,1);
z_ref = zeros(n*3,1);
for k = 1:size(p_ref,1)
    x_ref(1+(k-1)*3:k*3,1) = [p_ref(k,2);v_ref(k,2);a_ref(k,2)];
    y_ref(1+(k-1)*3:k*3,1) = [p_ref(k,3);v_ref(k,3);a_ref(k,3)];
    z_ref(1+(k-1)*3:k*3,1) = [p_ref(k,4);v_ref(k,4);a_ref(k,4)];
end
end