function f = ev(dim, b, rc, theta11, theta30, theta31, theta32, theta33)
% Value function iteration by Rust

% control
grid = 1; % (2,571) Rust(1987) p.p.1015
iteration = 10000; % 10000
error = 0.001; % 0.001

% state space for milage
s = (0 : grid : dim+150)'; % x = s*2571 
n = size(s,1); % number of discretize mile space

% state space for decision
d = [0,1]'; % 0: non-replace, 1: replace
m = size(d,1); %number of dicision space


% set preference params Rust(1987) p.p.1022
%b = 0.999; % discount factor 
%rc = 10.896; % replacement cost 
%theta11 = 1.1732; % c = 0.001*theta11*x maintenance cost func.
%theta30 = 0.1191; % prob x:+0 no mile up
%theta31 = 0.5762; % prob x:+1 2571 mile up
%theta32 = 0.2868; % prob x:+2
%theta33 = 0.0158; % prob x:+3


% Value Function Iteration
V = ones(n-3,m); % initialize
Vtemp = ones(n,m); 
Vnew = ones(n-3,m); % initialize

for i = 1:iteration
% % j = 1 (no-replacement in the previous period)
 Vnew(:,1) = theta30*log(exp(-0.001*theta11*s(1:n-3) + b*Vtemp(1:n-3,1))... 
           + exp(-rc*ones(n-3,1) + b*Vtemp(1:n-3,2)))...
           +theta31*log(exp(-0.001*theta11*(s(1:n-3)+1*ones(n-3,1)) + b*Vtemp(2:n-2,1))...
             + exp(-rc*ones(n-3,1) + b*Vtemp(2:n-2,2)))...
           +theta32*log(exp(-0.001*theta11*(s(1:n-3)+2*ones(n-3,1)) + b*Vtemp(3:n-1,1))...
             + exp(-rc*ones(n-3,1) + b*Vtemp(3:n-1,2)))...
           +theta33*log(exp(-0.001*theta11*(s(1:n-3)+3*ones(n-3,1)) + b*Vtemp(4:n,1))...
             + exp(-rc*ones(n-3,1) + b*Vtemp(4:n,2)));
 
 % j = 2 (replacement is done in the previous period)
 Vnew(:,2) = log(exp(b*ones(n-3,1)*Vtemp(1,1))...
            + exp(-rc*ones(n-3,1) + b*ones(n-3,1)*Vtemp(1,2)));



 if max((max(abs(Vnew-V)))') < error;
     break;
 end;
 
 % value function update
 V = Vnew;
 Vtemp(1:n-3,:) = V;
 Vtemp(n-2:n,:) = V(n-3,1)*ones(3,2);
 
 
end;

Vfinal = V(1:dim,:);

f = Vfinal;

         
diary off