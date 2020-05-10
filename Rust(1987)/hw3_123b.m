% hw3_123b
% Simulation of Harold Zuecher model and estimation
% Use function ev.m nLL.m

clear;
delete hw3_123b.out
diary hw3_123b.out

% data size
y = 3; % 1yr
T=52*y; % 1yr = 52w
N=10000;
dim = 175*y; % dimention for value function iteration

%% Value Function Iteration

% set preference params Rust(1987) p.p.1022
b = 0.999; % discount factor 
rc = 10.896; % replacement cost 
theta11 = 1.1732; % c = 0.001*theta11*x maintenance cost func.
theta30 = 0.1191; % prob x:+0 no mile up
theta31 = 0.5762; % prob x:+1 2571 mile up
theta32 = 0.2868; % prob x:+2
theta33 = 0.0158; % prob x:+3

% compute ev (Note that "EV(1,0)" is EV(0,1))
EV = ev(dim, b, rc, theta11, theta30, theta31, theta32, theta33);

disp('1. Value Function Iteration')
figure(1);
plot(EV(:,1));
hold on
plot(EV(:,2), 'r');
hold on
xlabel('milage (*2571)');
ylabel('EV (not normalized yer)')
title('EV')
hold off

%% data generation
% initialize X, D
X = zeros(T,N); % T*N
D = zeros(T,N); % decision i

c = 0; % counter

% For t=1
% Draw T1EV
y = rand(2,N);
e = 0.577 - log( -log(y));

% Draw milage 
p = rand(1,N);
for k=1:N
    if p(k) < theta30
        X(1,k) = 0;
    elseif p(k) <= theta31 + theta30
        X(1,k) = 1;
    elseif p(k) <= theta32 + theta31 + theta30
        X(1,k) = 2;
    elseif p(k) <= theta33 + theta32 + theta31 + theta30
        X(1,k) = 3;
    else % irregular
        X(1,k) = 1;
    end
  % Compute replacement decision
    if -0.001*theta11*X(1,k) + e(1,k) + b*EV(X(1,k)+1,1) < -rc + e(2,k) + b*EV(1,2)
        D(1,k) = 1;
        c = c+1;
    end
end

% For t,
for i=2:T
% Draw T1EV
y = rand(2,N);
e = 0.577 - log( -log(y));

% Draw milage 
p = rand(1,N);
for k=1:N
    if D(i-1,k) ==1
        X(i,k) = 0;
    elseif p(k) <= theta30
        X(i,k) = X(i-1,k);
    elseif p(k) <=  theta31 + theta30
        X(i,k) = X(i-1,k) + 1;
    elseif p(k) <= theta32 + theta31 + theta30
        X(i,k) = X(i-1,k) + 2;
    elseif p(k) <= theta33 + theta32 + theta31 + theta30
        X(i,k) = X(i-1,k) + 3;
    else % irregular
        X(i,k) = X(i-1,k) + 1;
    end
  % Compute replacement decision
    if -0.001*theta11*X(i,k) + e(1,k) + b*EV(X(i,k)+1,1) < -rc + e(2,k) + b*EV(1,2)
        D(i,k) = 1;
        c = c+1;
    end
end    
end

Num = c;

%% summary statistics (Rust(1987) pp. 1003?)
disp('2. Sample Summary statistics')
Rep = zeros(Num,2);
c = 1;
for i=1:T
for k=1:N
    if D(i,k) ==1
        Rep(c,:) = [X(i,k) i];
        c=c+1;
    end
end
end


% Milage
disp('Milage summary')
Mmax = max(Rep(:,1)');
Mmax = Mmax*2571
Mmin = min(Rep(:,1)');
Mmin = Mmin*2571
Mmean = mean(Rep(:,1));
Mmean = Mmean*2571
Mvar = var(Rep(:,1));
Mvar = Mvar*2571^2

% Month
disp('Week summary')
Wmax = max(Rep(:,2)')
Wmin = min(Rep(:,2)')
Wmean = mean(Rep(:,2))
Wvar = var(Rep(:,1))


%% 1st step estimation theta30-33
disp('start estimation.')
% initialization
c0 = 0; %counter for theta30
c1 = 0; %counter for theta31
c2 = 0; %counter for theta32
c3 = 0; %counter for theta33
ANum = T*N;


% For t=1
for k=1:N
    if D(1,k) == 0
     if X(1,k) == 0;
        c0 = c0 + 1;
     elseif X(1,k) == 1;
        c1 = c1 + 1;
     elseif X(1,k) == 2;
        c2 = c2 + 1;
     elseif X(1,k) == 3;
        c3 = c3 + 1;
     end
    end
end

% For t
for i = 2:T
for k=1:N
    if D(i,k) == 0
    
     if X(i,k) == X(i-1,k);
        c0 = c0 + 1;
     elseif X(i,k) == X(i-1,k) + 1;
        c1 = c1 + 1;
     elseif X(i,k) == X(i-1,k) + 2;
        c2 = c2 + 1;
     elseif X(i,k) == X(i-1,k) + 3;
        c3 = c3 + 1; 
     end
    end
end
end

disp('Estimation result for theta30-33')
theta30hat = c0/(ANum-Num)
theta31hat = c1/(ANum-Num)
theta32hat = c2/(ANum-Num)
theta33hat = c3/(ANum-Num)

%% 2nd step estimation

disp('3. Estimation result for rc and theta11')
beta = fminsearch(@(beta) nLL(dim, b, beta, theta30hat, theta31hat, theta32hat, theta33hat, X, D), [10.7; 1.2]);
rchat = beta(1)
theta11hat = beta(2)

diary off
