%% Mittpunktsmetoden och symplektisk euler



% ----- IMPLICITA MITTPUNKTSMETODEN -------------------------------------------
% y(n+1) = y(n) + h * f( [y(n) + y(n+1)]/2 ) 

% f är inte linjär i y. Vi behöver lösa ett olinjärt system varje tidssteg

% y(n+1)-y(n)-h*f( [y(n)+y(n+1)]/2 ) = 0: (4 ekvationer, 4 okända)

% y = (q1, q2, p1, p2), p = (p1, p2), q = (q1, q2)

% f(q1, q2, p1, p2) = (p1, p2, -q1/|q|^3, -q2/|q|^3) : |q|=sqrt(q1^2+q2^2)



% ----- Systemet ------

% Svårt att skriva ut hela systemet i matlab. 
% Vi kan dock ange systemets jakobian:

% J_11 = 1   J_12 = 0   J_13 = -h/2   J_14 = 0
% J_21 = 0   J_22 = 1   J_23 = 0   J_24 = -h/2

% Definiera Q1 = [q1(n+1)+q1(n)]/2, Q2 = [q2(n+1)+q2(n)]/2

% J_31 = - h * [Q1^2 - Q2^2/2]/[Q1^2 + Q2^2]^{5/2}
% J_42 = - h * [Q2^2 - Q1^2/2]/[Q1^2 + Q2^2]^{5/2}

% J_32 = - h * 3/2 * Q1 * Q2 / [Q1^2 + Q2^2]^{5/2}
% J_41 = - h * 3/2 * Q1 * Q2 / [Q1^2 + Q2^2]^{5/2}

% J_33 = 1  J_34 = 0
% J_43 = 0  J_44 = 1

clear;

% parametrar
a = 0.5; 
N = 1000; % antal delintervall
T = 100; % slutpunkt
h = T/N; % steglängd

q1 = zeros(1,N); q1(1) = 1-a;
q2 = zeros(1,N); q2(1) = 0;
p1 = zeros(1,N); p1(1) = 0;
p2 = zeros(1,N); p2(1) = sqrt(1+a)/sqrt(1-a);

maxiter = 20;
tol = 1e-5;
for i = 1:N

    % Vid varje steg ska vi lösa: system = 0.
    % Vi gör startgissning y(n+1)=y(n) 

    v = [q1(i); q2(i); p1(i); p2(i)]; % Startgissning. Variabel.
    y_i = [q1(i); q2(i); p1(i); p2(i)]; % Konstant.

    for j = 1:maxiter
        
        % Ställ upp förenklingar
        Q1 = (q1(i) + v(1))/2;
        Q2 = (q2(i) + v(2))/2;
        P1 = (p1(i) + v(3))/2;
        P2 = (p2(i) + v(4))/2;
        R = sqrt(Q1^2 + Q2^2);

        % residualer (systemekvationerna)
        r = v - y_i - h * [P1; P2; -Q1/R^3; -Q2/R^3];

        % Jakobian 
        J = zeros(4,4);
        J(1,1) = 1; 
        J(1,3) = -h/2;
        J(2,2) = 1;
        J(2,4) = -h/2;

        J(3,1) = - h * (Q1^2 - Q2^2/2)/R^5;
        J(3,2) = - h * 3/2 * Q1 * Q2 / R^5;
        J(3,3) = 1;

        J(4,1) = - h * 3/2 * Q1 * Q2 / R^5;
        J(4,2) = - h * (Q2^2 - Q1^2/2)/R^5;
        J(4,4) = 1;

        s = J \ r;
        v = v - s;
    
        if norm(s) < tol
            q1(i+1) = v(1);
            q2(i+1) = v(2);
            p1(i+1) = v(3);
            p2(i+1) = v(4);
            break
        end
        if j == maxiter
            q1(i+1) = v(1);
            q2(i+1) = v(2);
            p1(i+1) = v(3);
            p2(i+1) = v(4);
        end
    end
end
subplot(1,2,1);
plot(q1, q2);
xlabel("q1");
ylabel("q2");
title("Implicita mittpunktsmetoden")
xlim([-1.5, .5])
ylim([-1, 1])









% ------- EULERS SYMPLEKTISKA METOD -----------

% p(n+1) = p(n) - h * q(n) / |q(n)|^3
% q(n+1) = q(n) + h * p(n+1)

% clear;
% 
% % parametrar
% a = 0.5; 
% N = 4000; % antal delintervall
% T = 100; % slutpunkt
% h = T/N; % steglängd

q = zeros(2,N); q(1,1) = 1-a; q(2,1) = 0;
p = zeros(2,N); p(1,1) = 0; p(2,1) = sqrt(1+a)/sqrt(1-a);

for n = 1:N
    R_q = (q(1,n)^2 + q(2,n)^2)^(3/2);
    p(:,n+1) = p(:,n) - h * q(:,n) / R_q;
    q(:,n+1) = q(:,n) + h * p(:,n+1);
end
subplot(1,2,2);
plot(q(1,:), q(2,:))
xlabel("q1")
ylabel("q2")
title("Eulers symplektiska metod")
xlim([-1.5, .5])
ylim([-1, 1])