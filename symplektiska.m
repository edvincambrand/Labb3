function [y, H] = symplektiska(a, N, T)
    % ------- EULERS SYMPLEKTISKA METOD -----------
    
    % p(n+1) = p(n) - h * q(n) / |q(n)|^3
    % q(n+1) = q(n) + h * p(n+1)
    
    % clear;
    % 
    % % parametrar
    % a = 0.5; 
    % N = 4000; % antal delintervall
    % T = 100; % slutpunkt


    h = T/N; % steglängd
    
    q1 = zeros(N,1); q1(1) = 1-a;
    q2 = zeros(N,1); q2(1) = 0;
    p1 = zeros(N,1); p1(1) = 0;
    p2 = zeros(N,1); p2(1) = sqrt((1+a)/(1-a));
    
    H = zeros(N,1); % Dokumentera energin för mittpunkt
    H(1) = -0.5;
    
    % Iterera!
    for i = 1:N
        % Energi:
        H(i+1) = 1/2 * (p1(i)^2 + p2(i)^2) - 1/sqrt(q1(i)^2 + q2(i)^2);
    
        R = (q1(i)^2 + q2(i)^2)^(3/2); % R = |q|^3

        p1(i+1) = p1(i) - h * q1(i) / R;
        p2(i+1) = p2(i) - h * q2(i) / R;
        q1(i+1) = q1(i) + h * p1(i+1);
        q2(i+1) = q2(i) + h * p2(i+1);
    end
    
    y = [q1(:), q2(:), p1(:), p2(:)];
    
end
