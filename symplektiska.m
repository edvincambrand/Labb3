function [y, H] = symplektiska(a, N, T)

    h = T/N; % steglängd
    
    % Allokera uttrymme + Begynnelsevillkor
    q1 = zeros(N+1,1); q1(1) = 1-a;
    q2 = zeros(N+1,1); q2(1) = 0;
    p1 = zeros(N+1,1); p1(1) = 0;
    p2 = zeros(N+1,1); p2(1) = sqrt((1+a)/(1-a));
    
    % Energi
    H = zeros(1,N+1); 
    
    % Iterera!
    for i = 1:N

        % Energi:
        H(i) = 1/2 * (p1(i)^2 + p2(i)^2) - 1/sqrt(q1(i)^2 + q2(i)^2);
        
        % "Radie" (i någon mening...)
        R = (q1(i)^2 + q2(i)^2)^(3/2); 

        % Iterera explicit euler för p
        p1(i+1) = p1(i) + h * ( - q1(i))/R;
        p2(i+1) = p2(i) + h * ( - q2(i))/R;
        
        % Iterera implicit euler för q
        q1(i+1) = q1(i) + h * p1(i+1);
        q2(i+1) = q2(i) + h * p2(i+1);
    end
    H(i+1) = 1/2 * (p1(i+1)^2 + p2(i+1)^2) - 1/sqrt(q1(i+1)^2 + q2(i+1)^2);
    y = [q1(:), q2(:), p1(:), p2(:)]';
end
