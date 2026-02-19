function [y, H] = mittpunkt(a, N, T)
    
    % Steglängd
    h = T/N;
    
    % Allokera uttrymme + begynnelsevillkor
    y = zeros(4,N);
    y(1,1) = 1-a;
    y(2,1) = 0;
    y(3,1) = 0;
    y(4,1) = sqrt((1+a)/(1-a));

    % Energi
    H = zeros(1,N+1);
    
    % Options
    maxiter = 10;
    tol = 1e-5;
    
    for i = 1:N

        % Startgissning med explicit euler
        v = y(:,i);

        % Spara konstant startpunkt
        v_i = v;
    
        % Beräkna energi
        H(i) = (y(3,i)^2 + y(4,i)^2)/2 - 1/sqrt(y(1,i)^2 + y(2,i)^2);
        for j = 1:maxiter % Iterera newtons metod till nollställe r = 0

            % Ställ upp förenklingar
            Q1 = (y(1,i) + v(1))/2; 
            Q2 = (y(2,i) + v(2))/2;
            P1 = (y(3,i) + v(3))/2;
            P2 = (y(4,i) + v(4))/2;
            R = sqrt(Q1^2 + Q2^2);
    
            % residualer (systemekvationerna)
            r = v - v_i - h * [P1; P2; -Q1/R^3; -Q2/R^3];
    
            % Jakobian 
            J = zeros(4,4);
            J(1,1) = 1; 
            J(1,3) = -h/2;
    
            J(2,2) = 1;
            J(2,4) = -h/2;
    
            J(3,1) = h/2 * (R^2 - 3*Q1^2)/R^5;
            J(3,2) = - h/2 * 3 * Q1 * Q2 / R^5;
            J(3,3) = 1;
    
            J(4,1) = - h/2 * 3 * Q1 * Q2 / R^5;
            J(4,2) = h/2 * (R^2 - 3*Q2^2)/R^5;
            J(4,4) = 1;
    
            s = J \ r;
            v = v - s;
        
            if norm(s) < tol
                if j >= 4
                    disp(j)
                end
                y(:,i+1) = v;
                break
            end

            if j == maxiter
                disp("maxiter reached")
                y(:,i+1) = v;
            end
        end
    end
    % Returnera energi, q, p
    H(i+1) = (y(3,i+1)^2 + y(4,i+1)^2)/2 - 1/sqrt(y(1,i+1)^2 + y(2,i+1)^2);
end