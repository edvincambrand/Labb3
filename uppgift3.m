%% Uppgift 3

% Lös keplerproblemet med ODE45. = Referenslösning
% Jämför Symplektisk/Mittpunkt för olika $h$, både vad gäller bevarande av energi och tillståndet $(q_1,q_2)$ i slutpunkten $t=T$.
% Vad blir noggrannhetsordningen av Symplektisk/Implicit?
% Givet tolerans i sluttillståndet (jmf med ode45), vilken metod är mest effektiv?
% Beror svaret på om du vill ha hög noggrannhet (location wise), eller nöjer du dig med mindre noggrann lösning?
clear;

% Konstanter & span
a = 0.5;
N = 4000;
T = 100;
h = T/N;
tspan = (0:h:T)';

function dydt = odefunc(t, y)
    dydt = zeros(4,1);
    dydt(1) = y(3);
    dydt(2) = y(4);
    dydt(3) = -y(1)/(y(1)^2+y(2)^2)^(3/2);
    dydt(4) = -y(2)/(y(1)^2+y(2)^2)^(3/2);
end

% Startparametrar
y0 = [1-a, 0, 0, sqrt((1+a)/(1-a))];
rtol = 1e-12;
atol = 1e-12;
options = odeset(RelTol=rtol, AbsTol=atol);

% Iterera över olika h, jämför olika sluttillstånd

Nvals = [1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000, 256000]; % -> h halveras.

% Referenslösning
[t, y] = ode45(@(t,y) odefunc(t,y), tspan, y0, options);

Hvals = zeros(length(Nvals), 2);
E_m = zeros(length(Nvals),1);
E_s = zeros(length(Nvals),1);

for i = 1:length(Nvals)
    disp("----------------------------------------------------")
    disp("Tid: Mittpunkt, Symplektiska:")
    N = Nvals(i);
    tic
    % Mittpunkt, symplektisk : lösningar
    [y1, H1] = mittpunkt(a, N, T);
    toc
    tic
    [y2, H2] = symplektiska(a, N, T);
    toc
    
    % Spara slutenergier
    Hvals(i, 1) = H1(end);
    Hvals(i, 2) = H2(end);

    % Skriv ut differenser
    
    fprintf("Sluttillstånd. N = %.0f, T = %.0f \n "+ ...
        "ODE45: q1 = %.6f | q2 = %.6f \n "+ ...
        "Mittp: q1 = %.6f | q2 = %.6f \n " +...
        "Sympl: q1 = %.6f | q2 = %.6f \n ", ...
        N, T, y(end, 1), y(end, 2), y1(end, 1), y1(end, 2), y2(end, 1), y2(end, 2));

    % bestäm norm av felen
    E_m(i) = sqrt((y(end, 1) - y1(end, 1))^2 +(y(end, 2) - y1(end, 2))^2);
    E_s(i) = sqrt((y(end, 1) - y2(end, 1))^2 +(y(end, 2) - y2(end, 2))^2);
end

disp("--- Total avvikelse i sluttillstånd q1, q2 ---")
for i = 1:length(Nvals)
    fprintf("N = %.0f | Err M = %.6f | Err S = %.6f \n", ...
        Nvals(i), E_m(i), E_s(i))
end

disp("--- Noggrannhetsordningar ---")
for i = 1:length(Nvals)-1
    fprintf("p_m = %.6f | p_s = %.6f \n", ...
        log(E_m(i)/E_m(i+1))/log(2), log(E_s(i)/E_s(i+1))/log(2))
end

disp("--- Slutenergier ---")
for i = 1:length(Nvals)
    fprintf("N = %.0f | H_m = %.8f | H_s = %.8f \n", ...
        Nvals(i), Hvals(i,1), Hvals(i,2))
end


%% Extra.

t_M = [0.005603, 0.011072, 0.017700, 0.034865, 0.061087, 0.113703, 0.228974, 0.476079, 0.929533];
t_S = [0.000952, 0.000471, 0.000732, 0.001831, 0.002924, 0.005312, 0.010509, 0.024630, 0.044039];

t_kvot = t_M./t_S;
t_av = sum(t_kvot)/length(t_kvot);

Error_kvot = E_m./E_s;
E_av = sum(Error_kvot)/length(Error_kvot);
