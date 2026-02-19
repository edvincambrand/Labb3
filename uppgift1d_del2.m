%% Mittpunktsmetoden och symplektisk euler

% Vi tillkallar koden genom mittpunkt.m och symplektiska.m

clear;

% Parametrar
a = 0.5;
N = 16000;
T = 100;

% Ta fram lösningar + energier
[y1, H1] = mittpunkt(a, N, T);
[y2, H2] = symplektiska(a, N, T);

% Mittpunkt:
q1_m = y1(:,1);
q2_m = y1(:,2);

% Symplektiska:
q1_s = y2(:,1);
q2_s = y2(:,2);

% Plotta q1, q2 för mittpunkt:
subplot(1,2,1);
plot(q1_m, q2_m); grid on;
xlabel("q1");
ylabel("q2");
title("Implicita mittpunktsmetoden")
xlim([-1.5, .5])
ylim([-1, 1])

% Plotta q1/q2 för symplektiska:
subplot(1,2,2);
plot(q1_s, q2_s); grid on;
xlabel("q1")
ylabel("q2")
title("Eulers symplektiska metod")
xlim([-1.5, .5])
ylim([-1, 1])

% Plotta energier:
figure;
plot(H1); hold on; plot(H2); grid on;
legend("H : Mittpunktsmetoden", "H : Eulers symplektiska metod")
ylabel("H : Energi")
xlabel("Iteration, i")
xlim([0,N])
ylim([-0.55, -0.45])

% Skriv ut energiernas intervallgränser:
fprintf("Range (mittpunkt): [%.8f, %.8f] | Range (symplektisk): [%.8f, %.8f] \n", ...
    min(H1),max(H1),min(H2),max(H2))

