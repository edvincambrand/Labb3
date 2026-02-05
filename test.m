clear;

a = 0.5;
N = 4000;
T = 100;

[y, H] = mittpunkt(a, N, T);

q1_i = y(:,1);
q2_i = y(:,2);

[v, E] = symplektiska(a, N, T);
q1_s = v(:,1);
q2_s = v(:,2);

subplot(1,2,1);
plot(q1_i,q2_i); hold on;
plot(q1_s,q2_s);

subplot(1,2,2);
plot(H); hold on; plot(E)

