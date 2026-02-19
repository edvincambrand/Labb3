




Nvals = [4000, 8000, 16000, 32000, 64000, 128000];
T = 100;
a = 0.5;

for i = 1:length(Nvals)
    N = Nvals(i);
    tic
    [y, H] = mittpunkt(a, N, T);
    toc
    disp(N)
end
