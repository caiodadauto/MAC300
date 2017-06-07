1;

n = 1000;
for jay = 1:4
    beginjob = cputime;
    A = randn(n);
    x = randn(n, 1);
    t = cputime;
    b = A * x;
    endjob = cputime;
    printf('Tempo gasto -> %f para n -> %d\n', endjob - beginjob, n);
    n *= 2;
endfor
