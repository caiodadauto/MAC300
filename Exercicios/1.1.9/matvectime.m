1;

n = 800;
for jay = 1:4
    if jay > 1
        oldtime = time;
    endif
    A = randn(n);
    x = randn(n, 1);
    t = cputime;
    b = A * x;
    msize = n
    time = cputime - t;
    if jay > 1
        ratio = time/oldtime
    endif
    n *= 2;
endfor
