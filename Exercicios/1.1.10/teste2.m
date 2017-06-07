1;

n = 200;
for jay = 1:4
    if jay > 1
    %    oldtime = time;
    endif
    A = randn(n);
    x = randn(n, 1);
    b = zeros(n, 1);
    %t = cputime;
    for j = 1:n
        for i = 1:n
            b(i) += A(i, j) * x(j);
        endfor
    endfor
    msize = n
    %time = cputime - t;
    if jay > 1
    %    ratio = time/oldtime
    endif
    n *= 2;
endfor
