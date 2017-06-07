#include<stdio.h>
#include<stdlib.h>
#include<time.h>

int main() {
    clock_t
        begin,
        end;
    double
        **A,
        *b,
        *x;
    int
        i,
        j,
        k,
        n = 1000;

    srand(time(NULL));
    for(k = 0; k < 4; k++) {
        begin = clock();

        A = (double **)malloc(n * sizeof(double));
        for(i = 0; i < n; i++)
            A[i] = (double *)malloc(n * sizeof(double));

        x = (double *)malloc(n * sizeof(double));
        b = (double *)malloc(n * sizeof(double));

        for(i = 0; i < n; i++) {
            b[i] = 0;
            x[i] = (double)rand()/RAND_MAX;
            for(j = 0; j < n; j++)
                A[i][j] = (double)rand()/RAND_MAX;
        }
        
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++)
                b[i] += A[i][j] * x[j];
        }
        printf("Tempo gasto -> %f com n -> %d\n", (double)(end - begin)/CLOCKS_PER_SEC, n);

        for(i = 0; i < n; i++)
            free(A[i]);
        free(A);
        free(x);
        free(b);

        n *= 2;
    }
}
