#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/***
 * Metodos Numericos em Algebra linear MAC0300  EP2
 *
 * Nome: Caio Vinicius Dadauto
 * Nusp: 7994808
 *
 * Compilacao: gcc gerador.c -lm
***/

/***
 * Distribuicao uniforme em [-1, 1]
***/
double uniform() {
    return ((double)rand()/RAND_MAX * 2 - 1);
}

/***
 * Redimensiona os vetore A e JA
***/
int resize(double **A, int **JA, int *n) {
    int
        i,
        *new_JA;
    double
        *new_A;

    new_A  = (double *)malloc(2 * *n * sizeof(double));
    new_JA = (int *   )malloc(2 * *n * sizeof(int   ));

    if((new_A == NULL) || (new_JA == NULL))
        return -1;
    
    for(i = 0; i < *n; i++) {
        new_A[i]  = A[0][i];
        new_JA[i] = JA[0][i];
    }

    free(*A);
    free(*JA);
    *A  = new_A;
    *JA = new_JA;
    *n  = 2 * *n;

    return 0;
}

/***
 * Gerador de matriz esparsa simétrica e definida positiva no formato Yale
***/
int sparseMatrix(int n, double tau, int **IA, int **JA, double **A) {
    int
        i,
        j,
        m,
        p,
        q,
        k,
        n_nonzero = 0;
    double
        a;

    srand(time(NULL));
    m  = (int)n * n * (tau + 0.003);
    if(m < n)
        m = n;

    *A  = (double *)malloc(m       * sizeof(double));
    *JA = (int *   )malloc(m       * sizeof(int   ));
    *IA = (int *   )malloc((n + 1) * sizeof(int   ));

    if((*A == NULL) || (*JA == NULL) || (*IA == NULL))
        return -1;

    IA[0][0] = n_nonzero;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if(i == j) {
                JA[0][n_nonzero] = j;
                A[0][n_nonzero]  = 1;
                n_nonzero++;
            }
            else if(j > i) {
                a = uniform();
                if(fabs(a) <= tau) {
                    A[0][n_nonzero]  = a;
                    JA[0][n_nonzero] = j;
                    ++n_nonzero;
                }
            }
            else if(j < i) {
                //Sempre q > p
                q    = IA[0][j + 1];
                p    = IA[0][j];
                for(k = p; k < q; k++) {
                    if(JA[0][k] > i)
                        break;
                    if(JA[0][k] == i) {
                        A[0][n_nonzero]  = A[0][k];
                        JA[0][n_nonzero] = j;
                        ++n_nonzero;
                    }
                }
            }
            if(n_nonzero >= m)
                resize(A, JA, &m);
        }
        IA[0][i + 1] = n_nonzero;
    }
    
    return 0;
}

/***
 * Gera o vetor b a partir de valores aleatórios em [-1, 1] 
***/
int bVector(int n, double **b) {
    int
        i;

    srand(time(NULL));
    *b = (double *)malloc(n * sizeof(double));
    if(*b == NULL)
        return -1;

    for(i = 0; i < n; i++)
        b[0][i] = uniform();

    return 0;
}

/***
 * Imprime matriz
***/
void printMatrix(int n, int *IA, int *JA, double *A){
    int
        i,
        j,
        q,
        k,
        pass;

    printf("Formato Yale para matrizes esparsas:\n");
    for(k = 0; k < IA[n]; k++)
        printf("%g  ", A[k]);
    printf("\n");
    for(k = 0; k < IA[n]; k++)
        printf("%d  ", JA[k]);
    printf("\n");
    for(k = 0; k <= n; k++)
        printf("%d  ", IA[k]);
    printf("\n");
    printf("----------------------------------------------\n");

    printf("Matrizes esparsa:\n");
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            pass = 0;
            for(q = IA[i]; q < IA[i + 1]; q++) {
                if(j == JA[q]) {
                    printf("%.3g\t", A[q]);
                    pass = 1;
                }
            }
            if(!pass)
                printf("0\t");
        }
        printf("\n");
    }
}
