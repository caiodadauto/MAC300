#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"gerador.c"

/***
 * Metodos Numericos em Algebra linear MAC0300  EP2
 *
 * Nome: Caio Vinicius Dadauto
 * Nusp: 7994808
 *
 * Dependencia: gerador.c
 *
 * Compilacao: gcc gradienteConjugado.c -lm
***/

/***
 * Aloca vetor na memoria e verifica se houve erro
***/
double *allocateVector(int n) {
    double
        *v;

    v = (double *)malloc(n * sizeof(double));
    if(v == NULL) {
        printf("Memoria insuficiente\n");
        exit(1);
    }
    
    return v;
}

/***
 * Produto de um vetor por um escalar
***/
double *dotScalarVector(int n, double s, double *u) {
    int
        i;
    double
        *v;

    v = allocateVector(n);
    for(i = 0; i < n; i++)
        v[i] = s * u[i];

    return v;
}

/***
 * Soma dois vetores
***/
double *sumVectors(int n, double *u, double *v) {
    int
        i;
    double
        *w;

    w = allocateVector(n);
    for(i = 0; i < n; i++)
        w[i] = u[i] + v[i];

    return w;
}

/***
 * Subtrai dois vetores
***/
double *subVectors(int n, double *u, double *v) {
    int
        i;
    double
        *w;

    w = allocateVector(n);
    for(i = 0; i < n; i++)
        w[i] = u[i] - v[i];

    return w;
}

/***
 * Faz o produto escalar entre dois vetores
***/
double dot(int n, double *u, double *v) {
    int
        i;
    double
        sum = 0;

    for(i = 0; i < n; i++)
        sum += u[i] * v[i];

    return sum;
}

/***
 * Faz produto de uma matriz no armazenada no formato yale por um vetor
***/
double *dotMatrizSparseVector(int n, double *v, double *A, int *JA, int *IA) {
    int
        i,
        k,
        first_nonzero_row,
        first_nonzero_next_row;
    double
        *u;

    u = allocateVector(n);
    for(i = 0; i < n; i++)
        u[i] = 0;

    for(i = 0; i < n; i++) {
        first_nonzero_row       = IA[i];
        first_nonzero_next_row  = IA[i + 1];

        for(k = first_nonzero_row; k < first_nonzero_next_row; k++)
            u[JA[k]] += A[k] * v[JA[k]];   
    }

    return u;
}

int main(){
    int
        i,
        n,
        *IA,
        *JA;
    double
        tau,
        alpha,
        beta,
        *A,
        *b,
        *x,
        *d,
        *g,
        *x_old,
        *g_old,
        *d_old,
        *A_d,
        *beta_d,
        *alpha_d,
        *alpha_A_d;

    printf("Entre com o tamanho da matriz esparsa: ");
    scanf("%d", &n);
    printf("Entre com o parametro tau (0, 1) para gerar a matriz esparsa: ");
    scanf("%lf", &tau);

    if(sparseMatrix(n, tau, &IA, &JA, &A) == -1) {
        printf("Nao foi possivel alocar a matriz A\n");
        return -1;
    }
    if(bVector(n, &b) == -1) {
        printf("Nao foi possivel alocar o vetor b\n");
        return -1;
    }

    //Aloca x, g, e d, e os inicializa com os valores iniciais
    x = allocateVector(n);
    for(i = 0; i < n; i++)
        x[i] = 0;

    g = allocateVector(n);
    for(i = 0; i < n; i++)
        g[i] = b[i];
    
    d = allocateVector(n);
    for(i = 0; i < n; i++)
        d[i] = b[i];

    //Inicia o metodo de gradientes conjugados
    for(i = 0; i < 20; i++) {
        x_old    = x;
        g_old    = g;
        d_old    = d;

        A_d       = dotMatrizSparseVector(n, d, A, JA, IA);
        alpha     = dot(n, g, g) / dot(n, d, A_d);

        alpha_d   = dotScalarVector(n, alpha, d);
        x         = sumVectors(n, x, alpha_d);

        alpha_A_d = dotScalarVector(n, alpha, A_d);
        g         = subVectors(n, g, alpha_A_d);

        beta      = dot(n, g, g) / dot(n, g_old, g_old);

        beta_d    = dotScalarVector(n, beta, d);
        d         = sumVectors(n, g, beta_d);
        
        printf("Residuo(%d) = %g\n", i, sqrt(dot(n, g, g)));

        free(x_old    );
        free(g_old    );
        free(d_old    );
        free(A_d      );
        free(alpha_d  );
        free(alpha_A_d);
        free(beta_d   );
    }

    free(x );
    free(g );
    free(d );
    free(b );
    free(A );
    free(JA);
    free(IA);
    
    return 0;
}
