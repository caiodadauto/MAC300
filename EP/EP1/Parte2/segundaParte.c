#include <stdlib.h>
#include <stdio.h>
#include <math.h> 
#include <time.h>

/***
 * Metodos Numericos em Algebra linear MAC0300  EP1
 *
 * Nome: Caio Vinicius Dadauto
 * Nusp: 7994808
 *
 * Compilacao: gcc segundaParte.c -lm
***/

// Esturtura utilizada para armazenar o sistema linear PAx = Pb
struct _linearsystem {
    double
        **A,
         *b;
    int
        *p,
         n;
};

// Simplifica a notacao
typedef struct _linearsystem linearsystem;

// Prototipo de funcoes
int lucol (int n, double **A, int *p);
int lurow (int n, double **A, int *p);
int sscol (int n, double **A, int *p, double *b);
int ssrow (int n, double **A, int *p, double *b);
int readfile(linearsystem *linear, char filename[]);

/***
 * Le de um arquivo um sistema linear Ax = b,
 * onde as n^2 primeiras linhas sao A e as n restantes sao b.
 ***/
int readfile(linearsystem *linear, char filename[]) {
    double
        a;
    int
        i,
        j,
        k;
    FILE
        *input;

    // Abre o arquivo com nome especificado
    input = fopen(filename, "r");
    if(input == NULL) {
        printf("\tErro na leitura do arquivo\n");
        return -1;
    }
    fscanf(input, "%d", &linear->n);

    // Aloca a matrix nxn A e o vetor b e p
    linear->A = (double **)malloc(linear->n * sizeof(double *));
    for(k = 0; k < linear->n; k++)
        linear->A[k] = (double *)malloc(linear->n * sizeof(double));

    linear->b = (double *)malloc(linear->n * sizeof(double));
    linear->p = (int *)   malloc(linear->n * sizeof(int));

    // Le o arquivo e inicializa A e b
    for(k = 0; k < linear->n * linear->n; k++) {
        fscanf(input, "%d%d%lf", &i, &j, &a);
        linear->A[i][j] = a;
    }
    
    for(k = 0; k < linear->n; k++) {
        fscanf(input, "%d %lf", &i, &a);
        linear->b[i] = a;
    }

    fclose(input);
    return 0;
}

/***
 * Determina a decomposicao de PA em LU e
 * armazena LU em A e P em p.
 * A implementacao e´ feita orientada a coluna.
 ***/
int lucol(int n, double **A, int *p) {
    double
        tmp;
    int
        i,
        j,
        k,
        i_max;

    // Inicia a decomposicao LU
    for(k = 0; k < n; k++) {
        i_max = k;
        for(i = k + 1; i < n; i++)
            if(fabs(A[i][k]) > fabs(A[i_max][k]))
                i_max = i;

        p[k] = i_max;
        if(i_max != k) {
            for(j = 0; j < n; j++) {
                tmp         = A[k][j];
                A[k][j]     = A[i_max][j];
                A[i_max][j] = tmp;
            }
        }
        if(A[k][k] == 0)
            return -1;

        for(i = k + 1; i < n; i++)
            A[i][k] /= A[k][k];

        for(j = k + 1; j < n; j++)
            for(i = k + 1; i < n; i++)
                A[i][j] -= A[i][k] * A[k][j];
    }
    
    return 0;
}

/***
 * Determina a solucao de LUx = Pb e a
 * armazena em b.
 * A implementacao e´ feita orientada a coluna.
 ***/
int sscol(int n, double **A, int *p, double *b) {
    int
        i,
        j;
    double
        tmp;

    // Determina Pb
    for(i = 0; i < n - 1; i++) {
        tmp     = b[i];
        b[i]    = b[p[i]];
        b[p[i]] = tmp;
    }

    // Determina o subsistema Ly = Pb, onde a diagonal de L e´ unitaria
    for(j = 0; j < n; j++)
        for(i = j + 1; i < n; i++)
            b[i] -= A[i][j] * b[j];

    // Determina o subsistema Ux = y 
    for(j = n - 1; j >= 0; j--) {
        if(A[j][j] == 0)
            return -1;
        b[j] /= A[j][j];
        for(i = j - 1; i >= 0; i--)
            b[i] -= A[i][j] * b[j];
    }

    return 0;
}

/***
 * Determina a decomposicao de PA em LU e
 * armazena LU em A e P em p.
 * A implementacao e´ feita orientada a linha.
 ***/
int lurow(int n, double **A, int *p) {
    int
        i,
        j,
        k,
        i_max;
    double
        tmp;

    // Inicia a decomposicao LU
    for(k = 0; k < n; k++) {
        i_max = k;
        for(i = k + 1; i < n; i++)
            if(fabs(A[i][k]) > fabs(A[i_max][k]))
                i_max = i;
        p[k] = i_max;
        if(i_max != k) {
            for(j = 0; j < n; j++) {
                tmp         = A[k][j];
                A[k][j]     = A[i_max][j];
                A[i_max][j] = tmp;
            }
        }
        if(A[k][k] == 0)
            return -1;

        for(i = k + 1; i < n; i++) {
            A[i][k] /= A[k][k];
            for(j = k + 1; j < n; j++)
                A[i][j] -= A[i][k] * A[k][j];
        }
    }
    
    return 0;
}

/***
 * Determina a solucao de LUx = Pb e a
 * armazena em b.
 * A implementacao e´ feita orientada a linha.
 ***/
int ssrow(int n, double **A, int *p, double *b) {
    int
        i,
        j;
    double
        tmp;

    // Determina Pb
    for(i = 0; i < n - 1; i++) {
        tmp     = b[i];
        b[i]    = b[p[i]];
        b[p[i]] = tmp;
    }

    // Determina o subsistema Ly = Pb, onde a diagonal de L e´ unitaria
    for(i = 0; i < n; i++)
        for(j = 0; j <= i - 1; j++)
            b[i] -= A[i][j] * b[j];

    // Determina o subsistema Ux = y 
    for(i = n - 1; i >= 0; i--) {
        for(j = i + 1; j < n; j++)
            b[i] -= A[i][j] * b[j];
        if(A[i][i] == 0)
            return -1;
        b[i] /= A[i][i];
    }

    return 0;
}

int main() {
    double
        lucol_t,
        sscol_t,
        lurow_t,
        ssrow_t;
    clock_t
        begin,
        end;
    char
        filename[100];
    linearsystem
        linear;

    do{
        printf("\nEntre com o nome do arquivo que contem o sistema linear a ser resolvido: ");
        scanf(" %[^\n]", filename);
    } while(readfile(&linear, filename) == -1);

    /***
    * Resolvendo o sistema linear pela orientacao por coluna
    ***/
    begin = clock();
    if(lucol(linear.n, linear.A, linear.p) == -1){
        printf("Matriz A e´ singular!\n");
        return -1;
    }
    end = clock();

    // Determina o tempo gasto para decompor PA em LU orientado a coluna
    lucol_t = (double)(end - begin)/CLOCKS_PER_SEC;

    begin = clock();
    if(sscol(linear.n, linear.A, linear.p, linear.b) == -1){
        printf("Matriz triangular superior e´ singular!\n");
        return -1;
    }
    end = clock();

    // Determina o tempo gasto para resolver o sistema LUx = Pb
    sscol_t = (double)(end - begin)/CLOCKS_PER_SEC;

    printf("\nTempo para orientacao a coluna aplicada ao arquivo %s:\n", filename);
    printf("        PA = LU        LUx = Pb\n");
    printf("        %5f       %5f\n", lucol_t, sscol_t);
   
    // Releitura do arquivo de entrada para restaurar o sistema linear
    readfile(&linear, filename);

    /***
    * Resolvendo o sistema linear pela orientacao por linha
    ***/
    begin = clock();
    if(lurow(linear.n, linear.A, linear.p) == -1){
        printf("Matriz A e´ singular!\n");
        return -1;
    }
    end = clock();

    // Determina o tempo gasto para gerar a decomposicao LU orientado a linha
    lurow_t = (double)(end - begin)/CLOCKS_PER_SEC;

    begin = clock();
    if(ssrow(linear.n, linear.A, linear.p, linear.b) == -1){
        printf("Matriz triangular superior e´ singular!\n");
        return -1;
    }
    end = clock();

    // Determina o tempo gasto para resolver o sistema LUx = Pb
    ssrow_t = (double)(end - begin)/CLOCKS_PER_SEC;

    printf("\nTempos para Orientacao a linha aplicada ao arquivo %s:\n", filename);
    printf("        PA = LU        LUx = Pb\n");
    printf("        %5f       %5f\n", lurow_t, ssrow_t);

    printf("\nRazao entre os tempos orientado a linha e orientado a coluna:\n");
    printf("\tlurow/lucol -> %g\n\tssrow/sscol -> %g\n", lurow_t/lucol_t, ssrow_t/sscol_t);

    return 0;
}
