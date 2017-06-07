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
 * Compilacao: gcc primeiraParte.c -lm
***/

// Esturtura utilizada para armazenar o sistema linear Ax = b
struct _linearsystem {
    double
        **A,
         *b;
    int
        n;
};

// Simplifica a notacao
typedef struct _linearsystem linearsystem;

// Prototipo de funcoes
int cholcol (int n, double **A);
int cholrow (int n, double **A);
int forwcol (int n, double **A, double *b);
int forwrow (int n, double **A, double *b);
int backcol (int n, double **A, double *b, int trans);
int backrow (int n, double **A, double *b, int trans);
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

    // Aloca a matrix nxn A e o vetor n b
    linear->A = (double **)malloc(linear->n * sizeof(double *));
    for(k = 0; k < linear->n; k++)
        linear->A[k] = (double *)malloc(linear->n * sizeof(double));

    linear->b = (double *)malloc(linear->n * sizeof(double));

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
 * Determina o fator de Cholesky de A e o armazena na propria A.
 * O algoritimo esta´ orientado a coluna.
 ***/
int cholcol(int n, double **A) {
    int
        i,
        j,
        k;

    // Inicia a decomposicao de Cholesky
    for(k = 0; k < n; k++) {
        if(A[k][k] <= 0)
            return -1;
        A[k][k] = sqrt(A[k][k]);
        for(i = k + 1; i < n; i++)
            A[i][k] /= A[k][k];
        for(j = k + 1; j < n; j++) {
            for(i = k + 1; i < n; i++)
                A[i][j] -= A[i][k] * A[j][k];
        }
    }
    
    return 0;
}

/***
 * Resolve um sistema Ax = b com A triangular inferior e armarena
 * a solucao x em b. O algoritimo esta´ orientado a coluna.
 ***/
int forwcol(int n, double **A, double *b) {
    int
        i,
        j;

    // Inicia a resolucao do sistema triangular inferior orientado por coluna
    for(j = 0; j < n; j++) {
        if(A[j][j] == 0)
            return -1;
        b[j] /= A[j][j];
        for(i = j + 1; i < n; i++)
            b[i] -= A[i][j] * b[j]; 
    }

    return 0;
}


/***
 * Resolve um sistema Ax = b com A triangular superior ou com A^t triangular superior
 * e armarena a solucao x em b. O algoritimo esta´ orientado a coluna.
 ***/
int backcol(int n, double **A, double *b, int trans) {
    int
        i,
        j;

    // Inicia a resolucao do sistema triangular superior orientado por coluna para ambos os casos Ax = b e A^tx = b
    if(trans == 0) {
        for(j = n - 1; j >= 0; j--) {
            if(A[j][j] == 0)
                return -1;
            b[j] /= A[j][j];
            for(i = j - 1; i >= 0; i--)
                b[i] -= A[i][j] * b[j];
        }
    }
    else {
        for(j = n - 1; j >= 0; j--) {
            for(i = j + 1; i < n; i++)
                b[j] -= A[i][j] * b[i];
            if(A[j][j] == 0)
                return -1;
            b[j] /= A[j][j];
        }
    }
    
    return 0;
}

/***
 * Determina o fator de Cholesky de A e o armazena na propria A.
 * O algoritimo esta´ orientado a linha.
 ***/
int cholrow(int n, double **A) {
    int
        i,
        j,
        k;

    // Inicia a decomposicao de Cholesky
    for(k = 0; k < n; k++) {
        if(A[k][k] <= 0)
            return -1;
        A[k][k] = sqrt(A[k][k]);
        for(i = k + 1; i < n; i++)
            A[i][k] /= A[k][k];
        for(i = k + 1; i < n; i++) {
            for(j = k + 1; j < n; j++)
                A[i][j] -= A[i][k] * A[j][k];
        }
    }
    
    return 0;
}

/***
 * Resolve um sistema Ax = b com A triangular inferior e armarena
 * a solucao x em b. O algoritimo esta´ orientado a linha.
 ***/
int forwrow(int n, double **A, double *b) {
    int
        i,
        j;
    
    // Inicia a resolucao do sistema triangular inferior orientado por linha
    for(i = 0; i < n; i++) {
        for(j = 0; j <= i - 1; j++)
            b[i] -= A[i][j] * b[j];
        if(A[i][i] == 0)
            return -1;
        b[i] /= A[i][i];
    }

    return 0;
}


/***
 * Resolve um sistema Ax = b com A triangular superior ou com A^t triangular superior
 * e armarena a solucao x em b. O algoritimo esta´ orientado a linha.
 ***/
int backrow(int n, double **A, double *b, int trans) {
    int
        i,
        j;

    // Inicia a resolucao do sistema triangular superior orientado por linha para ambos os casos Ax = b e A^tx = b
    if(trans == 0) {
        for(i = n - 1; i >= 0; i--) {
            for(j = i + 1; j < n; j++)
                b[i] -= A[i][j] * b[j];
            if(A[i][i] == 0)
                return -1;
            b[i] /= A[i][i];
        }
    }
    else {
        for(i = n - 1; i >= 0; i--) {
            if(A[i][i] == 0)
                return -1;
            b[i] /= A[i][i];
            for(j = i - 1; j >= 0; j--)
                b[j] -= A[i][j] * b[i];
        }
    }
    
    return 0;
}

int main() {
    double
        cholcol_t,
        forwcol_t,
        backcol_t,
        cholrow_t,
        forwrow_t,
        backrow_t;
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
    if(cholcol(linear.n, linear.A) == -1){
        printf("Matriz A nao e´ definida positiva!\n");
        return -1;
    }
    end = clock();

    // Determina o tempo gasto para gerar o fator de Cholesky orientado a coluna
    cholcol_t = (double)(end - begin)/CLOCKS_PER_SEC;

    begin = clock();
    if(forwcol(linear.n, linear.A, linear.b) == -1){
        printf("Matriz triangular e´ singular!\n");
        return -1;
    }
    end = clock();

    // Determina o tempo gasto para resolver o sistema triangular inferior orientado a coluna
    forwcol_t = (double)(end - begin)/CLOCKS_PER_SEC;

    begin = clock();
    if(backcol(linear.n, linear.A, linear.b, 1) == -1){
        printf("Matriz triangular e´ singular!\n");
        return -1;
    }
    end = clock();

    // Determina o tempo gasto para resolver o sistema triangular superior orientado a coluna
    backcol_t = (double)(end - begin)/CLOCKS_PER_SEC;

    printf("\nTempo para orientacao a coluna aplicada ao arquivo %s:\n", filename);
    printf("        Cholesky        Forward        Backward\n");
    printf("        %5f        %5f       %5f\n", cholcol_t, forwcol_t, backcol_t);

    // Releitura do arquivo de entrada para restaurar o sistema linear
    readfile(&linear, filename);

    /***
    * Resolvendo o sistema linear pela orientacao por linha
    ***/
    begin = clock();
    if(cholrow(linear.n, linear.A) == -1){
        printf("Matriz A nao e´ definida positiva!\n");
        return -1;
    }
    end = clock();

    // Determina o tempo gasto para gerar o fator de Cholesky orientado a linha
    cholrow_t = (double)(end - begin)/CLOCKS_PER_SEC;

    begin = clock();
    if(forwrow(linear.n, linear.A, linear.b) == -1){
        printf("Matriz triangular e´ singular!\n");
        return -1;
    }
    end = clock();

    // Determina o tempo gasto para resolver o sistema triangular inferior orientado a linha
    forwrow_t = (double)(end - begin)/CLOCKS_PER_SEC;

    begin = clock();
    if(backrow(linear.n, linear.A, linear.b, 1) == -1){
        printf("Matriz triangular e´ singular!\n");
        return -1;
    }
    end = clock();

    // Determina o tempo gasto para resolver o sistema triangular superior orientado a linha
    backrow_t = (double)(end - begin)/CLOCKS_PER_SEC;

    printf("\nTempos para Orientacao a linha aplicada ao arquivo %s:\n", filename);
    printf("        Cholesky        Forward        Backward\n");
    printf("        %5f        %5f       %5f\n", cholrow_t, forwrow_t, backrow_t);

    printf("\nRazao entre os tempos orientado a linha e orientado a coluna:\n");
    printf("\tcholrow/cholcol -> %g\n\tforwrow/forwcol -> %g\n\tbackrow/backcol -> %g\n",
            cholrow_t/cholcol_t, forwrow_t/forwcol_t, backrow_t/backcol_t);

    return 0;
}
