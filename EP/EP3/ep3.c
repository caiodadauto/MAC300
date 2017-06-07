#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<float.h>

/***
 * Metodos Numericos em Algebra linear MAC0300  EP3
 *
 * Nome: Caio Vinicius Dadauto
 * Nusp: 7994808
 *
 * Compilacao: gcc ep3.c -lm
***/

struct _linearsystem {
    double
        **A,
         *b,
         factor,
         middle;
    int
        n,
        m;
};

typedef struct _linearsystem linearsystem;

/*-----------------------------------Funcoes Auxiliares------------------------------------*/
int     *allocateIntVector(int n);
double  *allocateVector(int n); 
double **allocateMatrix(int row, int col);
/*-----------------------------------Funcoes para Teste------------------------------------*/
int      readFile(linearsystem *linear, char filename[], int opt);
void     systemParse(FILE *input, linearsystem *linear);              
void     pointsParse(FILE *input, linearsystem *linear);
/*----------------------Funcoes que Implementam decomposicao QR----------------------------*/
void     columnNorms(int n, int m, double *norms, double **A, double *b);
int      largerNorm(int n, int m, int step, int *P, double *norms, double **A);
void     reflectorGenerator(int n,int step, double norm, double *gammas, double **A);
void     reflectorProductMatrix(int n, int m, int step, double gamma, double **A);
/*----------------------Funcoes que Implementam MMQ com refletores-------------------------*/
void     reflectorProductVector(int n, int step, double gamma, double *b, double **A);
void     backrow(int rank, int m, int *P, double *b, double **A);
/*-----------------------------------------------------------------------------------------*/


/***
 * Aloca vetor de double na memoria
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
 * Aloca vetor de int na memoria
***/
int *allocateIntVector(int n) {
    int
        *v;

    v = (int *)malloc(n * sizeof(int));
    if(v == NULL) {
        printf("Memoria insuficiente\n");
        exit(1);
    }
    
    return v;
}

/***
 * Aloca Matriz na memoria
***/
double **allocateMatrix(int row, int col) {
    int
        i;
    double
        **A;
    
    A  = (double **)malloc(row * sizeof(double *));
    if(A == NULL) {
        printf("Memoria insuficiente\n");
        exit(1);
    }

    for(i = 0; i < row; i++) {
        A[i] = (double *)malloc(col * sizeof(double));
        if(A[i] == NULL) {
            printf("Memoria insuficiente\n");
            exit(1);
        }
    }
    
    return A;
}

/***
 * Analisa o arquivo de entrada mediante a opcao escolhida
 * pelo usuario.
***/
int readFile(linearsystem *linear, char filename[], int opt) {
    FILE
        *input;

    input = fopen(filename, "r");
    if(input == NULL) {
        printf("\tErro na leitura do arquivo\n");
        return -1;
    }

    if(opt == 1) {
        systemParse(input, linear);
    }
    else if(opt == 2) {
        pointsParse(input, linear);
    }

    fclose(input);
    return 0;
}

/***
 * Analisa o arquivo de entrada considerando que foi
 * passado uma matriz A e um vetor b como entrada.
***/
void systemParse(FILE *input, linearsystem *linear) {
    double
        a;
    int
        i,
        j,
        k;

    fscanf(input, "%d%d", &linear->n, &linear->m);

    linear->A = allocateMatrix(linear->n, linear->m);
    linear->b = allocateVector(linear->n);

    for(k = 0; k < linear->n * linear->n; k++) {
        fscanf(input, "%d%d%lf", &i, &j, &a);
        linear->A[i][j] = a;
    }
    
    for(k = 0; k < linear->n; k++) {
        fscanf(input, "%d %lf", &i, &a);
        linear->b[i] = a;
    }
}

/***
 * Analisa o arquivo de entrada considerando que foi
 * passado um conjunto de pontos (x, y) e o grau do
 * polinomio a ser ajustado aos pontos. A base de
 * polinomios escolhida, é baseada nos valores maximo
 * e minimo das entradas em x de forma a garantir que
 * a matriz A seja bem condicionada.
***/
void pointsParse(FILE *input, linearsystem *linear) {
    int
        i,
        j,
        n,
        degree;
    double
        s,
        *x,
        *y,
        max,
        min,
        middle,
        factor;

    fscanf(input, "%d", &n);
    fscanf(input, "%d", &degree);
    
    x = allocateVector(n);
    y = allocateVector(n);

    linear->n = n;
    linear->m = degree + 1;
    linear->A = allocateMatrix(linear->n, linear->m);
    linear->b = allocateVector(linear->n);

    for(i = 0; i < n; i++)
        fscanf(input, "%lf%lf", &x[i], &y[i]);
    
    max = DBL_MIN;
    min = DBL_MAX;
    for(i = 0; i < n; i++) {
        if(x[i] > max)
            max = x[i];
        if(x[i] < min)
            min = x[i];
    }

    middle = (max + min)/ 2;
    factor = fabs(max - middle);
    for(i = 0; i < linear->n; i++) {
        s = (x[i] - middle)/factor;
        for(j = 0; j < linear->m; j++)
            linear->A[i][j] = pow(s, j);
    }

    for(i = 0; i < linear->n; i++)
        linear->b[i] = y[i];

    linear->middle = middle;
    linear->factor = factor;

    free(x);
    free(y);
}

/***
 * Determina o vetor auxiliar que armazena as normas das
 * colunas de A. Alem disso, normaliza a matriz A e o 
 * vetor b a partir da razao de todas as entradas de A 
 * e de b pela a entrada de maior modulo de A.
***/
void columnNorms(int n, int m, double *norms, double **A, double *b) {
    int
        i,
        j;
    double
        mod,
        max;

    max = 0;
    for(i = 0; i < n; i++) {
        for(j = 0; j < m; j++) {
            mod = fabs(A[i][j]);
            if(mod > max)
                max = mod;
        }
    }

    if(max == 0) {
        printf("Matriz nula.\n");
        exit(1);
    }

    for(i = 0; i < n; i++) {
        b[i] /= max;
        for(j = 0; j < m; j++)
            A[i][j] /= max; 
    }

    for(i = 0; i < n; i++) {
        for(j = 0; j < m; j++)
            norms[j] += A[i][j] * A[i][j];
    }
}

/***
 * Determina a coluna com maior norma em uma matriz 
 * n - step x m - step e permuta a coluna do passo 
 * atual com a de maior norma, alem de permutar as 
 * normas e armazenar a permutacao no vetor P. 
 * As normas sao corrigidas para cada passo.
 * Caso a maior norma seja 0, retorna -1.
***/
int largerNorm(int n, int m, int step, int *P, double *norms, double **A) {
    int
        j,
        row,
        j_max;
    double
        tmp;
    
    if(step == m)
        return -1;

    if(step > 0) {
        row = step - 1;
        for(j = step; j < m; j++)
            norms[j] -= A[row][j] * A[row][j];
    }

    j_max = step;
    for(j = step; j < m - 1; j++) {
        if(norms[j] < norms[j + 1])
            j_max = j + 1;
    }

    if(norms[j_max] < DBL_EPSILON)
        return -1;

    if(j_max != step) {
        for(j = 0; j < n; j++) {
            tmp = A[j][step];
            A[j][step]  = A[j][j_max];
            A[j][j_max] = tmp;
        }

        tmp = norms[step];
        norms[step]  = norms[j_max];
        norms[j_max] = tmp;
    }

    P[step] = j_max;
    return 0;
}

/***
 * Determina o vetor u associado a uma dada coluna
 * de A. Armazena u e tau sobre a mesma coluna de A
 * e armazena gamma em um vetor auxiliar. Note que
 * a primeira coordenada de u e´ 1, por isso nao e´
 * armazenada. Alem disso, tau nao precisa ser 
 * renormalizado pois se dividiu todas as entradas
 * de A e de b por max, logo o problema ainda e´ o
 * mesmo.
***/
void reflectorGenerator(int n, int step, double norm, double *gammas, double **A) {
    int
        i;
    double
        tau;
    
    tau = sqrt(norm);

    if(A[step][step] < 0)
        tau *= -1;

    A[step][step] += tau;
    gammas[step]   = A[step][step]/tau;

    for(i = step + 1; i < n; i++)
        A[i][step] /= A[step][step];

    A[step][step] = -1 * tau;
}

/***
 * Determina o produto de uma matriz A (n x m)
 * por uma matriz de reflexao Q = I - gamma * u * u^T.
 * O resultado e´ armazenado sobre a propria A.
 * Esse produto e´ implementado de forma a ser
 * orientado a linhas. Note que a primeira entrada
 * nao nula de u vale 1 e nao esta armazenado em A.
***/
void reflectorProductMatrix(int n, int m, int step, double gamma, double **A) {
    int
        i,
        j;
    double
        *sum;

    sum  = allocateVector(m);

    for(j = step + 1; j < m; j++)
        sum[j] = gamma * A[step][j];

    for(i = step + 1; i < n; i++)
        for(j = step + 1; j < m; j++)
            sum[j] += gamma * A[i][step] * A[i][j];
 
    for(j = step + 1; j < m; j++)
        A[step][j] -= sum[j];

    for(i = step + 1; i < n; i++) 
        for(j = step + 1; j < m; j++)
            A[i][j] -= A[i][step] * sum[j];

    free(sum);
}

/***
 * Determina o produto de uma vetor b
 * por uma matriz de reflexao Q = I - gamma * u * u^T.
 * O resultado e´ armazenado sobre o propria b.
 * Esse produto e´ implementado de forma a ser
 * orientado a linhas. Note que a primeira entrada
 * nao nula de u vale 1 e nao esta armazenado em A.
***/
void reflectorProductVector(int n, int step, double gamma, double *b, double **A) {
    int
        i;
    double
        sum;

    sum = gamma * b[step];

    for(i = step + 1; i < n; i++)
        sum += gamma * A[i][step] * b[i];
 
    b[step] -= sum;

    for(i = step + 1; i < n; i++)
        b[i] -= A[i][step] * sum;
}
/***
 * Supondo R = [R_1 R_2] e y^T = [y_1 y_2], determina o valor
 * de y em,
 *              R_1 y_1 = c_1 - R_2 y_2
 * no caso de posto completo R_2 e y_2 nao existem, pois possuem
 * dimensao 0. No caso do posto ser deficiente, y_2 e´ arbitrario,
 * portanto toma-se y_2 = 0. Assim, basta resolver,
 *              R_1 y_1 = c_1
 * onde y = [y_1 0], com x = P^T y. Em seguida determina o valor de
 * x que minimiza o problema Ax = b.
***/
void backrow(int rank, int m, int *P, double *b, double **A) {
    int
        i,
        j;
    double
        tmp;

    for(i = rank - 1; i >= 0; i--) {
        for(j = i + 1; j < rank; j++)
            b[i] -= A[i][j] * b[j];
        b[i] /= A[i][i];
    }

    for(i = rank; i < m; i++)
        b[i] = 0;
    for(i = 0; i < rank; i++) {
        tmp  = b[i]; 
        b[i] = b[P[i]];
        b[P[i]] = tmp;
    }
}

int main() {
    int
        i,
        *P,
        opt,
        rank;
    double
        *norms,
        *gammas;
    char
        filename[50];
    linearsystem
        linear;

    printf("\tEste programa resolve o problema de minimos quadrados para\n");
    printf("\tpolinomios de grau m - 1 utilizando matrizes de reflexao.\n");
    printf("\tO programa aceita como entrada um arquivo texto, este arquivo\n");
    printf("\tpode conter uma matriz A e um vetor b de um sistema Ax = b, ou\n");
    printf("\tentao, pode conter dados explicitos de  um problema de minimizacao,\n");
    printf("\tou seja, as coordenadas (x, y) de pontos a serem minimizados.\n");
    printf("\tPara isso, entre com ( 1 ) se desejar inserir a matriz A e o vetor b\n");
    printf("\tdiretamente, caso contrario entre com ( 2 ) para especificar os pontos\n");
    printf("\ta serem minimizados.\n");
    do {
        printf("\n\tEntre com a opcao desejada: ");
        scanf("%d", &opt);
    } while(!(opt == 1 || opt == 2));

    do {
        printf("\tEntre com nome do arquivo a ser lido: ");
        scanf(" %[^\n]", filename);
    } while(readFile(&linear, filename, opt) == -1);

    P      = allocateIntVector(linear.m);
    norms  = allocateVector(linear.m);
    gammas = allocateVector(linear.m);

    for(i = 0; i < linear.m; i++)
        norms[i] = 0;

    columnNorms(linear.n, linear.m, norms, linear.A, linear.b );
    for(rank = 0; largerNorm(linear.n, linear.m, rank, P, norms, linear.A) != -1; rank++) {
        reflectorGenerator(linear.n, rank, norms[rank], gammas, linear.A);
        reflectorProductMatrix(linear.n, linear.m, rank, gammas[rank], linear.A);
    }

    for(i = 0; i < rank; i++)
        reflectorProductVector(linear.n, i, gammas[i], linear.b, linear.A);
    
    backrow(rank, linear.m, P, linear.b, linear.A);

    if(rank == linear.m) {
        if(opt == 2) {
            printf("\n\tA possue posto completo, ou seja, posto de A e´ %d.\n", rank);
            printf("\tComo solucao do problema dos quadrados minimos, obteve-se\n");
            printf("\to seguinte polinomio,\n\n\t\t");
            for(i = 0; i < linear.m - 1; i++)
                printf("%g * s^%d + ", linear.b[i], i);
            printf("%g * s^%d\n\n", linear.b[linear.m - 1], linear.m - 1);
            printf("\tonde s = (t - %g)/%g \n", linear.middle, linear.factor);
        }
        else {
            printf("\n\tA possue posto completo, ou seja, posto de A é %d.\n", rank);
            if(rank <= 20) {
                printf("\tA solucao x para o MMQ e´ igual a,\n\n\t[");
                for(i = 0; i < linear.m - 1; i++)
                    printf("%g ", linear.b[i]);
                printf("%g]\n", linear.b[linear.m - 1]);
            }
        }
    }
    else {
        printf("\n\tA não possui posto completo, a saber o posto e´ %d.\n", rank);
        printf("\tComo ha´ infinitas solucoes possiveis ao problema de\n");
        printf("\tquadrados minimos, por praticidade escolhida uma solucao\n");
        printf("\tespecifica que toma y_2 como vetor nulo.\n");
        if(rank <= 20) {
            printf("\tEssa  solucao exemplo para o MMQ e´ igual a,\n\n\t\t[");
            for(i = 0; i < linear.m - 1; i++)
                printf("%g ", linear.b[i]);
            printf("%g]\n", linear.b[linear.m - 1]);
        }
    }

    for(i = 0; i < linear.m; i++)
        free(linear.A[i]);
    free(linear.A);
    free(linear.b);
    free(norms);
    free(gammas);
    free(P);

    return 0;
}
