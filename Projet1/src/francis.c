#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void descente(int *b, int *a, int taille, char op);
void final(int *b, int *a, int taille, char op);
void montee(int *tab, int taille, char op);
int operation(int, int, char);

int elem_neutre(char op)
{
    switch (op)
    {
    case '+':
        return 0;

    case '-':
        return 0;

    case '*':
        return 1;

    case '/':
        return 1;

    case 'm':
        return -1410065407;

    default:
        fprintf(stderr, "Error. Invalid opeator !\nPrograme ended with code error -1\n\a");
        exit(EXIT_FAILURE);
    }
}

// O(1) sans parallèle
int operation(int val1, int val2, char op)
{
    switch (op)
    {
    case '+':
        return val1 + val2;
        break;
    case '-':
        return val1 - val2;
        break;
    case '*':
        return val1 * val2;
        break;
    case '/':
        if (val2 != 0)
            return val1 / val2;
        else
        {
            fprintf(stderr, "Error. Divding by zero !\nPrograme ended with code error -1\n\a");
            exit(EXIT_FAILURE);
        }
        break;
    
    case 'm':
        return val1>val2 ? val1:val2; 
    default:
        fprintf(stderr, "Error. Invalid opeator !\nPrograme ended with code error -1\n\a");
        exit(EXIT_FAILURE);
    }
}
// O(log(N)) en parallel et O(Nlog(N)) en normal
void descente(int *b, int *a, int taille, char op)
{
    int m = log2(taille);

    b[0] = 0;

    for (int l = 1; l < m + 1; l++)
    {
        //printf("l = %d\n", l);
        int maxi = pow(2, l + 1);
#pragma omp parallel for
        for (int j = pow(2, l) - 1; j < maxi - 1; j++)
        {
            //printf("\tj = %d\n",j);
            if ((j + 1) % 2 == 0)
            {
                b[j] = b[(j - 1) / 2];
            }
            else
            {

                b[j] = operation(b[(j - 1) / 2], a[j - 1], op);
            }
        }
    }
}

int isPowerOfTwo(int x)
{
    return (x != 0) && ((x & (x - 1)) == 0);
}

void montee(int *tab, int taille, char op)
{
    int m = log2(taille);
    //printf("m = %d\n",m);
    for (int l = m - 1; l >= 0; l--)
    {
        //printf("l = %d:\n ",l);
        int max = pow(2, l + 1) - 1;
#pragma omp parallel for
        for (int j = pow(2, l) - 1; j < max; j++)
        {
            //printf("\tj = %d\n",j);
            //printf("\ttab[%d], tab[%d] = %d, %d",2*j+1,2*j+2, tab[2*j+1], tab[2*j+2]);
            tab[j] = operation(tab[2 * j + 1], tab[2 * j + 2], op);
            //printf("\ttab[%d] = %d\n", j, tab[j]);
        }
        //printf("\n");
    }
}

int *create_copieTable(int *a, int taille, char op)
{

    int *tab_to_ret = (int *)malloc(sizeof(int *) * (pow(2, log2(taille) + 1) - 1));
#pragma omp parallel for
    for (int i = 0; i < 2 * taille - 1; i++)
        tab_to_ret[i] = elem_neutre(op);
#pragma omp parallel for
    for (int i = 0; i < taille; i++)
    {

        tab_to_ret[(taille - 1) + i] = a[i];
    }

    return tab_to_ret;
}

void final(int *b, int *a, int taille, char op)
{

    int m = log2(taille);
    int maxi = pow(2, m + 1) - 1;
#pragma omp parallel for
    for (int j = pow(2, m) - 1; j < maxi; j++)
    {
        //printf("j = %d\t",j);
        b[j] = operation(b[j], a[j], op);
        //printf("b[%d] = %d\n", j, b[j]);
    }
}

void remplir(int *t1, int *t2, int *t, int n1, int n2)
{
    int index_on_t = 0;
#pragma omp parallel for
    for (int i = 0; i < n1; i++)
    {
        t1[i] = t[index_on_t++];
    }
#pragma omp parallel for
    for (int i = 0; i < n2; i++)
    {
        t2[i] = t[index_on_t++];
    }
}

int digit_val(char c){
    return 0;
}


void copy(int *b, int *a, int taille)
{
#pragma omp parallel for
    for (int i = 0; i < 2 * taille - 1; i++)
        b[i] = a[i];
}

int main(int argc, char **argv)
{

    FILE *fp;
    fp = fopen(argv[1], "r");
    if (!fopen)
    {
        fprintf(stderr, "coudln't open file name %s", argv[1]);
        exit(EXIT_FAILURE);
    }
    else
    {
        char ch;

        char buff[50];

        int MAX = 16 ;//2^10

        //we take the current time and store it in start
        //gettimeofday(&start, NULL);

        char op = '+';
        int n2;
        int value_test_log2 = (int)log2(MAX);
        int reste = MAX - pow(2, value_test_log2);
        //printf("reste = %d \n", reste);
        if (reste != 0)
        {
            if (isPowerOfTwo(reste) == 1)
            {
                n2 = pow(2, (int)log2(reste));
                //printf("n2 = %d\n", n2);
            }
            else
            {
                n2 = pow(2, (int)log2(reste) + 1);
                //printf("n2 = %d\n", n2);
            }
        }
        else
        {
            n2 = 0;
        }
        int n1 = pow(2, value_test_log2);
        int tab[MAX];
        int i = 0, j = 0;
        do
        {
            ch = fgetc(fp);
            buff[j++]=ch;
            if(ch==' ' || ch=='\t' || ch == '\n'){
                buff[j] = '\0';
                tab[i++] = atoi(buff);
                j = 0;
            }
        }while(ch!=EOF);    
        while(i<MAX)tab[i++]=elem_neutre(op);        

        printf("tab : \t");
        for(int i = 0;i<MAX; i++){
            printf("%d ", tab[i]);
        }
        printf("\n");
        int tab1[n1];
        int tab2[n2];

        remplir(tab1, tab2, tab, n1, n2);
        printf("tab1 : \t");
        for(int i = 0; i<n1; i++){
            printf("%d ", tab1[i]);
        }
        printf("\n");
        printf("tab2 : \t");
        for(int i = 0; i<n2; i++){
            printf("%d ", tab2[i]);
        }
        printf("\n");
        int size_a1 = pow(2, log2(n1) + 1) - 1;
        int *a1 = malloc(sizeof(int) * size_a1); //create_copieTable(tab1, n1, op);
        #pragma omp parallel for
        for (int i = 0; i < 2 * n1 - 1; i++)
            a1[i] = elem_neutre(op);
        #pragma omp parallel for
        for (int i = 0; i < n1; i++)
        {
            a1[(n1 - 1) + i] = tab1[i];
        }

        montee(a1, n1, op);
        int size_a2 = 0;
        if (reste != 0)
            size_a2 = pow(2, log2(n2) + 1) - 1;
        int *a2;
        if (reste != 0)
        {
            a2 = malloc(sizeof(int) * size_a2);
            #pragma omp parallel for
            for (int i = 0; i < 2 * n2 - 1; i++)
                a2[i] = elem_neutre(op);
            #pragma omp parallel for
            for (int i = 0; i < n2; i++)
            {
                a2[(n2 - 1) + i] = tab2[i];
            }
            montee(a2, n2, op);
        }
        printf("a1 :\t");
        for(int i = 0; i<2*n1-1; i++){
            printf("%d ", a1[i]);
        }
        printf("\n");
        printf("a2 :\t");
        for(int i = 0; i<2*n2-1; i++){
            printf("%d ", a2[i]);
        }
        printf("\n");
        int size_b1 = pow(2, log2(n1) + 1) - 1;
        int *psum = malloc(sizeof(int) * size_b1); //create_copieTable(tab1, n1, op);
        #pragma omp parallel for
        for (int i = 0; i < 2 * n1 - 1; i++)
            psum[i] = elem_neutre(op);
        #pragma omp parallel for
        for (int i = 0; i < n1; i++)
        {
            psum[(n1 - 1) + i] = tab1[i];
        }
        int size_b2 = 0;
        if (reste != 0)
            size_b2 = pow(2, log2(n2) + 1) - 1;
        int *b2;
        if (reste != 0)
        {
            b2 = malloc(sizeof(int) * size_b2);
            //b2 = create_copieTable(tab2, n2, op);
            #pragma omp parallel for
            for (int i = 0; i < 2 * n2 - 1; i++)
                b2[i] = elem_neutre(op);
            #pragma omp parallel for
            for (int i = 0; i < n2; i++)
            {
                b2[(n2 - 1) + i] = tab2[i];
            }
            copy(b2, a2, n2);
        }
        copy(psum, a1, n2);
        printf("psum : \t");
        for(int i = 0; i<2*n1-1; i++){
            printf("%d ", psum[i]);
        }
        printf("\n");
        
        printf("psum : \t");
        for(int i = 0; i<2*n2-1; i++){
            printf("%d ", b2[i]);
        }
        printf("\n");

        descente(psum, a1, n1, op);
        printf("après descente : psum :\t");
        for(int i = 0; i<2*n1-1; i++){
            printf("%d ", psum[i]);
        }
        final(psum, a1, n1, op);
        if (reste != 0)
        {
            descente(b2, a2, n2, op);
            printf("\naprès descente : b2 :\t");
            for(int i = 0; i<2*n2-1; i++){
                printf("%d ", b2[i]);
            }
            final(b2, a2, n2, op);
        }

        
        printf("\naprès finale \psum : ");
        for(int i = 0; i<2*n1-1; i++){
            printf("%d ", psum[i]);
        }
        printf("\n");

        printf("\naprès finale \tb2 : ");
        for(int i = 0; i<2*n2-1; i++){
            printf("%d ", b2[i]);
        }
        printf("\n");

        free(a1);
        free(psum);
        if (reste != 0)
        {
            free(b2);
            free(a2);
        }
        //we  store the current time in end

        //timeval is a struct with 2 parts for time, one in seconds and the other in
        //microseconds. So we convert everything to microseconds before computing
        //the elapsed time
        /*printf("%ld\n", ((end.tv_sec * 1000000 + end.tv_usec)
		  - (start.tv_sec * 1000000 + start.tv_usec)));*/

        long sum_to_ret;
        if (reste != 0)
        {
            sum_to_ret = psum[size_b1 - 1] + b2[size_b2 - 1];
        }
        else
        {
            sum_to_ret = psum[size_b1 - 1];
        }
        fclose(fp);
        return sum_to_ret;
    }
}