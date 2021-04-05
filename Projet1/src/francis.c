#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <limits.h>

// The implementation of the dynamic table was found on the
// web site of professor EricK GALLESIO and has been modified to be conform to this algorithm :
// http://users.polytech.unice.fr/~eg/PP/Tds/Td08/corrige.html
// Everything else is written by :
// © Anas FRANCIS
// © 04-04-2021
// © Copyright 2021 Anas FRANCIS


#define ALLOC_SIZE 4

struct t_table {                // Implémentation avec un tableau
    long *tab;
    int size;
    int nb_elem;
};

typedef struct t_table *Table;

// prototype of functions :
void descente(long*, long *, int, char);
void final(long*, long *, int, char);
void montee(long*, int, char);
long operation(long, long, char);
long elem_neutre(char);
void reverse(long*, int);
void padding(Table*, char);
void fill_up(Table*, Table*, int);
void initialize(Table*, Table*,int, char);
void copy(long*, long*, int);

//
// Creation d'une table vide
//
Table creer_table(int alloc_size) {
    Table t         = malloc(sizeof(struct t_table));
    long *tab = malloc(alloc_size * sizeof(long));

    if (!t || !tab) {
        fprintf(stderr, "cannot allocate memory\n");
        return 0;
    }
    t->tab = tab;
    t->size = alloc_size;
    t ->nb_elem = 0;
    return t;
}


//
// Ajout dans un tableau (ajustable)
// (Le résultat est le nombre d'occurences de elt)

static int add_word(Table *table, int index, int my_nb) {
    Table T           = *table;
    long *numbers = T->tab;
    int nb        = T->nb_elem;
    int size         = T->size;

    if (nb == size) {
        // La table est pleine, la "rallonger" avant d'essayer d'insérer my_nb
        size *= 2;
        numbers = realloc(numbers, size*sizeof(long));

        if (!numbers) {
            fprintf(stderr, "cannot reallocate memory\n");
            return 0;
        }

        // conserver les nouvelles valeurs dans la table
        T->tab = numbers;
        T->size  = size;
    }

    // Décaler l'intervalle [i .. nb[ d'une case à droite en commençant par la dérnière valeur
    //Donc on ne peut pas parallèliser cette boucle.
    for (int i=nb; i > index; i--) {
        numbers[i] = numbers[i-1];
    }

    // InsÃ©rer le nouveau number à  la position index
    numbers[index]= my_nb;

    // On a un number de plus dans la table
    T->nb_elem += 1;

    return 1;                     // car ce number apparaît une fois
}

//
// Destruction d'une table
//
void detruire_table(Table *table) {
    Table T = *table;

    // libérer le tableau de numbers (le tableau allouée dans le t_table)
    free(T->tab);

    //libérer la table elle-même
    free(T);

    *table = NULL; // pas vraiment utile, mais pourquoi pas?
}


int ajouter_table(Table *table, int elt, int index) {
    return add_word(table, index, elt);
}

void reverse(long *table, int n){
    #pragma omp parallel for
    for(int i = 0; i<n/2;i++){
        int tmp = table[i];
        table[i] = table[n-1-i];
        table[n-1-i] = tmp;
    }
}

// O(1)
long elem_neutre(char op)
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
            return LONG_MIN;

        default:
            fprintf(stderr, "Error. Invalid operator !\nPrograme ended with code error -1.\n\a");
            exit(EXIT_FAILURE);
    }
}

// O(1)
void padding(Table* table, char op){
    Table t = *table;
    #pragma omp parallel for
    for(int i = t->nb_elem;i<=t->size;i++){
        t->tab[i]=elem_neutre(op);
    }
}

// O(1)
long operation( long val1, long val2, char op)
{
    switch (op)
    {
        case '+':
            return val1 + val2;
        case '-':
            return val1 - val2;
        case '*':
            return val1 * val2;
        case '/':
            if (val2 != 0)
                return val1 / val2;
            else
            {
                fprintf(stderr, "Error. Dividing by zero !\nProgramme ended with code error -1.\n\a");
                exit(EXIT_FAILURE);
            }
        case 'm':
            if(val1>val2)return val1;
            else return val2;

        default:
            fprintf(stderr, "Error. Invalid operator !\nProgramme ended with code error -1.\n\a");
            exit(EXIT_FAILURE);
    }
}

// O(log(N)) in parallel & O(Nlog(N)) in normal
void descente( long *b, long *a, int taille, char op)
{
    int m = log2(taille);
    b[0] = elem_neutre(op);

    for (int l = 1; l < m + 1; l++)
    {
        int maxi = pow(2, l + 1);
        #pragma omp parallel for
        for (int j = pow(2, l) - 1; j < maxi - 1; j++)
        {
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


//O(log(N) N is the size of input tab (here it is taille)
void montee( long *tab, int taille, char op)
{
    int m = log2(taille);
    for (int l = m - 1; l >= 0; l--)
    {
        int max = pow(2, l + 1) - 1;
        #pragma omp parallel for
        for (int j = pow(2, l) - 1; j < max; j++)
            tab[j] = operation(tab[2 * j + 1], tab[2 * j + 2], op);
    }
}

//O(1) because one loop in parallel each process compute one case of the table b.
void final(long *b, long *a, int taille, char op)
{
    int m = log2(taille);
    int maxi = pow(2, m + 1) - 1;

    #pragma omp parallel for
    for (int j = pow(2, m) - 1; j < maxi; j++)
        b[j] = operation(b[j], a[j], op);

}

// Amortize in O(1)
// O(N)
void fill_up(Table *t1, Table *tab, int n1)
{
    Table t = *tab;
    int index_on_t = 0;
    for (int i = 0; i < n1; i++)
    {
        ajouter_table(t1, t->tab[index_on_t++], i);
    }
}

// O(N)
void initialize(Table *a_b_i, Table * tabi,int ni, char op){

    Table t = *tabi;
    Table t_a_b_i = *a_b_i;
    int double_size= 2 * ni - 1;

    for (int i = 0; i <double_size ; i++){
        if(t_a_b_i->size<double_size){
            ajouter_table(a_b_i, elem_neutre(op), i);
        }
        else{
            t_a_b_i->tab[i] = elem_neutre(op);
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < ni; i++)
    {
        t_a_b_i->tab[(ni - 1) + i] =  t->tab[i];
    }
}

//O(1) because in parallel.
void copy( long *b, long *a, int taille)
{
    #pragma omp parallel for
    for (int i = 0; i < 2 * taille - 1; i++)
        b[i] = a[i];
}

int main(int argc, char **argv)
{
    if(argc!=2){
        if(argc<2)fprintf(stderr, "Not enough arguments ... \n\nProgramme exit with code error -1.\n");
        else fprintf(stderr, "Too many arguments ... \n\nProgramme exit with code error -1.\n");
        exit(EXIT_FAILURE);
    }

    FILE *fp = fopen(argv[1], "r");

    if (!fp){
        fprintf(stderr, "couldn't open file name %s...\n\nProgramme exit with code error -1.\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    else{

        // Start reading file and fill the dynamic table with numbers in the file fp !
        char ch, buff[50], op = '+';
        int i, j = 0, nb_elements = 0;

        Table dyn_table = creer_table(ALLOC_SIZE);

        do
        { // filling dyn_table with elements of fp file is O(N).
            ch = fgetc(fp);
            buff[j++]=ch;

            if(ch==' ' || ch=='\t' || ch == '\n'){
                buff[j] = '\0';
                int value = atoi(buff);
                ajouter_table(&dyn_table, value, nb_elements++);// ajouter_table will expand the dynamic table
                // if there is no more place.
                j = 0;
            }

        }while(ch!=EOF);
        fclose(fp);

        int integer_part_of_log2 = (int)log2(nb_elements);
        int reste = nb_elements - pow(2, integer_part_of_log2);

        int n1 = pow(2, integer_part_of_log2);
        if (reste != 0) // Padding dyn_table if nb_elements is not a power of 2.
            // size of dyn_table has been doubled we need to do a padding to the end
            // with elem_neutre(op) function.
            // O(1) because loop in parallel.
        {
            padding(&dyn_table, op);
            nb_elements = dyn_table->size;
            n1 = dyn_table->size;
        }
        int size_b1 = pow(2, log2(n1) + 1) - 1;
        Table tab1_bis = creer_table(ALLOC_SIZE);

        fill_up(&tab1_bis, &dyn_table , n1);
        Table a1 = creer_table(ALLOC_SIZE);

        initialize(&a1,&tab1_bis, n1, op);
        /**
         * Compute sum-prefix of Q and store them in array PSUM
        */
        montee(a1->tab, n1, op);
        Table b1 = creer_table(ALLOC_SIZE);
        initialize(&b1,&tab1_bis,n1,op);
        copy(b1->tab, a1->tab, n1);
        descente(b1->tab, a1->tab, n1, op);
        final(b1->tab, a1->tab, n1, op);
        Table psum = creer_table(ALLOC_SIZE);
        i = 0;

        // O(N) fill the psum table.
        while(i<nb_elements){
            ajouter_table(&psum, b1->tab[size_b1-1-i],i);
            i++;
        }

        /**
         * Compute sum-suffix of Q and store them in array SSUM
        */

        Table ssum= creer_table(ALLOC_SIZE);
        reverse(tab1_bis->tab, nb_elements);
        initialize(&a1, &tab1_bis, n1, op);
        montee(a1->tab, n1, op);
        initialize(&b1, &dyn_table, n1, op);
        copy(b1->tab, a1->tab, n1);
        descente(b1->tab, a1->tab, n1, op);
        final(b1->tab, a1->tab, n1, op);
        i = 0;

        // O(N) fill the ssum table.
        while(i<nb_elements){
            ajouter_table(&ssum, b1->tab[size_b1-1-i],i);
            i++;
        }
        op = 'm'; // change the operation to Max
        /**
         * Compute max-suffix of PSUM and store them in array SMAX
         */

        Table smax= creer_table(ALLOC_SIZE);
        initialize(&a1, &psum, n1, op);
        montee(a1->tab, n1, op);
        initialize(&b1, &psum, n1, op);
        copy(b1->tab, a1->tab, n1);
        descente(b1->tab, a1->tab, n1, op);
        final(b1->tab, a1->tab, n1, op);
        i=0;

        // O(N) fill the smax table.
        while(i<nb_elements){
            ajouter_table(&smax,b1->tab[size_b1-1-i],i);
            i++;
        }
        /**
         * Compute max-prefix of SSUM and store them in array PMAX
        */

        Table pmax= creer_table(ALLOC_SIZE);
        initialize(&a1, &ssum, n1, op);
        montee(a1->tab, n1, op);
        initialize(&b1, &ssum, n1, op);
        copy(b1->tab, a1->tab, n1);
        descente(b1->tab, a1->tab, n1, op);
        final(b1->tab, a1->tab, n1, op);
        i=0;

        // O(N) fill the pmax table.
        while(i<nb_elements){
            ajouter_table(&pmax, b1->tab[size_b1-1-i],i);
            i++;
        }
        // The next 2 reverse operations are in O(1).
        // Reverse the psum and pmax table to compute M.
        // Because the ajouter_table(...) function add the elements in the reverse order.
        reverse(psum->tab, nb_elements);
        reverse(pmax->tab, nb_elements);

        // Create a dynamic table of size nb_elements which is a power of 2, to compute M.
        Table M = creer_table(nb_elements);

        long maximum_val = elem_neutre('m');
        int index_of_max_start = -1, end_index = 0; // index of start and index of end in the Q table

        /**
         * for 1 <= i <= n do in parallel
         * Ms[i] = pmax[i] - ssum[i] + Q[i]
         * Mp[i] = smax[i] - psum[i] + Q[i]
         * M[i] = Ms[i] + Mp[i] - Q[i]
         */

        #pragma omp parallel for
        for(int i = 0;i<nb_elements;i++){ // O(1) loop in parallel
            M->tab[i] = pmax->tab[i] - ssum->tab[i]+smax->tab[i] - psum->tab[i] + dyn_table->tab[i];
        }
        // Search the Max in O(N)
        for(int i = 0;i<nb_elements;i++){
            if(M->tab[i] > maximum_val){
                maximum_val = M->tab[i];
                index_of_max_start = i;
                end_index = index_of_max_start;
                if(maximum_val<0) {
                    end_index++;
                }
                else{
                    while(M->tab[end_index] == maximum_val){
                        end_index++;
                    }
                }
                i = end_index-1;
            }
        }

        // Print the maximum sub_table O(K); K size of sub_table
        // 1<= K <= N
        printf("%ld ",M->tab[index_of_max_start]);
        for (i=index_of_max_start;i<end_index;i++)
        {
            if(i<end_index-1)printf("%ld ", dyn_table->tab[i]);
            else{
                printf("%ld\n", dyn_table->tab[i]);
            }
        }

        // Free the heap.
        detruire_table(&a1);
        detruire_table(&b1);
        detruire_table(&tab1_bis);
        detruire_table(&dyn_table);
        detruire_table(&psum);
        detruire_table(&ssum);
        detruire_table(&pmax);
        detruire_table(&smax);
        detruire_table(&M);

        return 0;
    }
}