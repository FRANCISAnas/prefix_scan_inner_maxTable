# prefix_scan_inner_maxTable

See on github : https://github.com/FRANCISAnas/prefix_scan_inner_maxTable

By the time you read this readme.txt file the github repository will be public.

Project to find the maximum inner array of an array of 1D, this can be extended to be used for 2D arrays too.

By Anas Francis e-mail : anas.francis@etu.unice.fr

Don't use this code only for education purposes.

To compile and execute the code you must use openmp option :
"gcc -g -o output.mp francis.c -Wall -std=gnu99 -lm -fopenmp"
And the to execute the programme juste do : "./output.mp your_input_data_file"
Note that you must specify exactly ONE input file otherwise the code will exit with code error -1.

This Work contains a Makefile, if you're under a machine that has linux OS just launch "make".

Methode in parallel :
void descente(long*, long *, int, char); Inner loop in parallel.
void final(long*, long *, int, char); Inner loop in parallel.
void montee(long*, int, char); Inner loop in parallel.
void reverse(long*, int); In parallel
void padding(Table*, char); In parallel
void initialize(Table*, Table*,int, char); Contains 2 for loop the 2nd is in parallel
void copy(long*, long*, int); In parallel

And also the section to compute the elements of the final table M, this loop is in parallel.
