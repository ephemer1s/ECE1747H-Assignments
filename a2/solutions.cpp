#include <stdio.h> // for printf
#include <stdlib.h> // for exit
#include <errno.h> // for perror
#include <math.h> // for sqrt
#include "hrtime.h"
#include <iostream>


#define STRIDE 1024
#define ITER 1000000
#define NUM_ELEM ((STRIDE)*(ITER))
#define CACHELINE_SIZE 64
#define WORDS_PER_CACHELINE  ((CACHELINE_SIZE)/sizeof(int))

int error(const char * msg)
{
    perror(msg);
    exit(1);
}

// Function that is applied to an element of the array to transform it
// After that, the transformed alement is copied into another element of the array by the caller
// In the simplest form this function just returns the element unmodified.
// Use some math function to give it more weight
int transform(int a)
{
    int n = (int)round(sqrt(a));
    return n*n;
}

void sequentialSolution(int *A)
{
    for (int n = STRIDE; n < NUM_ELEM; ++n)
        A[n] = transform(A[n - STRIDE]);
}

void wrongSolution(int *A)
{
    // Wrong solution that does not take care of dependencies
    #pragma omp parallel for
    for (int n = STRIDE; n < NUM_ELEM; ++n)
        A[n] = transform(A[n - STRIDE]);
}

//replace with your parallel solution here

void parallelSolution1(int *A){
    for (int i = 0; i < ITER; i++)
        #pragma omp parallel for
        for (int j = i * STRIDE; j < (i + 1) * STRIDE; j++)
            A[j + STRIDE] = transform(A[j]);
}

void parallelSolution2(int *A){
    #pragma omp parallel for
    for (int i = 0; i < STRIDE; i++) {
        for (int j = 0; j < ITER; j++) {
            A[i + STRIDE * (j + 1)] = transform(A[i + j * STRIDE]);
        }
    }
}


// We tried the third approach from the slides, but we were not able to implement
// recursion in the innermost loop, which caused failure.
/*
void parallelSolution3(int *A){
    #pragma omp parallel for
    for (int i = 1; i < ITER; i++) {
        int k = i * STRIDE;
        for (int j = 0; j < STRIDE; j++) {
            A[k + j] = transform(A[j]);
            for(int l = 1; l < i; l++){
                A[k + j] = transform(A[k + j]);
            }
        }
    }
}
*/


int main(int argc, char *argv[])
{
    int solution_select = atoi(argv[1]);
    
    double start_time, end_time;
  
    int *A = (int*)malloc(sizeof(int) * NUM_ELEM);
    if (!A)
        error("cannot allocate array");
    // Initiallize array
    #pragma omp parallel for
    for (int n = 0; n < NUM_ELEM; ++n)
        A[n] = n;
    // print values for the last 1024 block of elements
    for (int n = 0; n < STRIDE; ++n)
        printf("%d ", A[n + (ITER-1)*STRIDE]);
    printf("\n");

    start_time = getElapsedTime();
    
    // sequentialSolution(A);
    // call your parallel solution here
    //parallelSolution1(A);
    //parallelSolution2(A);
    //parallelSolution3(A);

    if(solution_select == 0){
        sequentialSolution(A);
    }else if(solution_select == 1){
        parallelSolution1(A);
    }else if(solution_select == 2){
        parallelSolution2(A);
    }else if(solution_select == 3){
        parallelSolution3(A);
    }else{
        std::cout << "No such solution\n";
        exit(-1);
    }

    end_time = getElapsedTime();
    
    // verify results for the last 1024 block of elements
    for (int n = 0; n < STRIDE; ++n)
        printf("%d ", A[n + (ITER-1)*STRIDE]);
    printf("\n");

    printf("Solution %d Code took %lf seconds \n", solution_select, end_time - start_time);
    
    return 0;
}

	
	
	
