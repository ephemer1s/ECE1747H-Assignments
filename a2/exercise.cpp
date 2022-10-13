#include <stdio.h> // for printf
#include <stdlib.h> // for exit
#include <errno.h> // for perror
#include <math.h> // for sqrt
#include "hrtime.h"


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



int main()
{
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
    
    sequentialSolution(A);
    // call your parallel solution here
    end_time = getElapsedTime();
    
    // verify results for the last 1024 block of elements
    for (int n = 0; n < STRIDE; ++n)
        printf("%d ", A[n + (ITER-1)*STRIDE]);
    printf("\n");

    printf("Code took %lf seconds \n", end_time - start_time);
    
    return 0;
}

	
	
	
