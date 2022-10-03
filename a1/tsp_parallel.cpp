//////////////////////////////////////////////////////////////////////////
// Traveling Salesman Problem.
// Note: this is a C++ program: compile it with g++
//
// Starting city is assumed to be city 0.
//////////////////////////////////////////////////////////////////////////

// to compile, run:
// g++ -c tsp_parallel.cpp
// g++ -pthread -o parallel.exe tsp_parallel.cpp

#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <chrono>
#include <mutex>
#include <cstdlib>
#include <pthread.h>
#include <vector>
using namespace std;

const int MAX_CITIES = 20;
const int MAX_THREADS = 4;
mutex mutex_shortest;
int **Dist;			// Dist[i][j] =  distance from  i to j

// pthread_mutex_t queue = PTHREAD_MUTEX_INITIALIZER;
// pthread_mutex_t best = PTHREAD_MUTEX_INITIALIZER;
// pthread_t tid[NUM_THREADS];   // thread IDs

///////////////////////////////////////////////////////////////////////////

class Path 
{
 public:
  int numVisited;	// Number of cities in the partial path
  int length;			// Current length of partial path
  int city[MAX_CITIES];
  int numCities;

  // Array city[] is a permutation of all cities.
  // city[0]..city[numVisited-1] is the current partial path;
  // city[numVisited]..city[numCities-1] are the cities not yet in the path 

  Path(int n): length(0), numVisited(1), numCities(n)			// Initialize a Path with city 0 as visited
  { 
    for (int i = 0; i < numCities; i++) 
      city[i] = i;
  }
  Path (const Path& other) : numVisited(other.numVisited), length(other.length), numCities(other.numCities)
  {
    for (int i = 0; i < numCities; i++) 
      city[i] = other.city[i];
  }

  void AddCity (int i)		// Extends path by adding the ith city not yet visited. 
  {
    assert(numVisited <= i && i < numCities);
    // In order to keep the path as a permutation of all cities
    // we "add" a city on the path by swapping 
    swap(city[i], city[numVisited]);
 
    // visit city[numVisited]
    length += Dist[city[numVisited-1]][city[numVisited]];
    numVisited++;
  }

  void Print()
  {
    for (int i = 0; i < numVisited; i++) 
	    cout << ' ' << city[i];
    cout << "; length = " << length << endl;
  }
};

class Queue 
{
  static const int MAXQSIZE = 100000;
  int size;
  Path *path[MAXQSIZE]; // keeps pointers to objects of type Path

 public:
  Queue(): size(0) {};
  void Put(Path *P);
  Path *Get();
  int isEmpty() { return size==0; };
};

///////////////////////////////////////////////////////////////////////////

void Queue::Put(Path *P) 
{ 
  assert(size < MAXQSIZE); 
  path[size++] = P;
}

Path *Queue::Get()
{
  if (isEmpty())
    return NULL;
  // to decrease the size of the queue, move the last element
  Path *p = path[0]; 
  path[0] = path[--size];
  return p;
}

///////////////////////////////////////////////////////////////////////////

void Fill_Dist(int numCities)
{
  // create Distance Matrix, Fill with random distances, symmetrical.
  cout << "Distance matrix:\n";
  Dist = new int*[numCities];
  for (int i=0; i<numCities; i++) {
    Dist[i] = new int[numCities];
    for (int j=0; j<numCities; j++) {
      // Not necessary, but let's make Dist simmetric:
      if (i==j) 
	      Dist[i][j] = 0;
      else if (j<i) 
	      Dist[i][j] = Dist[j][i];
      else
	      Dist[i][j] = rand()%1000; 
      cout << Dist[i][j] << '\t';
    }
    cout << endl << endl;
  }
}

struct Params
{
  int numCities;
  Queue *Q;
  Path *shortestPath;
};

void* tsp (void* arg)
{ 
  Params *params = (Params *)arg;
  int numCities = params->numCities;
  Queue *Q = params->Q;
  Path *shortestPath = params->shortestPath;

  // For every partial path in the queue, extend the partial path with all
  // of the unvisited cities (one by one) and add these newly formed partial
  // paths back to the queue
  while (!Q->isEmpty()) 
  {
    Path *p = Q->Get(); 

    // For each city not yet visited, extend the path:
    for (int i = p->numVisited; i < numCities; i++) 
    {
      Path *p1 = new Path(*p);	// initially p1 is a clone of p
      p1->AddCity(i);

      // decide what to do with new path
      if (p1->numVisited == numCities) // If we visited them all
      {
	      // Add last hop to return to origin:
	      p1->length += Dist[p1->city[numCities-1]][0];

	      // update shortestPath, if p1 is better
	      if (p1->length < shortestPath->length) {
          lock_guard<mutex> guard(mutex_shortest);
          *shortestPath = *p1;
        }
	      delete p1;
      } 
      else 
      {
        if (p1->length > shortestPath->length)
          delete p1; // This path is bad; just discard it
        else
	        Q->Put(p1);
      }
    } // end for
    delete p;			// p is not needed any more
  }
  return NULL;
}


int main(int argc, char *argv[])
{
  /*========= read args =========*/
  if (argc!=2) {
    cout << "Usage: " << argv[0] << " num_cities\n";
    exit(-1);
  }  

  int num_cities = atoi(argv[1]);
  assert(num_cities <= MAX_CITIES);
  
  Fill_Dist(num_cities);			// initialize Distance matrix

  /*====== create threads ======*/
  int num_threads = min(MAX_THREADS, num_cities - 1);
  pthread_t threads[num_threads];


  // Path *P = new Path(num_cities);
  // Queue Q;
  // Q.Put(P);			// initialize Q with one zero-length path
  // Path Shortest(num_cities);
  // Shortest.length = INT_MAX;    // The initial Shortest path must be bad
  // Params params = {num_cities, &Q, &Shortest};


  vector<Queue> queues(num_threads);
  vector<Params> params(num_threads);
  Path shortest = Path(num_cities);


  auto startTime = chrono::steady_clock::now();
  
  // tsp(&params);
  for (int i = 0; i < num_threads; i++) {
    int loop_cnt = (num_cities - 1) / MAX_THREADS + ((((num_cities - 1) % MAX_THREADS) > i) ? 1 : 0);
    for (int j = 0; j < loop_cnt; j++) { 
      Path *P = new Path(num_cities);
      P -> AddCity(i + j * MAX_THREADS + 1);
      queues[i].Put(P);
    }

    shortest.length = INT_MAX;
    params[i] = {num_cities, &queues[i], &shortest};
    
    int rc = pthread_create(&threads[i], NULL, tsp, (void *)&params[i]);
    if (rc) {
      cout << "Error:unable to create thread," << rc << endl;
      exit(-1);
    } 
  }   
  for (int i = 0; i < num_threads; i++) {
    pthread_join(threads[i], NULL);
  }
  
  
  auto endTime = chrono::steady_clock::now();
  auto ms = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);


  /*====== print solution ======*/
  cout << "Shortest path:";
  shortest.Print();
  
  cout << "TSP parallel solution took " << ms.count() << " ms\n";
  return 0;
}


// note to self: 
// to compile: https://www.geeksforgeeks.org/compiling-with-g-plus-plus/