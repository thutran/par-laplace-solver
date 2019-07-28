#include <upc_relaxed.h>
// #include <upc_collective.h>
#include <bupc_collectivev.h>
// #include <upc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <ctype.h>
#include <unistd.h> // getopt

#ifndef MIN
#define MIN(x,y)  ((x)<(y)?(x):(y))
#endif

#if !defined(NPES)
#define NPES   4        // number of processors
#endif

// max size of plate
#define TOTAL_COLUMNS 1000
#define MAX_ROWS 1000

#define ROWS (MAX_ROWS/THREADS)  
#define LOCAL_SIZE MIN(ROWS*TOTAL_COLUMNS, UPC_MAX_BLOCK_SIZE) 
#define LOCAL_ROWS LOCAL_SIZE/TOTAL_COLUMNS // number of local rows at each thread
#define TOTAL_ROWS THREADS*LOCAL_ROWS

#define MAX_TEMP_ERROR 0.01

// temperature arrays
shared [LOCAL_SIZE] double Temperature[TOTAL_ROWS][TOTAL_COLUMNS]; // extra bottom row for thread 0, bottom row as horizontal heating element at the last thread
shared [LOCAL_SIZE] double Temperature_last[TOTAL_ROWS][TOTAL_COLUMNS]; // extra left column for all cells furthest to the heating element, right column as the vertical heating element
shared double *temp_last_global = &Temperature_last;

// heating elements
shared [LOCAL_ROWS] double Heating_vertical[TOTAL_ROWS]; // to the right of the plate
// double Heating_horizontal[TOTAL_COLUMNS]; // will be initialized only at the bottom plate
// shared [0] double Heating_horizontal[TOTAL_COLUMNS]; // give affinity to thread 0
shared [] double *Heating_horizontal;

shared int max_iterations = 100;
shared double dt_global=100.0;


// to be called by only 1 thread (thread 0)
void initialize_globally(){
  for (int i = 0; i < TOTAL_ROWS; i++){
    for (int j = 0; j < TOTAL_COLUMNS; j++){
      Temperature[i][j] = 0.0;    
      Temperature_last[i][j] = 0.0;    
    }
  }
}

void initialize_locally(){
  // Local boundry condition endpoints
  double tMin = (MYTHREAD)*100.0/THREADS;
  double tMax = (MYTHREAD+1)*100.0/THREADS;

  // Left and right boundaries
  for (int i = 0; i < LOCAL_ROWS; i++) {
    Heating_vertical[MYTHREAD*LOCAL_ROWS + i] = tMin + ((tMax-tMin)*(i+1)/(LOCAL_ROWS + 0.0));
  }

  // Horizontal heating element to be used by the last thread only
  if (MYTHREAD == THREADS - 1){
    Heating_horizontal = upc_alloc(TOTAL_COLUMNS * sizeof(double));
    for (int i = 0; i < TOTAL_COLUMNS; i++){
      Heating_horizontal[i] = (100.0/(TOTAL_COLUMNS+0.0)) * (i+1);
    } 
  }
}

// only called by last thread
void track_progress(int iteration) {
    printf("---------- Iteration number: %d ------------\n", iteration);

    // output global coordinates so user doesn't have to understand decomposition
    for(int i = MIN(5, TOTAL_ROWS); i > 0; i--) {
      printf("[%d,%d]: %5.2f  ", TOTAL_ROWS-i,TOTAL_COLUMNS-i, Temperature[TOTAL_ROWS-i][TOTAL_COLUMNS-i]);
    }
    printf("\n");
}


double temp_neighbor_below_at(int column_index)  {
  if (MYTHREAD == 0)
    return 0.0;
  
  int index_in_global_array = (TOTAL_COLUMNS*THREADS) * (LOCAL_ROWS-1) + column_index*THREADS + MYTHREAD -1;

  return temp_last_global[index_in_global_array];
}

double temp_neighbor_above_at(int column_index)  {
  if (MYTHREAD == THREADS - 1)
    return Heating_horizontal[column_index];
  
  int index_in_global_array = column_index*THREADS + MYTHREAD + 1;
  
  return temp_last_global[index_in_global_array];
}

//-------------------------------------MAIN--------------------------------------//

int main(int argc, char **argv) 
{
  // verify only NPES PEs are being used
  if(NPES != THREADS) {
    if(MYTHREAD == 0) {
      printf("This code must be run with %d PEs\n", NPES);
    }
    exit(1);
  }
  
  struct timeval start_time, stop_time, elapsed_time;

  char *max_inter_arg = NULL;
    char *sum_name = NULL;
    int c;
    int quiet =0;

    // BEGIN arguments--------------------
    opterr = 0; // already declared in the environment

    while ((c = getopt(argc, argv, "qm:s:")) != -1)
        switch (c)
        {
        case 'q':
          quiet = 1;
          break;
        case 'm':
          max_inter_arg = optarg;
          break;
        case 's':
            sum_name = optarg;
            break;
        case '?':
            if (optopt == 'm')
                fprintf(stderr, "Option -%c requires an argument for max number of iterations.\n", optopt);
            else if (optopt == 's')
                fprintf(stderr, "Option -%c requires an argument for name of the summary file.\n", optopt);
            else if (isprint(optopt))
                fprintf(stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf(stderr,
                        "Unknown option character `\\x%x'.\n",
                        optopt);
            return 1;
        default:
            abort();
        }

  max_iterations = max_inter_arg ? atoi(max_inter_arg) : 0;

  if (MYTHREAD == 0 && !max_iterations)
  {
    // initialize_globally();
    int max_iter;
    printf("Maximum iterations [100-4000]?\n");
    fflush(stdout); // Not always necessary, but can be helpful
    // scanf("%d", max_iterations);
    scanf("%d", &max_iter);
    max_iterations = max_iter;
  }
  // initialize_locally();
  // upc_barrier; // to make sure max_iterations is updated

  FILE *fsum = sum_name ? fopen ( sum_name, "a" ) : NULL;
  // END arguments--------------------


  int iteration=1;
  double *temp_local = (double*) &Temperature ; //[MYTHREAD*LOCAL_ROWS][0];
  double *temp_local_last = (double*) &Temperature_last; //[MYTHREAD*LOCAL_ROWS][0];
  double *temp_heating_vertial_local = (double*) &Heating_vertical;

  double dt_local = 0.0;
  double dt_sum;
  

  upc_barrier;

  if (MYTHREAD == 0)
  {
    // printf("UPC_MAX_BLOCK_SIZE %d\n", UPC_MAX_BLOCK_SIZE);
    gettimeofday(&start_time, NULL);
    initialize_globally();
  }
  initialize_locally();
  

  while (dt_global > MAX_TEMP_ERROR && iteration <= max_iterations ){
    
    if (LOCAL_ROWS == 1)
    {
      for (int j = 0; j < TOTAL_COLUMNS; j++)
      {
        int cell_index = j;

        /* furthest column from the vertical heating element*/
        if (j == 0)
        {
          temp_local[cell_index] = 0.25 * (temp_local_last[cell_index + 1] + temp_neighbor_below_at(j) + temp_neighbor_above_at(j));
        }
        /* closest column to the vertical heating element*/
        else if (j == TOTAL_COLUMNS - 1)
        {
          // temp_local[cell_index] = 0.25 * (temp_local_last[j] + temp_heating_vertial_local[j] + temp_neighbor_above_at(j) );
          temp_local[cell_index] = 0.25 * (temp_local_last[cell_index - 1] + temp_heating_vertial_local[0] + temp_neighbor_below_at(j) + temp_neighbor_above_at(j) );
        }
        /* middle columns */
        else
        {
          temp_local[cell_index] = 0.25 * (temp_local_last[cell_index - 1] + temp_local_last[cell_index + 1] + temp_neighbor_below_at(j) + temp_neighbor_above_at(j) );
        }
      }
    } else
    {
      for (int i = 0; i < LOCAL_ROWS; i++)
      {
        for (int j = 0; j < TOTAL_COLUMNS; j++)
        {
          int cell_index = i * TOTAL_COLUMNS + j;

          /* first row, needs temp from the below neighbor */
          if (i == 0)
          {
            /* furthest column from the vertical heating element*/
            if (j == 0)
            {
              temp_local[cell_index] = 0.25 * (temp_local_last[cell_index + 1] + temp_local_last[cell_index + TOTAL_COLUMNS] + temp_neighbor_below_at(j));
            }
            /* closest column to the vertical heating element*/
            else if (j == TOTAL_COLUMNS - 1)
            {
              // temp_local[cell_index] = 0.25 * (temp_local_last[j] + temp_heating_vertial_local[j] + temp_neighbor_above_at(j) );
              temp_local[cell_index] = 0.25 * (temp_local_last[cell_index - 1] + temp_heating_vertial_local[i] + temp_local_last[cell_index + TOTAL_COLUMNS] + temp_neighbor_below_at(j));
            }
            /* middle columns */
            else
            {
              temp_local[cell_index] = 0.25 * (temp_local_last[cell_index - 1] + temp_local_last[cell_index + 1] + temp_local_last[cell_index + TOTAL_COLUMNS] + temp_neighbor_below_at(j));
            }
          }
          /* last row, needs temperature from the above neighbor */
          else if (i == LOCAL_ROWS - 1)
          {
            /* furthest column from the vertical heating element*/
            if (j == 0)
            {
              temp_local[cell_index] = 0.25 * (temp_local_last[cell_index + 1] + temp_local_last[cell_index - TOTAL_COLUMNS] + temp_neighbor_above_at(j));
            }
            /* closest column to the vertical heating element*/
            else if (j == TOTAL_COLUMNS - 1)
            {
              temp_local[cell_index] = 0.25 * (temp_local_last[cell_index - 1] + temp_heating_vertial_local[i] + temp_local_last[cell_index - TOTAL_COLUMNS] + temp_neighbor_above_at(j));
            }
            /* middle columns */
            else
            {
              temp_local[cell_index] = 0.25 * (temp_local_last[cell_index - 1] + temp_local_last[cell_index + 1] + temp_local_last[cell_index - TOTAL_COLUMNS] + temp_neighbor_above_at(j));
            }
          }
          else
          {
            /* furthest column from the vertical heating element*/
            if (j == 0)
            {
              // temp_local[cell_index] = 0.25 * (temp_local_last[j] + temp_neighbor_below_at(j) + temp_neighbor_above_at(j) );
              temp_local[cell_index] = 0.25 * (temp_local_last[cell_index + 1] + temp_local_last[cell_index - TOTAL_COLUMNS] + temp_local_last[cell_index + TOTAL_COLUMNS]);
            }
            /* closest column to the vertical heating element*/
            else if (j == TOTAL_COLUMNS - 1)
            {
              // temp_local[cell_index] = 0.25 * (temp_local_last[j] + temp_heating_vertial_local[j] + temp_neighbor_above_at(j) );
              temp_local[cell_index] = 0.25 * (temp_local_last[cell_index - 1] + temp_heating_vertial_local[i] + temp_local_last[cell_index - TOTAL_COLUMNS] + temp_local_last[cell_index + TOTAL_COLUMNS]);
            }
            /* middle columns */
            else
            {
              temp_local[cell_index] = 0.25 * (temp_local_last[cell_index - 1] + temp_local_last[cell_index + 1] + temp_local_last[cell_index - TOTAL_COLUMNS] + temp_local_last[cell_index + TOTAL_COLUMNS]);
            }
          }
        }
      }
    }
    
    dt_local = 0.0;

    upc_barrier;
    // calculate temperature changes
    for (int i = 0; i < LOCAL_SIZE; i++){
      dt_local = fmax(fabs(temp_local[i] - temp_local_last[i]), dt_local);
      
      // printf(">> thread %d iteration %d update temp cell %d\n", MYTHREAD, iteration, i);

      temp_local_last[i] = temp_local[i];
      
    }

    dt_sum = bupc_allv_reduce(double, dt_local, 0, UPC_MAX); // max dt across all threads
    if (MYTHREAD ==0)
    {
      dt_global = dt_sum; 
    }

    // periodically print test values - only for thread in lower corner
    if ((iteration % 100) == 0)
    {
      if (MYTHREAD == THREADS - 1 && !quiet)
      {
        track_progress(iteration);
      }
    }

    ++iteration;
    // printf("-------- thread %d END iteration %d --------\n", MYTHREAD, iteration-1);
  }

  upc_barrier; // for more accurate time calculation


  if (MYTHREAD == 0){
    // printf("iteration %d\t dt_local dt %f\t dt_global dt %f\n", iteration, dt_local, dt_global);
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time);
    
    if (!quiet)
    {
      printf("\nMax error at iteration %d was %f\n", iteration - 1, dt_global);
      printf("Total time was %f seconds.\n", elapsed_time.tv_sec + elapsed_time.tv_usec / 1000000.0);
    }

    // Printing summary data
    // version   nproc   max_iterations stop_at_iteration   gloabl_dt   time
    if (fsum)
      fprintf(fsum, "upc\t%d\t%d\t%d\t%f\t%g\n", THREADS, max_iterations, iteration - 1, dt_global, elapsed_time.tv_sec + elapsed_time.tv_usec / 1000000.0);

  }

  if (MYTHREAD == THREADS - 1 && !quiet)
  {
    track_progress(iteration-1);
  }

  // Clearing space
  if (fsum)
    fclose(fsum);

  return(0); 
}

