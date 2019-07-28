/*************************************************
 * Laplace Serial C Version
 *
 * Temperature is initially 0.0
 * Boundaries are as follows:
 *
 *      0         T         0
 *   0  +-------------------+  0
 *      |                   |
 *      |                   |
 *      |                   |
 *   T  |                   |  T
 *      |                   |
 *      |                   |
 *      |                   |
 *   0  +-------------------+ 100
 *      0         T        100
 *
 *  John Urbanic, PSC 2014
 *
 ************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <ctype.h>
#include <unistd.h> // getopt

// size of plate
#define COLUMNS    1000
#define ROWS       1000

// largest permitted change in temp (This value takes about 3400 steps)
#define MAX_TEMP_ERROR 0.01

double Temperature[ROWS+2][COLUMNS+2];      // temperature grid
double Temperature_last[ROWS+2][COLUMNS+2]; // temperature grid from last iteration

//   helper routines
void initialize();
void track_progress(int iter);


int main(int argc, char **argv)
{
    int i, j;                                            // grid indexes
    // int max_iterations;                                  // number of iterations
    int iteration=1;                                     // current iteration
    double dt=100;                                       // largest change in t
    struct timeval start_time, stop_time, elapsed_time;  // timers

    
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

    int max_iterations = max_inter_arg ? atoi(max_inter_arg) : 0;

    FILE *fsum = sum_name ? fopen ( sum_name, "a" ) : NULL;
    // END arguments--------------------

    
    if (!max_iterations)
    {
        printf("Maximum iterations [100-4000]?\n");
        scanf("%d", &max_iterations);
    }

    gettimeofday(&start_time,NULL); // Unix timer

    initialize();                   // initialize Temp_last including boundary conditions

    // do until error is minimal or until max steps
    while ( dt > MAX_TEMP_ERROR && iteration <= max_iterations ) {

        // main calculation: average my four neighbors
        for(i = 1; i <= ROWS; i++) {
            for(j = 1; j <= COLUMNS; j++) {
                Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] +
                                            Temperature_last[i][j+1] + Temperature_last[i][j-1]);
            }
        }
        
        dt = 0.0; // reset largest temperature change

        // copy grid to old grid for next iteration and find latest dt
        for(i = 1; i <= ROWS; i++){
            for(j = 1; j <= COLUMNS; j++){
	      dt = fmax( fabs(Temperature[i][j]-Temperature_last[i][j]), dt);
	      Temperature_last[i][j] = Temperature[i][j];
            }
        }

        // periodically print test values
        if ((iteration % 100) == 0 && !quiet)
        {
            track_progress(iteration);
        }

        iteration++;
    }

    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

    if (!quiet)
    {
        printf("\nMax error at iteration %d was %f\n", iteration - 1, dt);
        printf("Total time was %f seconds.\n", elapsed_time.tv_sec + elapsed_time.tv_usec / 1000000.0);
    }


    // Printing summary data
    // version   nproc   max_iterations stop_at_iteration   gloabl_dt   time
    if( fsum) 
        fprintf(fsum,"serial\t1\t%d\t%d\t%f\t%g\n", max_iterations, iteration-1, dt, elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0); 
    // Clearing space
    if( fsum )
        fclose( fsum );   


    return 0;
}


// initialize plate and boundary conditions
// Temp_last is used to to start first iteration
void initialize(){

    int i,j;

    for(i = 0; i <= ROWS+1; i++){
        for (j = 0; j <= COLUMNS+1; j++){
            Temperature_last[i][j] = 0.0;
        }
    }

    // these boundary conditions never change throughout run

    // set left side to 0 and right to a linear increase
    for(i = 0; i <= ROWS+1; i++) {
        Temperature_last[i][0] = 0.0;
        Temperature_last[i][COLUMNS+1] = (100.0/ROWS)*i;
    }
    
    // set top to 0 and bottom to linear increase
    for(j = 0; j <= COLUMNS+1; j++) {
        Temperature_last[0][j] = 0.0;
        Temperature_last[ROWS+1][j] = (100.0/COLUMNS)*j;
    }
}


// print diagonal in bottom right corner where most action is
void track_progress(int iteration) {

    int i;

    printf("---------- Iteration number: %d ------------\n", iteration);
    for(i = ROWS-5; i <= ROWS; i++) {
        printf("[%d,%d]: %5.2f  ", i, i, Temperature[i][i]);
    }
    printf("\n");
}
