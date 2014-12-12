/*    module load openmpi/gcc/64/1.8      //to setup mpicc
 * Test this example with 
 * 	mpirun -n $X ./single_MPI
 * where for X in {2,4,8,16}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>


#define COLS (5)
#define ROWS (8)

extern double get_change(double *x, double *y, int n);
extern void relax(double *dest,double *srce, int cols, int rows, int rowStart);
extern void init_grid(double **,double **,int cols, int rows);
extern void init_boundaries(double *,int,int rows);
extern void print_buffer(double *p, int cols, int rowStart, int rowEnd);


int main(int argc, char **argv) {

	double  *p,*p_new;
	init_grid(&p,&p_new,COLS,ROWS+2);

	init_boundaries(p,COLS,ROWS+2); 
	memmove(p_new,p, COLS*(ROWS+2) * sizeof(double) );
	
	int rank, size;   // MPI_START
 	MPI_Status status;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size); 

	int rowStart = (ROWS+2)/size*rank -1;            if (rank == 0) rowStart =0;
	int rowFinish =(ROWS+2)/size + (ROWS+2)/size*rank;	if (rank == size-1)rowFinish = ROWS+1;
	
//	printf ("%d %d %d \n", rank, rowStart, rowFinish);

//***********************************************************************
        double *send_end = malloc(COLS*sizeof(double));
        double *recv_end = malloc(COLS*sizeof(double));
        double *send_begin = malloc(COLS*sizeof(double));
        double *recv_begin = malloc(COLS*sizeof(double));

        int tag1 = 50;
        int tag2 = 51;


// Do some calculations.
	int count = 0;
	while (count++ < 300) {

               for (int j = 0; j < COLS; j++) send_end[j]=p[(rowFinish-1)*COLS+j];
               for (int j = 0; j < COLS; j++) send_begin[j]=p[(rowStart+1)*COLS+j];

//sending the second from the end to the next processor
	       if (rank != size-1){
                          MPI_Send(send_end, COLS, MPI_DOUBLE, rank+1, tag1, MPI_COMM_WORLD);
                          MPI_Recv(recv_begin, COLS, MPI_DOUBLE,rank+1, tag2, MPI_COMM_WORLD, &status); 
                          for (int j = 0; j < COLS; j++) p[(rowFinish)*COLS+j] =recv_begin[j] ;
		} 
	 
	      if (rank != 0) {
          		 MPI_Recv(recv_end, COLS, MPI_DOUBLE, rank-1, tag1, MPI_COMM_WORLD, &status); 
			 MPI_Send(send_begin, COLS, MPI_DOUBLE, rank-1, tag2, MPI_COMM_WORLD);		
			 for (int j = 0; j < COLS; j++) p[(rowStart)*COLS+j] =recv_end[j] ; 
		}
  
              relax(p_new,p,COLS,rowFinish+1, rowStart+1);
	      memmove(p,p_new, COLS*(ROWS+2) * sizeof(double) );

	}

//gathering info from different processors and printing the result
	double *p_res = malloc (COLS*(ROWS+2)*sizeof(double));
	MPI_Reduce(p,p_res, COLS*(ROWS+2) ,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if (rank ==0) {
		p_res[COLS/2] =p_res[COLS/2]/size;
		print_buffer(p_res, COLS, 0, ROWS+1);
        }

//**************************************************************
         MPI_Finalize();

	free (send_end);
	free (send_begin);
	free (recv_end);
	free (recv_begin);
	free(p_new);
	free(p); 

	return 0; 
}












void relax(double *new,double *old, int cols, int rows, int rowStart){
	for ( int j = rowStart ; j < (rows)-1; j++) {
		for ( int i = 1 ; i < cols-1; i++) {
			new[i+j*cols] = 0.25*(old[i-1+j*cols]+old[i+1+j*cols]+old[i+(j-1)*cols]+old[i+(j+1)*cols]);
		}
	}


}

void init_boundaries(double *l_p,int cols,int rows){
/* 
	for ( int i = 0 ; i < cols; i++) {
		l_p[i] = -1;
		l_p[i + cols*(rows -1)] = 1;
	}

	for ( int i = 0 ; i < rows; i++) {
		l_p[i*cols] = -1;
		l_p[(i+1)*cols -1] = 1; //column
		// l_p[i*(cols-1)+rows-1] = -1; //diagonal
	}
*/
	for (int i = 0 ; i < cols*rows; i++ )
		l_p[i] = 0;
	l_p[cols/2] = 1;
}


void init_grid(double **p, double **p_new,int cols, int rows){
	if (NULL == (*p = malloc(  cols*rows * sizeof(double) )  )  ) {
		puts("Allocation Error.");
		exit(99);
	}
	if (NULL == (*p_new = malloc(cols*rows * sizeof(double)))) {
		puts("Allocation Error.");
		exit(99);
	}
	for ( int i = 0 ; i < cols*rows; i++){

		(*p)[i] = 0. ;
		(*p_new)[i] = 0. ;
	}
}

void print_buffer(double *p, int cols, int rowStart, int rowEnd ){
	for (int i = rowStart; i < rowEnd+1; i++) {
		for (int j = 0; j < cols; j++) {
			printf("%f\t",p[i*cols+j]);
		}
		puts(" ");
	}


}

double get_change(double *x, double *y, int n){
	double result=0.;
	for (int i=0; i<n;i++)
		result+=x[i]*x[i]-y[i]*y[i];
	return result;
}

