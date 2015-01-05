#include "proc.h"
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include "fox.h"
#include <sys/time.h>

#define MO 1000000

// Initialize block type and arrays used in scatter and gather
void init_type(MPI_Datatype *blocktype, MPI_Datatype *blocktype2, int block_size, int msize, int nb_blocks, int* disps, int* counts)
{
  int i,jj;
  MPI_Type_vector(block_size, block_size, msize, MPI_DOUBLE, blocktype2);
  MPI_Type_create_resized( *blocktype2, 0, sizeof(double), blocktype);
  MPI_Type_commit(blocktype);

  for (i=0; i<nb_blocks; i++) {
    for (jj=0; jj<nb_blocks; jj++) {
      disps[i*nb_blocks+jj] = i*msize*block_size+jj*block_size;
      counts [i*nb_blocks+jj] = 1;
    }
  }
}

// scatter one matrix A to each process in A_local
void scatter(MPI_Comm final, int myrank, MPI_Datatype blocktype, int* disps, int* counts, struct matrix* A, struct matrix* A_local, int block_size)
{
  if (myrank == 0){
    A_local->size =  block_size;
  }
  MPI_Bcast(&A_local->size,1,MPI_INT,0,final);
  
  A_local->data = malloc(block_size * block_size * sizeof(double));
   
  MPI_Scatterv(A->data, counts, disps, blocktype, A_local->data, block_size*block_size, MPI_DOUBLE, 0, final);
}

void scatter_all(MPI_Comm final, int myrank, MPI_Datatype blocktype, int* disps, int* counts, struct matrix* A, struct matrix* B, struct matrix* C, struct matrix* A_local, struct matrix* B_local, struct matrix* C_local, int block_size)
{
  scatter(final,myrank,blocktype,disps,counts,A,A_local,block_size);
  scatter(final,myrank,blocktype,disps,counts,B,B_local,block_size);
  scatter(final,myrank,blocktype,disps,counts,C,C_local,block_size);
}

void gather_all(struct matrix* A_local, struct matrix* B_local, struct matrix* C_local, int block_size, struct matrix* A, struct matrix* B, struct matrix* C, int* counts, int* disps, MPI_Datatype blocktype, MPI_Comm final)
{
  MPI_Gatherv(A_local->data,  block_size*block_size, MPI_DOUBLE, A->data, counts, disps, blocktype, 0, final);
  MPI_Gatherv(B_local->data, block_size*block_size , MPI_DOUBLE, B->data, counts, disps, blocktype, 0, final);
  MPI_Gatherv(C_local->data,  block_size*block_size, MPI_DOUBLE, C->data, counts, disps, blocktype, 0, final);
}


// main function calling mpi functions
int proc(char *argv0, int gen_matrix, int block_size, int msize)
{
  int i,j;
  int myrank, nb_proc, nb_blocks, grid_rank[DIMENSION] = {42, 42}, src_rank[DIMENSION] = {42, 42}, dest_rank[DIMENSION] = {42, 42}, size;
  double initialTime = 0.0;
  char buffer[1024], buffer2[1024];
  char *spawn_params[4] = {"0", buffer, buffer2, NULL};
  struct matrix *A, *B;
  FILE* fd = NULL;
  
  MPI_Datatype blocktype;
  MPI_Datatype blocktype2;
  
  struct timeval t1, t2;
  unsigned long time_difference = 0;
  
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
    
  if (myrank == 0)
    {
      //if (gen_matrix == 1)
      //	{
	  A = rand_matrix(msize);
	  //write_matrix_file("A.mat", A);
	  B = rand_matrix(msize);
	  //write_matrix_file("B.mat", B);
	  /*	}
      else
	{
	  A = read_file("A.mat");
	  if (A == NULL){
	    fprintf(stderr, "Cannot load A.mat\n");
	    exit(EXIT_FAILURE);
	  }
	  B = read_file("B.mat");
	  if (B == NULL){
	    fprintf(stderr, "Cannot load B.mat\n");
	    exit(EXIT_FAILURE);
	  }
	  msize = A->size;
	  if (msize % block_size) {
	    fprintf(stderr, "The matrix size (%d) must be a multiple of block_size(%d)\n", msize, block_size);
	    exit(EXIT_FAILURE);
	  }
	  }*/
    }
  else
    {
      A = init_matrix(msize);
      B = init_matrix(msize);
    }
  struct matrix * C = init_matrix(msize);
      
  nb_blocks = msize / block_size; // number of blocks on one line
  nb_proc = nb_blocks * nb_blocks; // number of processus to be created = total number of blocks
    
  MPI_Comm grid, world, final = MPI_COMM_WORLD;
  sprintf(buffer, "%d", block_size);
  sprintf(buffer2, "%d", msize);
 
  MPI_Comm_get_parent(&world);
  
  if (size < nb_proc){
    if (world == MPI_COMM_NULL) {
      MPI_Comm_spawn(argv0, spawn_params, nb_proc-size, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &world, MPI_ERRCODES_IGNORE);
      MPI_Intercomm_merge(world, 0, &final);
    }
    else{
      MPI_Intercomm_merge(world, 1, &final);
    }
  }
   
  MPI_Comm_size(final, &size);
  MPI_Comm_rank(final, &myrank);
  
  int dims[DIMENSION] = {nb_blocks, nb_blocks}, periods[DIMENSION] = {1, 0}, reorder = 1;
  MPI_Cart_create(final, DIMENSION, dims, periods, reorder, &grid);

  MPI_Cart_coords(grid, myrank, DIMENSION, grid_rank);

  if (myrank == 0)
    initialTime = MPI_Wtime();

  // Calculs divers et variés sur la grille
  MPI_Cart_shift(grid, Y_DIMENSION, 1, src_rank, dest_rank);
  MPI_Cart_coords(grid, src_rank[0], DIMENSION, src_rank);
  MPI_Cart_coords(grid, dest_rank[0], DIMENSION, dest_rank);

  struct matrix A_local, B_local, C_local;
  int disps[nb_blocks*nb_blocks];
  int counts[nb_blocks*nb_blocks];
    
  if (myrank == 0)
    fprintf(stdout, "%d ", msize);
  
  init_type(&blocktype, &blocktype2, block_size, msize, nb_blocks, disps, counts);
  gettimeofday(&t1, NULL);
  scatter_all(final, myrank, blocktype, disps, counts, A, B, C, &A_local, &B_local, &C_local, block_size);
  gettimeofday(&t2, NULL);
  time_difference = (t2.tv_sec-t1.tv_sec)*1000000 +(t2.tv_usec-t1.tv_usec);
  if (myrank == 0)
    fprintf(stdout, "%lu ", time_difference);
  
  gettimeofday(&t1, NULL);
  fox_multiplication(grid, nb_blocks, &A_local, &B_local, &C_local);
  gettimeofday(&t2, NULL);
  time_difference = (t2.tv_sec-t1.tv_sec)*1000000 +(t2.tv_usec-t1.tv_usec);
  if (myrank == 0)
    fprintf(stdout, "%lu ", time_difference);
  /*
    print_matrix(&A_local);
    print_matrix(&B_local);
    print_matrix(&C_local);
    /* */
  gettimeofday(&t1, NULL);
  gather_all(&A_local, &B_local, &C_local, block_size, A, B, C, counts, disps, blocktype, final);
  gettimeofday(&t2, NULL);
  time_difference = (t2.tv_sec-t1.tv_sec)*MO +(t2.tv_usec-t1.tv_usec);

  if (myrank == 0)
    fprintf(stdout, "%lu ", time_difference);
  /*if (myrank  == 0){
    print_matrix(A);
    print_matrix(B);
    print_matrix(C);
    }*/
  
  if (myrank == 0)
    fprintf(stdout, "%lf\n", MO*(MPI_Wtime() - initialTime));
 
  MPI_Finalize();
    
  free(A_local.data);
  free(B_local.data);
  free(C_local.data);

  free_matrix(A);
  free_matrix(B);
  free_matrix(C);
    
  return 0;
}
