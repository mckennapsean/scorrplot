
#ifndef SHMGYROSCOPE_H
#define SHMGYROSCOPE_H

#ifndef NULL
#define NULL 0
#endif



#include <stdio.h>
#include <semaphore.h>
#include <sys/mman.h>
#include <errno.h>
#include <fcntl.h> 
#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <list>
#include <string.h>





extern "C" {

//mutex for shared memory access
#define	FILE_MODE (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)
#define SEM_NAME "sem_gyroscope"
#define SHM_NAME_COR_INDEX "shm_gyroscope_cor_index"
//#define SHM_NAME_COR_VALUE "shm_gyroscope_cor_value"
sem_t *mutex;

//status indicator for shared memory
int *shm_request;
#define SHM_REQUEST_NONE 0
#define SHM_REQUEST_PASS_NAME 1
#define SHM_REQUEST_PASS_INDEX 2
#define SHM_REQUEST_GET_NAME 3
#define SHM_REQUEST_GET_INDEX 4
#define SHM_REQUEST_GET_COR_INDEX 5
#define SHM_REQUEST_GET_ACOR_INDEX 6
#define SHM_REQUEST_GET_COR_VALUE 7
#define SHM_REQUEST_CLOSE 8
#define SHM_REQUEST_PASS_PRIMARY 9
#define SHM_REQUEST_PASS_SECONDARY 10
#define SHM_REQUEST_GET_PRIMARY 11
#define SHM_REQUEST_GET_SECONDARY 12
#define SHM_REQUEST_GET_HIGHLIGHTED 13
#define SHM_REQUEST_GET_DENSITY 14
#define SHM_REQUEST_GET_SIZE 15

// shared memory for size of the dataset
int *shm_size;
//shared memory for passing current data point index
int *shm_index;
//shared memory for passing data point name
char *shm_name;
//shared memory for getting top/bottom n correlation of current projection
int *shm_cor_n;
int *shm_cor_indices;
double *shm_cor_value;
double *shm_vector;
double *shm_density;


void createShmCorIndex(int n){
  int fd = open(SHM_NAME_COR_INDEX, O_RDWR | O_CREAT, FILE_MODE);
  int *neg = new int[n];
  for(int i=0; i<n; i++){ neg[i] = -1; } 
  write(fd, neg, sizeof(int)*n);
  shm_cor_indices =(int*) mmap(NULL, sizeof(int)*n, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if(shm_cor_indices == MAP_FAILED){
    std::cout << "mmap index failed: " << errno << std::endl;
  }
  close(fd);
};

void closeShmCorIndex(int n){
  int res = munmap(shm_cor_indices, sizeof(int)*n);
  if(res != 0){
    std::cout << "munmap failed: " << errno << std::endl;
  }
  remove(SHM_NAME_COR_INDEX); 
};







}//end extern C

#endif
