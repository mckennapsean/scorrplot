#ifndef NULL
#define NULL 0
#endif

#define R_NO_REMAP


#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <list>
#include <string.h>



//needs to be included before other scorr stuff
#include "Shmscorr.h"

#include "GlutStrokeFont.h"
#include "sCorr.h"
#include "Data.h"
#include "DenseVector.h"
#include "DenseMatrix.h"


extern "C" {

//dimensionality
int m;


//shm message passing R <-> scorr


SEXP scorrHighlightIndex(SEXP Rindices, SEXP Rn){
  int *indices = INTEGER(Rindices);
  int n = *INTEGER(Rn);

  for(int i=0; i<n; i++){
    sem_wait(mutex);
    *shm_index = indices[i]-1;
    *shm_request = SHM_REQUEST_PASS_INDEX;
    sem_post(mutex);
  
    while(*shm_request != 0){
      usleep(10);
    }
  }
  
  return R_NilValue;
};


SEXP scorrHighlightName(SEXP Rnames, SEXP Rn){
  int n = *INTEGER(Rn);

  for(int i=0; i<n; i++){
    const char *name = CHAR(STRING_ELT(Rnames,i));
    sem_wait(mutex);
    strcpy(shm_name, name);
    *shm_request = SHM_REQUEST_PASS_NAME;
    sem_post(mutex);
  
    while(*shm_request != 0){
      usleep(10);
    }
  }
  
  return R_NilValue;
};




SEXP scorrGetSize(){
   
   sem_wait(mutex);
   *shm_request = SHM_REQUEST_GET_SIZE;
   sem_post(mutex);
   
   while(*shm_request != 0){
     usleep(10);
   }
   
   SEXP Rsize;
   PROTECT(Rsize = Rf_allocVector(INTSXP, 1));
   *INTEGER(Rsize) = *shm_size;
   UNPROTECT(1);
   
   return Rsize;
};




SEXP scorrGetSelected(){

   sem_wait(mutex);
   *shm_request = SHM_REQUEST_GET_INDEX;
   sem_post(mutex);
  
   while(*shm_request != 0){
     usleep(10);
   }

   SEXP Rindex;
   PROTECT(Rindex = Rf_allocVector(INTSXP, 1));
   *INTEGER(Rindex) = *shm_index + 1;
   UNPROTECT(1);
   return Rindex;
};


SEXP scorrClose(){

   sem_wait(mutex);
   *shm_request = SHM_REQUEST_CLOSE;
   sem_post(mutex);
  
   while(*shm_request != 0){
     usleep(10);
   }
   return R_NilValue;
};


SEXP scorrGetCorIndex(SEXP Rn){

   int n = *INTEGER(Rn);
   bool top = n>0;
   n=abs(n);

   
   sem_wait(mutex);
   *shm_cor_n = n;
   if(top){
     *shm_request = SHM_REQUEST_GET_COR_INDEX;
   }
   else{
     *shm_request = SHM_REQUEST_GET_ACOR_INDEX;
   }
   createShmCorIndex(n);
   sem_post(mutex);
  
   while(*shm_request != 0){
     usleep(10);
   }
     
   SEXP Rindices;
   PROTECT( Rindices= Rf_allocVector(INTSXP, n) );
   for(int i=0; i<n; i++){
     INTEGER(Rindices)[i] = shm_cor_indices[i] + 1;
   }
   UNPROTECT(1);

   closeShmCorIndex(n);
   return Rindices;
};






SEXP scorrGetName(SEXP Rindices, SEXP Rn){
  int *indices = INTEGER(Rindices);
  int n = *INTEGER(Rn);

  SEXP list;
  PROTECT( list = Rf_allocVector(STRSXP, n));

  for(int i=0; i<n; i++){
    sem_wait(mutex);
    *shm_index = indices[i]-1;
    *shm_request = SHM_REQUEST_GET_NAME;
    sem_post(mutex);
  
    while(*shm_request != 0){
      usleep(2);
    }
    SET_STRING_ELT(list, i, Rf_mkChar(shm_name) );
    
  }
  UNPROTECT(1);
  
  return list;
};




SEXP scorrGetDensity(SEXP Rindices, SEXP Rn){
  int *indices = INTEGER(Rindices);
  int n = *INTEGER(Rn);

  SEXP list;
  PROTECT( list = Rf_allocVector(REALSXP, n));

  for(int i=0; i<n; i++){
    sem_wait(mutex);
    *shm_index = indices[i]-1;
    *shm_request = SHM_REQUEST_GET_DENSITY;
    sem_post(mutex);
  
    while(*shm_request != 0){
      usleep(2);
    }
    REAL(list)[i] = *shm_density;
    
  }
  UNPROTECT(1);
  
  return list;
};




SEXP scorrGetCorValue(SEXP Rindices, SEXP Rn){
  int *indices = INTEGER(Rindices);
  int n = *INTEGER(Rn);

  SEXP list;
  PROTECT( list = Rf_allocVector(REALSXP, n));

  for(int i=0; i<n; i++){
    sem_wait(mutex);
    *shm_index = indices[i]-1;
    *shm_request = SHM_REQUEST_GET_COR_VALUE;
    sem_post(mutex);
  
    while(*shm_request != 0){
      usleep(2);
    }
    REAL(list)[i] = *shm_cor_value;
    
  }
  UNPROTECT(1);
  
  return list;
};


SEXP scorrSetProjection(SEXP Rv, SEXP Rmv, SEXP Rtype){
  double *v = REAL(Rv);
  int type = *INTEGER(Rtype);
  int mv = *INTEGER(Rmv); 
  if(mv != m){
    return R_NilValue;
  }


  sem_wait(mutex);
  for(int i=0; i<m; i++){
    shm_vector[i] = v[i];
  }
  if(type == 0){
    *shm_request = SHM_REQUEST_PASS_PRIMARY;
  }
  else{
    *shm_request = SHM_REQUEST_PASS_SECONDARY;
  }
  sem_post(mutex);
  
  while(*shm_request != 0){
    usleep(2);
  }
  


  return R_NilValue;
};



SEXP scorrGetProjection(SEXP Rtype){
  int type = *INTEGER(Rtype);


  sem_wait(mutex);
  if(type == 0){
    *shm_request = SHM_REQUEST_GET_PRIMARY;
  }
  else{
    *shm_request = SHM_REQUEST_GET_SECONDARY;
  }
  sem_post(mutex);
  
  while(*shm_request != 0){
    usleep(2);
  }
  
  SEXP Rv;
  PROTECT( Rv = Rf_allocVector(REALSXP, m));
  double *v=REAL(Rv);
  for(int i=0; i<m; i++){
    v[i] = shm_vector[i];
  }
  UNPROTECT(1);

  return Rv;
};
//rendering scorr methods


//#define MAKE_STRING_(x) #x
//#define MAKE_STRING(x) MAKE_STRING_(x)


int mainWindow;
sCorr<double> *D_projection = NULL;


void display1(void){
  D_projection->display();
};

void idle1(void){
  D_projection->idle();
};

void mouse1(int button, int state, int x, int y){
  D_projection->mouse(button, state, x, y);
}

void motion1(int x, int y){
  D_projection->motion(x, y);
}

void passive1(int x, int y){
  D_projection->passive(x, y);
}

void keyboard1(unsigned char key, int x, int y){
  switch(key){
        case 27:   
	case 'q':
	case 'Q':
          exit(0);
          break;
  }
  D_projection->keyboard(key, x, y);
}

void special1(int key, int x, int y){
  D_projection->special(key, x, y);
}

void reshape1(int w, int h){
  D_projection->reshape(w, h);
}



void visible(int vis)
{
    if (vis == GLUT_VISIBLE)
        glutIdleFunc(idle1);
    else
        glutIdleFunc(NULL);
}




void printHelp(){
  D_projection->printHelp();	
}


SEXP scorr(SEXP Rm, SEXP Rn, SEXP Rx, SEXP Rl, SEXP Rln, SEXP Rperms, SEXP
    Rcolors, SEXP RuseDensity, SEXP RshowProfile, SEXP RshowPatch) {
      using namespace FortranLinalg;
  m = *INTEGER(Rm);

  //create shared memory and mutex
  shm_request = (int *)  mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);
  if(shm_request == MAP_FAILED){
    std::cout << "mmap index failed: " << errno << std::endl;
  }
  (*shm_request) = 0;
  
  shm_size = (int *)  mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);
  if(shm_size == MAP_FAILED){
    std::cout << "mmap index failed: " << errno << std::endl;
  }
  (*shm_size) = 0;

  shm_index = (int *)  mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);
  if(shm_index == MAP_FAILED){
    std::cout << "mmap index failed: " << errno << std::endl;
  }
  (*shm_index) = 0;

  shm_name = (char *) mmap(NULL, sizeof(char)*10000, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);
  if(shm_name == MAP_FAILED){
    std::cout << "mmap name failed: " << errno << std::endl;
  }

  shm_cor_n = (int *)  mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);
  if(shm_cor_n == MAP_FAILED){
    std::cout << "mmap index failed: " << errno << std::endl;
  }
  (*shm_cor_n) = 0;
  
  shm_density = (double *)  mmap(NULL, sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);
  if(shm_density == MAP_FAILED){
    std::cout << "mmap value failed: " << errno << std::endl;
  }
  (*shm_density) = -1;

  shm_cor_value = (double *)  mmap(NULL, sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);
  if(shm_cor_value == MAP_FAILED){
    std::cout << "mmap value failed: " << errno << std::endl;
  }
  (*shm_cor_value) = -2;

  shm_vector = (double *)  mmap(NULL, sizeof(double)*m, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);
  if(shm_vector == MAP_FAILED){
    std::cout << "mmap vector failed: " << errno << std::endl;
  }
  (*shm_vector) = -2;



  mutex = sem_open(SEM_NAME, O_CREAT | O_EXCL, FILE_MODE, 1);
  if(mutex == SEM_FAILED){
    std::cout << "sem_open failed: " << errno << std::endl;
  }
  if(sem_unlink(SEM_NAME) != 0 ){
    std::cout << "sem_unlink failed: " << errno << std::endl;
  }
   
  //create scorr process 
  pid_t pID = fork();
  if(pID == 0){//child

    SEXP rl, cl;
    const char *rn, *cn;
    Rf_GetMatrixDimnames(Rx, &rl, &cl, &rn, &cn);
    double *x = REAL(Rx);
    double *cols = REAL(Rcolors);
    int *l = INTEGER(Rl);
    int n = *INTEGER(Rn);
    int nPerms = *INTEGER(Rperms);
 
    bool useDensity = *INTEGER(RuseDensity);
    bool showProfile = *INTEGER(RshowProfile);
    bool showPatch = *INTEGER(RshowPatch);

    std::vector<std::string> cnames;
    for(int i=0; i<m; i++){
      const char *name = CHAR(STRING_ELT(rl,i));
      cnames.push_back(name);
    }

    std::vector<std::string> rnames;
    for(int i=0; i<n; i++){
      const char *name = CHAR(STRING_ELT(cl,i));
      rnames.push_back(name);
    }

    DenseMatrix<double> X(m, n, x);
    X = Linalg<double>::Copy(X);
    DenseVector<int> L(n, l);
    L = Linalg<int>::Copy(L);
    DenseVector<double> colors(n, cols);


    int ml = Linalg<int>::Max(L);
    std::vector<std::string> lnames;
    for(int i=0; i<=ml; i++){
      const char *name = CHAR(STRING_ELT(Rln,i));
      lnames.push_back(name);
    }

    GlutStrokeFont font;
    Data<double> data(X, L, rnames, cnames, lnames, nPerms);
    data.setColoring(colors);

    int argc = 1;
    char *argv = "scorr";

    D_projection = new sCorr<double>(font, data);
      
    glutInit(&argc, &argv);
    glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE);
    glutInitWindowSize(1100, 710);

  
 
    mainWindow = glutCreateWindow(D_projection->title().c_str());
    

    glutDisplayFunc(display1);
    glutVisibilityFunc(visible);
    glutReshapeFunc(reshape1);
    glutMouseFunc(mouse1);
    glutMotionFunc(motion1);
    glutPassiveMotionFunc(passive1);
    glutKeyboardFunc(keyboard1);
    glutSpecialFunc(special1);

    D_projection->init(useDensity, showProfile, showPatch);

    glutShowWindow();
    glutMainLoop(); 

  }

  return R_NilValue;
};


}//end extern C
