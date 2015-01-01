#ifndef GYROSCOPE_H
#define GYROSCOPE_H

#include "Display.h"
#include "Data.h"
#include "DenseVector.h"
#include "DenseMatrix.h"

#include "Animator.h"
#include "Rotation.h"
#include "DisplayElement.h"
#include "ProjectionDEL.h"
#include "DensityProjectionDEL.h"
//#include "GMRADensityProjectionDEL.h"
#include "PCADEL.h"
#include "ProfileDEL.h"
#include "TextDEL.h"
#include "SearchDEL.h"
#include "LabelsDEL.h"
#include "PatchDEL.h"

#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <list>

#include "Font.h"

bool whiteBackground = false;

template<typename TPrecision>
class Gyroscope : public Display{

  private:

#define BUFSIZE 512

    int width, height;
    int xM, yM;
    int pickH, pickW; 
    Data<TPrecision> &data; 


    Animator animator;
    std::list<DisplayElement *> elements;

    Font &font;

    Rotation<TPrecision> Rone;
    Rotation<TPrecision> Rtwo;


  public:
    
    // density projection (for grabbing density values)
    DensityProjectionDEL<TPrecision> *dpd;

    Gyroscope(Font &f, Data<TPrecision> &d) : font(f), data(d), Rone(d, Rotation<TPrecision>::Primary),
    Rtwo(d, Rotation<TPrecision>::Secondary){ 
      pickH = 5;
      pickW = 5;
    };



    std::string title(){
      return "Gyroscope";
    };



    void reshape(int w, int h){
      width = w;
      height = h;

      int size = std::min(w, h)/2 - 5;
      glViewport(0, 0, w, h);       
      glMatrixMode(GL_PROJECTION);  
      glLoadIdentity();
      glOrtho (0, width, height, 0, 0, 1);

      // for each element in display (ignoring PCA graph), resize to the window size
      // please note that this is manually using pixel values defined in each class
      // not the best, but do-able until more elegant variables exist to store the values
      for(std::list<DisplayElement *>::iterator it=elements.begin(); 
          it != elements.end(); ++it){
        (*it)->reshape(w, h);
      }
    };



    void setupDisplay(bool useDensity = true, bool showProfile = true, bool
        showPatch = false){
      //setup display elements

      for(std::list<DisplayElement *>::iterator it = elements.begin(); it !=
          elements.end(); ++it){
        delete (*it);
      }

      elements.clear();

 
      if(useDensity){
        //GMRADensityProjectionDEL<TPrecision> *pd = new GMRADensityProjectionDEL<TPrecision>(data, animator, font);
        dpd = new DensityProjectionDEL<TPrecision>(data, animator, font);
        dpd->location(151, 10, 630, 630);
        elements.push_back(dpd);
      }
      else{
        ProjectionDEL<TPrecision> *pd = new ProjectionDEL<TPrecision>(data, animator, font);
        pd->location(151, 10, 630, 630);
        elements.push_back(pd);
      }

      int w=150, h=50;
      PCADEL<TPrecision> *b1 = new PCADEL<TPrecision>(-1, data, animator, 0.5, 0.5, 0.5);
      b1->location(10, 10, w, h); 
      elements.push_back(b1);
      if(data.lPCA.size() > 1){
        int nl=0;
        for(typename std::map<int, PCA<TPrecision>* >::iterator it = data.lPCA.begin();
            it!=data.lPCA.end(); ++it){
          int l = it->first;
          PCADEL<TPrecision> *b = new PCADEL<TPrecision>(l, data, animator, 
              data.colLabelR(l), data.colLabelG(l), data.colLabelB(l));
          b->location(10, 10+(nl+1)*(h+5), w, h); 
          elements.push_back(b);
          nl++;
        }
      }


      if(showProfile){
        ProfileDEL<TPrecision> *pr=new ProfileDEL<TPrecision>(data, animator, font);
        pr->location(790, 10, 280, 150);
        elements.push_back(pr);
      }

      TextDEL<TPrecision> *tx=new TextDEL<TPrecision>(data, animator, font);
      tx->location(790, 210, 230, 770);
      elements.push_back(tx);

      /*
         SearchDEL<TPrecision> *se=new SearchDEL<TPrecision>(data, animator, font);
         se->location(790, 170, 230, 50);
         elements.push_back(se);
       */

      
      LabelsDEL<TPrecision> *la=new LabelsDEL<TPrecision>(data, animator, font);
      la->location(151, 650, 630, 40);
      elements.push_back(la);

      // new patch-by-patch image viewer
      if(showPatch){
        PatchDEL<TPrecision> *pa = new PatchDEL<TPrecision>(data);
        pa->location(10, 300, 100, 100);
        elements.push_back(pa);
      }

    }


    void init(){
      init(true, true, false);
    };

    void init(bool useDensity, bool showProfile, bool showPatch) {  
      glClearColor(0.15, 0.15, 0.15, 0);
      glEnable(GL_POLYGON_SMOOTH);
      glEnable(GL_LINE_SMOOTH);
      glEnable( GL_POINT_SMOOTH );
      glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_BLEND);


      glutSetCursor(GLUT_CURSOR_CROSSHAIR);

      setupDisplay(useDensity, showProfile, showPatch); 

    };



    void printHelp(){

      std::cout << "Anywhere" << std::endl;
      std::cout << "------------------" << std::endl;
      std::cout << "Press q to quit and return order list of correlated genes" << std::endl  << std::endl;
      std::cout << "Press i to search for genes in the file search.txt in the current working directory" << std::endl  << std::endl;
      std::cout << "Press o to ouput the 20 most correlated genes in cor.txt and the 20 most anti-correlated genes in acor.txt in the current working directory" << std::endl  << std::endl;

      std::cout << "Projection Display" << std::endl;
      std::cout << "------------------" << std::endl;
      std::cout << "Click on point to select current set of genes as active (shown in text list when mouse pointer is not in projection display)" << std::endl;
      std::cout << "Shift click to add to profile display" << std::endl;
      std::cout << "Up and down arrow to increase / decrease point size" << std::endl ; 
      std::cout << "Left and right arrow to adjust transparency" << std::endl ; 
      std::cout << "Press t for computation of p-Values for current projection" << std::endl ;
      std::cout << "Press s for Spearman's correlation" << std::endl ;
      std::cout << "Press p for Pearson's correlation" << std::endl ; 
      std::cout << "Press c to clear searched" << std::endl ; 
      std::cout << "Press + to increase selection square" << std::endl ;
      std::cout << "Press - to decrease selection square" << std::endl << std::endl;
      std::cout << "Press > to increase bandwidth for density" << std::endl ;
      std::cout << "Press < to decrease bandwidth for density" << std::endl << std::endl;




      std::cout << "Profile Display" << std::endl;
      std::cout << "--------------" << std::endl;
      std::cout << "Hover and click (ctrl) click to color/align by coordinates" << std::endl;
      std::cout << "Press 1,2,3 to color by clicked, ctrl clicked or labels" << std::endl; 
      std::cout << "Press c to clear selected profiles" << std::endl  << std::endl;


      std::cout << "PCA Display" << std::endl;
      std::cout << "-----------" << std::endl;
      std::cout << "Click on bar to align primary direction to principal component" << std::endl;
      std::cout << "Ctrl click to align secondary direction" << std::endl  << std::endl;

      std::cout << "Label Display" << std::endl;
      std::cout << "-------------" << std::endl;
      std::cout << "Click to highlight label" << std::endl << std::endl;

      std::cout << "Text Display" << std::endl;
      std::cout << "-------------" << std::endl;
      std::cout << "Click to align first axis of projection" << std::endl << std::endl;
      std::cout << "Ctrl click to align second axis of projection" << std::endl << std::endl;
      std::cout << "Shift click to add to selected profiles" << std::endl << std::endl;
    };




    void display(void){

      glMatrixMode(GL_MODELVIEW); 	
      glLoadIdentity();

      glClear(GL_COLOR_BUFFER_BIT);

      // set background to white if user changes it
      if(whiteBackground)
        glClearColor(1.00, 1.00, 1.00, 0);
      else
        glClearColor(0.15, 0.15, 0.15, 0);

      for(std::list<DisplayElement *>::iterator it=elements.begin(); 
          it != elements.end(); ++it){
        (*it)->display();
      }

      glutSwapBuffers();

    };




    void idle(){
      using namespace FortranLinalg;
#ifdef SHMGYROSCOPE_H
      sem_wait(mutex);
      bool cor = (*shm_request) == SHM_REQUEST_GET_COR_INDEX; 
      switch( (*shm_request) ){
        case SHM_REQUEST_PASS_INDEX:
          if(*shm_index >=0 && *shm_index < data.raw.N()){
            data.showRowname.push_back(*shm_index);
            data.showProfile.push_back(*shm_index);
          }
          (*shm_request) = SHM_REQUEST_NONE;
          glutPostRedisplay();
          break;
        case SHM_REQUEST_PASS_NAME:
          data.addShowed( std::string(shm_name) );
          (*shm_request) = SHM_REQUEST_NONE;
          glutPostRedisplay();
          break;
        case SHM_REQUEST_GET_SIZE:
          *shm_size = data.raw.N();
          (*shm_request) = SHM_REQUEST_NONE;
          break;
        case SHM_REQUEST_GET_INDEX:
          *shm_index = data.getSelected();
          (*shm_request) = SHM_REQUEST_NONE;
          break;
        case SHM_REQUEST_GET_COR_INDEX:
        case SHM_REQUEST_GET_ACOR_INDEX:
          createShmCorIndex(*shm_cor_n);
          data.getCorIndicies(*shm_cor_n, shm_cor_indices, cor);
          (*shm_request) = SHM_REQUEST_NONE;
          break;
        case SHM_REQUEST_GET_NAME:
          if(*shm_index >=0 && *shm_index < data.raw.N()){
            strcpy(shm_name, data.rowname(*shm_index).c_str() );
          }
          else{
            strcpy(shm_name, "NA");
          }
          (*shm_request) = SHM_REQUEST_NONE;
          break;
        case SHM_REQUEST_GET_DENSITY:
          if(*shm_index >=0 && *shm_index < data.raw.N()){
            (*shm_density) = dpd->getDensityFromIndex(*shm_index);
          }else{
            (*shm_density) = -1;
          }
          (*shm_request) = SHM_REQUEST_NONE;
          break;
        case SHM_REQUEST_GET_COR_VALUE:
          if(*shm_index >=0 && *shm_index < data.raw.N()){
            (*shm_cor_value) = data.P(0, *shm_index);
          }
          else{
            (*shm_cor_value) = -2;
          }
          (*shm_request) = SHM_REQUEST_NONE;
          break;
        case SHM_REQUEST_PASS_PRIMARY:
          {
            DenseVector<TPrecision> v(data.raw.M(), shm_vector);
            Rone.setTarget(v);
            animator.setAnimation(&Rone);
            (*shm_request) = SHM_REQUEST_NONE;
          }
          break;
        case SHM_REQUEST_PASS_SECONDARY:
          {
            DenseVector<TPrecision> v(data.raw.M(), shm_vector);
            Rtwo.setTarget(v);
            animator.setAnimation(&Rtwo);
            (*shm_request) = SHM_REQUEST_NONE;
          }
          break;
        case SHM_REQUEST_GET_PRIMARY:
          for(int i=0; i<data.V.M(); i++){
            shm_vector[i] = data.V(i, 0);
          }
          (*shm_request) = SHM_REQUEST_NONE;
          break;
        case SHM_REQUEST_GET_SECONDARY:
          for(int i=0; i<data.V.M(); i++){
            shm_vector[i] = data.V(i, 1);
          }
          (*shm_request) = SHM_REQUEST_NONE;
          break;
        case SHM_REQUEST_CLOSE:
          (*shm_request) = SHM_REQUEST_NONE;
          exit(0);
      }
      sem_post(mutex);
#endif

      animator.step();
    };





    void keyboard(unsigned char key, int x, int y){
      if(animator.isRunning()){
        return;
      }      
      switch(key){

        case 'h':
        case 'H':
          printHelp();
          break;

          // toggle background color between white and black
        case 'b':
        case 'B':
          if(whiteBackground){
            whiteBackground = false;
          }
          else{
            whiteBackground = true;
          }
          break;
      }
      for(std::list<DisplayElement *>::iterator it=elements.begin(); 
          it != elements.end(); ++it){
        (*it)->keyboard(key, x, y);
      }
      glutPostRedisplay();
    };



    void special(int key, int x, int y){
      if(animator.isRunning()){
        return;
      }

      for(std::list<DisplayElement *>::iterator it=elements.begin(); 
          it != elements.end(); ++it){
        (*it)->special(key, x, y);
      }
      glutPostRedisplay();
    };




    void mouse(int button, int state, int x, int y){
      if(animator.isRunning()){
        return;
      }
      xM = x;
      yM = y;
      for(std::list<DisplayElement *>::iterator it=elements.begin(); 
          it != elements.end(); ++it){
        (*it)->mouse(button, state, x, y);
      }
      glutPostRedisplay();

    };



    void motion(int x, int y){
      if(animator.isRunning()){
        return;
      }
      for(std::list<DisplayElement *>::iterator it=elements.begin(); 
          it != elements.end(); ++it){
        (*it)->motion(x, y);
      }
      glutPostRedisplay();
    };




    void passive(int x, int y){
      if(animator.isRunning()){
        return;
      }
      xM = x;
      yM = y;

      for(std::list<DisplayElement *>::iterator it=elements.begin(); 
          it != elements.end(); ++it){
        (*it)->passive(x, y);
      }
      glutPostRedisplay();

    };









};

#endif
