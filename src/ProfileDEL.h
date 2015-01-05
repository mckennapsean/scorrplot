#ifndef ProfileDEL_H
#define ProfileDEL_H

#include "Display.h"
#include "DenseVector.h"
#include "DenseMatrix.h"
#include "Rotation.h"
#include "Animator.h"
#include "Font.h"

#include <sstream>
#include <iomanip>

template<typename TPrecision>
class ProfileDEL : public DisplayElement{

  private:


#define BUFSIZE 512



    int xM, yM;
    int pickH, pickW;  

    FortranLinalg::DenseVector<TPrecision> v;
    Data<TPrecision> &data;
    
    Rotation<TPrecision> Rone;
    Rotation<TPrecision> Rtwo;
    Animator &animator;

    Font &font;
   
    int selected;
    FortranLinalg::DenseVector<TPrecision> coloring;
    FortranLinalg::DenseVector<TPrecision> previous;
    FortranLinalg::DenseVector<TPrecision> agg1;
    FortranLinalg::DenseVector<TPrecision> agg2;
    int colormapping;
    bool colorOnly;
    std::set<int> aggIndex1;
    std::set<int> aggIndex2;

    FortranLinalg::DenseMatrix<TPrecision> A;
    FortranLinalg::DenseMatrix<TPrecision> b;

  public:

    ProfileDEL(Data<TPrecision> &d, Animator &a, Font &f) : data(d), Rone(d,
Rotation<TPrecision>::Primary), Rtwo(d, Rotation<TPrecision>::Secondary), animator(a), font(f){ 
      using namespace FortranLinalg;
   
      pickW = 7;
      pickH = 7;

      selected = -1;

      coloring = DenseVector<TPrecision>(d.raw.N());
      previous = DenseVector<TPrecision>(d.raw.N());
      colormapping = -1;
      colorOnly = false;
      agg1 = DenseVector<TPrecision>(d.raw.N());
      agg2 = DenseVector<TPrecision>(d.raw.N());
      Linalg<TPrecision>::Zero(agg1);
      Linalg<TPrecision>::Zero(agg2);
      A = DenseMatrix<TPrecision>(d.raw.N(), d.raw.M()+1);
      b = DenseMatrix<TPrecision>(d.raw.N(), 1);
    };


    ~ProfileDEL(){
      A.deallocate();
      b.deallocate();
      agg1.deallocate();
      agg2.deallocate();
      coloring.deallocate();
      previous.deallocate();
      v.deallocate();
    };

    void location(int xPos, int yPos, int w, int h){
      width = w;
      height = h;
      xLeft = xPos;
      yTop = yPos;
    };


    void init(){  

    };


    void display(void){
      using namespace FortranLinalg;

      glMatrixMode(GL_MODELVIEW); 	
      glLoadIdentity();
     
      glLineWidth(2);
     
      //draw coordinates
      glPushName(-1); 
      int n= data.raw.M();
      glColor4f(0.5, 0.5, 0.5, 1);
      glBegin(GL_LINES);
      for(int i=0; i<n; i++){
        double xc = xLeft + 10+ (i*(width-20))/(n-1.0);
        glLoadName(i);
        glVertex2f(xc, yTop);
        glVertex2f(xc, yTop+height);
      }
      glEnd();
      glLoadName(-1);
  
      //draw magnitude of point
      font.setSize(8);
      glColor4f(0.7, 0.7, 0.7, 1);
      for(int i=0; i<n; i++){
        double xc = xLeft + 10+ (i*(width-20))/(n-1.0);
        double yc = yTop+height;
        font.renderString(data.colname(i), xc-2, yc, 0, true);
      }
      

      //draw selected coloring aggregates
      glBegin(GL_LINES);
      glColor4f(0.420, 0.682, 0.839, 0.5);
      for(std::set<int>::iterator it=aggIndex1.begin(); it != aggIndex1.end(); ++it){
        double xc = xLeft + 10+ ((*it)*(width-20))/(n-1.0);
        glVertex2f(xc, yTop);
        glVertex2f(xc, yTop+height);
        
      }
      glColor4f(0.992, 0.553, 0.235, 0.5);
      for(std::set<int>::iterator it=aggIndex2.begin(); it != aggIndex2.end(); ++it){
        double xc = xLeft + 10+ ((*it)*(width-20))/(n-1.0);
        glVertex2f(xc, yTop);
        glVertex2f(xc, yTop+height);
      }
      glEnd();


      //draw selected point
      int index = data.getSelected();
      if(index >= 0){

       TPrecision xMin = Linalg<TPrecision>::MinColumn(data.raw, index);
       TPrecision xMax = Linalg<TPrecision>::MaxColumn(data.raw, index);
       // TPrecision xD = xMax-xMin;

       glColor4f(0.5, 0.5, 0.5, 1);
       font.setSize(9);
    
       std::stringstream s1;
       s1 << std::setiosflags(std::ios::fixed) << std::setprecision(2);
       s1 << xMin;
       font.renderString(s1.str(), xLeft-40, yTop+height);
     
       std::stringstream s2;
       s2 << std::setiosflags(std::ios::fixed) << std::setprecision(2);
       s2 << xMax;
       font.renderString(s2.str(), xLeft-40, yTop+11);
       
       drawProfile(index);
     
     }

     
     //draw profiles of other selected points
     for(int i=0; i<data.showProfile.size(); i++){
        drawProfile(data.showProfile[i], 0.5);
     }

     if(isInside(xM, yM)){
        //draw pick rectangle
        glColor4f(0.75, 0.75, 0.75, 0.75);
        glLineWidth(2);
        double pw = pickW/2.0;
        double ph = pickH/2.0;
        glBegin(GL_LINE_LOOP);
        glVertex2f(xM-pw, yM-ph);
        glVertex2f(xM+pw, yM-ph);
        glVertex2f(xM+pw, yM+ph);
        glVertex2f(xM-pw, yM+ph);
        glEnd();
     }

     glPopName();
    };

// resize profiles to the window size (manual default values set)
void reshape(int w, int h){
  int startLocX = (int) (790. * (double) w / 1100. );
  int startLocY = (int) (10. * (double) h / 710. );
  int endLocX = (int) (280. * (double) w / 1100. );
  int endLocY = (int) (150. * (double) h / 710. );
  location(startLocX, startLocY, endLocX, endLocY);
}    

    void idle(){

    };


    void keyboard(unsigned char key, int x, int y){
       if(!isInside(x,y)){ return;}
       if(key=='1'){
         data.setColoring(agg1);
         data.colormapping = true;
       }
       else if(key=='2'){
         data.setColoring(agg2);
         data.colormapping = true;
       }
       else if(key=='3'){
         data.colormapping = false;
         colormapping = -1;
       }

       glutPostRedisplay();

       if(!isInside(x, y)){ return; };
       if(key =='r'){
         colorOnly = ! colorOnly;
       }
       else if(key == 'c'){
         FortranLinalg::Linalg<TPrecision>::Zero(agg1);
         aggIndex1.clear();
         FortranLinalg::Linalg<TPrecision>::Zero(agg2);
         aggIndex2.clear();
         data.colormapping = false;
         data.showProfile.clear();
         data.showRowname.clear();
       }

    };



    void special(int key, int x, int y){
    };




    void mouse(int button, int state, int x, int y){
      using namespace FortranLinalg;
      xM = x;
      yM = y;

      if(!isInside(x, y)){ return; };

      int mod = glutGetModifiers();
      if ( selected >= 0 && button == GLUT_LEFT_BUTTON ){
	  if(state == GLUT_DOWN){
             bool rotate=false;
             if(mod == GLUT_ACTIVE_CTRL){
                add(aggIndex2, agg2);
                if(aggIndex2.size() > 0){
                  rotate=true;
                  Linalg<TPrecision>::SetColumn(b, 0, agg2);
                }
                  
             }
             else{
                add(aggIndex1, agg1);
                if(aggIndex1.size() > 0){
                  rotate=true;
                  Linalg<TPrecision>::SetColumn(b, 0, agg1);
                }

             }
             data.colormapping = true;
             colormapping = -1; 
          
             if(rotate) {
               Linalg<TPrecision>::Transpose(data.centered, A);
               for(int i=0; i<A.M(); i++){
                 A(i, A.N()-1) = 1;
               }
               DenseMatrix<TPrecision> sol = Linalg<TPrecision>::LeastSquares(A, b);
	       if(mod == GLUT_ACTIVE_CTRL){
                 Rtwo.setTarget(sol, 0);
                 animator.setAnimation(&Rtwo);
               }else{
                 Rone.setTarget(sol, 0);
                 animator.setAnimation(&Rone);
               }
               sol.deallocate(); 
              
            }

          } 
      }
      selected=-1;
      glutPostRedisplay();
    };



    // catch mouse move events
    void motion(int x, int y){

    };




    void passive(int x, int y){
      xM = x;
      yM = y;

      if( !isInside(x, y) ){ 
        if(colormapping != -1){
          data.colormapping = colormapping;
          data.setColoring(previous);
          colormapping = -1;
        }
        return;
      };

      if(colormapping == -1){
        colormapping = data.colormapping;
        FortranLinalg::Linalg<TPrecision>::Copy(data.coloring, previous);
      }

      GLint vp[4];
      glGetIntegerv(GL_VIEWPORT, vp);
      GLuint selectBuf[BUFSIZE];
      glSelectBuffer(BUFSIZE, selectBuf);
      glRenderMode(GL_SELECT);
      glInitNames();

      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      gluPickMatrix(x, vp[3]-y, pickW, pickH, vp);
      glOrtho (0, vp[2], vp[3], 0, 0, 1);

      display();
      GLint hits = glRenderMode(GL_RENDER);
      selected = -1;
      for(int i=0; i<hits; i++){
        int tmp = selectBuf[i*4 + 3];
        if(tmp != -1){
          selected = tmp;
        }
      }

      if(selected >= 0){
        FortranLinalg::Linalg<TPrecision>::ExtractRow(data.centered, selected, coloring); 
        data.setColoring(coloring);
        data.colormapping = true;
      }
      else{
        data.colormapping = colormapping;
        data.setColoring(previous);
      }

      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glutPostRedisplay();

    };



private:
   

  void add(std::set<int> &aggIndex, FortranLinalg::DenseVector<TPrecision> &agg){
      using namespace FortranLinalg;
    std::set<int>::iterator it = aggIndex.find(selected);
    if( it == aggIndex.end()){
      Linalg<TPrecision>::Add(agg, coloring, agg);
      aggIndex.insert(selected);
    }
    else{
      Linalg<TPrecision>::Subtract(agg, coloring, agg);
      aggIndex.erase(it);
    }
    data.setColoring(agg);
  };


  void drawProfile(int index, double alpha=1){
      using namespace FortranLinalg;

     TPrecision xMin = Linalg<TPrecision>::MinColumn(data.raw, index);
     TPrecision xMax = Linalg<TPrecision>::MaxColumn(data.raw, index);
     TPrecision xD = xMax-xMin;
     int n= data.raw.M();

     glColor4f(data.labelR(index), data.labelG(index), 
	    data.labelB(index), alpha);

     glLineWidth(2);
     glBegin(GL_LINE_STRIP);
     for(int i=0; i<n; i++){
       double xc =  (i*(width-20))/(n-1.0)+xLeft+10;
       glVertex2f(xc, yTop + height - height*(data.raw(i, index)-xMin)/xD);
     }
     glEnd();
  };


};

#endif
