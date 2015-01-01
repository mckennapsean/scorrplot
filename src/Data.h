#ifndef DATA_H
#define DATA_H

#include "DenseMatrix.h"
#include "Linalg.h"
#include "PCA.h"
#include "colormapper.h"
//#include "IO.h"
#include "Geometry.h"
#include "SquaredEuclideanMetric.h"

#include <vector>
#include <map>
#include <set>
#include <algorithm>


template<typename TPrecision>
class Pair{
  public:
    int index;
    TPrecision value;

    Pair(){};
    Pair(int i, TPrecision v):index(i), value(v){};

    friend bool operator< (const Pair<TPrecision> &p1, const Pair<TPrecision> &p2){
      return p1.value < p2.value;
    };
};




template<typename TPrecision>
class Data{

 

  private:
      //NULL-distribution (permutation testing)
      FortranLinalg::DenseMatrix<TPrecision> pP0; //(positive correlations)
      FortranLinalg::DenseMatrix<TPrecision> pN0; //(negative correlations)
  
  public:
    int nPerms;
 
    //unit vectors for correlation sphere
    //DenseMatrix<TPrecision> units; 

    //raw data
    FortranLinalg::DenseMatrix<TPrecision> raw;
    //centered/ranked data
    FortranLinalg::DenseMatrix<TPrecision> centered;
    //label id for each data point
    FortranLinalg::DenseVector<int> labels;

    //coloring values
    FortranLinalg::DenseVector<TPrecision> coloring;

    //names of rows and columns
    std::vector<std::string> rownames;
    std::vector<std::string> colnames;
    //names for label ids
    std::vector<std::string> labelnames;
    
    //Spearman correlation coefficients
    bool useSpearman;

    
    //PCA of complete data set 
    PCA<TPrecision> *pca; 
    //PCA for each label group
    std::map<int, PCA<TPrecision>* > lPCA;
    



    //Projection matrix
    FortranLinalg::DenseMatrix<TPrecision> V;
    //store secondary directon (primary is in first column)
    FortranLinalg::DenseVector<TPrecision> secondary;

    //Projected data
    FortranLinalg::DenseMatrix<TPrecision> P;
    //Projected unit vectros
    //DenseMatrix<TPrecision> PU;


    //p-values for current projection    
    FortranLinalg::DenseVector<TPrecision> pValsP;
    FortranLinalg::DenseVector<TPrecision> pValsN;


    
    //highlight for labels
    int labelHighlight;
    //list of label ids to draw in foreground
    FortranLinalg::DenseVector<int> labelSelected;



    //Selection of points under cursor (clicked or currently under mouse
    //pointer)
    std::vector<int> selectedIndicies;
    //selected point: index into of above list
    int selectedIndex;


    //list of points for which names are shown
    std::vector<int> showRowname;

    //list of points for parallel coordinates plot
    std::vector<int> showProfile;

   
    bool colormapping;
    ColorMapper<float> cmap;


    Data( FortranLinalg::DenseMatrix<TPrecision> &r, FortranLinalg::DenseVector<int> l,
          std::vector<std::string> &rn, std::vector<std::string> &cn, 
          std::vector<std::string> &ln, int perms =0 ){
      using namespace FortranLinalg;
      pca = NULL;
      useSpearman = false;

      pValsP = DenseVector<TPrecision>(10);
      pValsN = DenseVector<TPrecision>(10);
      Linalg<TPrecision>::Set(pValsP, -1); 
      Linalg<TPrecision>::Set(pValsN, -1);
      nPerms = perms;
      pP0 = DenseMatrix<TPrecision>(10, nPerms);
      pN0 = DenseMatrix<TPrecision>(10, nPerms);
    
      setData(r, l, rn, cn, ln);  
      labelHighlight = -1;
      selectedIndex = -1;
     /* if(!IO<TPrecision>::fileExists(SEARCHFILE ) ){
        std::ofstream file;
        file.open(SEARCHFILE);
        file.close();
      }*/
    };


    void cleanup(){
      //units.deallocate();
      raw.deallocate();
      centered.deallocate();
      labels.deallocate();
      V.deallocate();
      secondary.deallocate();
      P.deallocate();
      //PU.deallocate();
      labelSelected.deallocate();
      
      if(pca != NULL){
        pca->cleanup();
        delete pca;
        pca = NULL;
      }
      for(typename std::map<int, PCA<TPrecision>* >::iterator it = lPCA.begin();
          it!=lPCA.end(); ++it){
        it->second->cleanup();
      } 
      lPCA.clear();

      
      pP0.deallocate();
      pN0.deallocate();
      pValsP.deallocate();
      pValsN.deallocate();
     

    };



    void setSpearman(bool spearman){
       useSpearman = spearman;
       updateCentered();
       updateProjection();
    };



    void updateCentered(){
      using namespace FortranLinalg;
      
      centered.deallocate();
      centered = Linalg<TPrecision>::Copy(raw);
      if(useSpearman){
        //rank each row
        std::vector< Pair<TPrecision> > pairs(centered.M());

        for(int i=0; i<centered.N(); i++){

          for(int j=0; j<centered.M(); j++){
            pairs[j].index = j;
            pairs[j].value = centered(j, i);
          }
 
          std::sort(pairs.begin(), pairs.end());
          for(int j=0; j<centered.M(); j++){
            centered(pairs[j].index, i) = j;
          }

        }
      }
 


      DenseVector<TPrecision> means = Linalg<TPrecision>::SumRows(centered);
      Linalg<TPrecision>::Scale(means, 1.0f/raw.M(), means);
      Linalg<TPrecision>::SubtractRowwise(centered, means, centered);
      means.deallocate();   

      for(int i=0; i<centered.N(); i++){
        TPrecision l = Linalg<TPrecision>::LengthColumn(centered, i);
        Linalg<TPrecision>::ScaleColumn(centered, i, 1.0f/l);
      } 


      //per label pca
      if(pca != NULL){
        pca->cleanup();
        delete pca;
        pca = NULL;
      }
      for(typename std::map<int, PCA<TPrecision>* >::iterator it = lPCA.begin();
          it!=lPCA.end(); ++it){
        it->second->data.deallocate();
        it->second->cleanup();
  
      } 
      lPCA.clear();


      pca = new PCA<TPrecision>(centered, centered.M()-1, false);



      std::set<int> nlabels;
      for(int i=0; i<labels.N(); i++){
        nlabels.insert(labels(i));
      };

      labelSelected.deallocate();
      labelSelected = DenseVector<int>(nlabels.size());
      Linalg<int>::Zero(labelSelected);

      for(std::set<int>::iterator it=nlabels.begin(); it != nlabels.end();
          ++it){
        int lID = *it;
        std::vector<int> index;
        for(int i=0; i<labels.N(); i++){
          if(labels(i) == lID){
            index.push_back(i); 
          } 
        }
        DenseMatrix<TPrecision> tmp(raw.M(), index.size());
        for(int i=0; i<index.size(); i++){
          Linalg<TPrecision>::SetColumn(tmp, i, centered, index[i]); 
        }
        lPCA[lID] = new PCA<TPrecision>(tmp, tmp.M()-1, false);
      }
 
    };


    void setData(FortranLinalg::DenseMatrix<TPrecision> &r, FortranLinalg::DenseVector<int> l,
          std::vector<std::string> &rn, std::vector<std::string> &cn, 
          std::vector<std::string> &ln){
      using namespace FortranLinalg;

      //cleanup();

      raw = r;
      labels=l;
      rownames = rn;
      colnames = cn;
      labelnames = ln;

      /*
      units = DenseMatrix<TPrecision>(raw.M(), raw.M());
      double ul = sqrt(1.0/raw.M());
      Linalg<TPrecision>::Set(units, -ul);
      for(int i=0; i<units.N(); i++){
        units(i,i) = 1.0 - ul;
      };      
 */

      updateCentered();
      computeNull(); 
 

      V = DenseMatrix<TPrecision>(raw.M(), 2);
      secondary = Linalg<TPrecision>::ExtractColumn(pca->ev, pca->ev.N()-2);
      Linalg<TPrecision>::SetColumn(V, 0, pca->ev, pca->ev.N()-1 );
      Linalg<TPrecision>::SetColumn(V, 1, pca->ev, pca->ev.N()-2);
      
      P = project(centered);
      //PU = project(units);
      pValues();


      selectedIndex = 0;
      colormapping = false;
      cmap.set(0.1, 0.996, 0.941, 0.929, 0.698, 0.231, 0.627, 0.298, 0.125);
      coloring = Linalg<TPrecision>::ExtractRow(raw, 0);
      setColoring(coloring);

    };






    void setColoring(FortranLinalg::DenseVector<TPrecision> &v){
      using namespace FortranLinalg;
       Linalg<TPrecision>::Copy(v, coloring);
       TPrecision mi = Linalg<TPrecision>::Min(coloring);
       TPrecision ma = Linalg<TPrecision>::Max(coloring);
       Linalg<TPrecision>::Subtract(coloring, mi, coloring);
       Linalg<TPrecision>::Scale(coloring, 1.0/(ma-mi), coloring);
    };


    void clearSelection(){
      selectedIndicies.clear();
    };

    void addSelected(int index){
      selectedIndicies.push_back(index);
    };

    void addSelected(std::vector<int> indices){
      for(int i=0; i<indices.size(); i++){
        addSelected(indices[i]);
      }
    };


   
    int getSelected(){
      int size = selectedIndicies.size();
      if(size == 0){
        return -1;
      }
      
      if(selectedIndex < 0){
        selectedIndex = 0;
      }
      else if(selectedIndex >= size){
        selectedIndex = size-1;
      }
      return selectedIndicies[selectedIndex];
    };
   




    void setProjection(FortranLinalg::DenseVector<TPrecision> &v, int index){
      using namespace FortranLinalg;
      Linalg<TPrecision>::SetColumn(V, index, v);
    };



    void updateProjection(){
      project(P, centered);
      //project(PU, units);
      pValues();
    };

   FortranLinalg:: DenseMatrix<TPrecision> project(FortranLinalg::DenseMatrix<TPrecision> &X){
      using namespace FortranLinalg;
      DenseMatrix<TPrecision> pr(2, X.N());
      project(pr, X);
      return pr;

    };


    void project(FortranLinalg::DenseMatrix<TPrecision> &pr, FortranLinalg::DenseMatrix<TPrecision> &X){
      using namespace FortranLinalg;
      Linalg<TPrecision>::Multiply(V, X, pr, true);
      Linalg<TPrecision>::Set(pValsP, -1); 
      Linalg<TPrecision>::Set(pValsN, -1);
    };



    //Compute p values /  permutation test
    void pValues(){
      using namespace FortranLinalg;
 
      //empirical
      DenseVector<TPrecision> pP(10);
      DenseVector<TPrecision> pN(10);
        
      Linalg<TPrecision>::Zero(pValsP);
      Linalg<TPrecision>::Zero(pValsN);
        
      Linalg<TPrecision>::Zero(pP);
      Linalg<TPrecision>::Zero(pN);

      for(int k=0; k < P.N(); k++){
        double dot = P(0, k);
        for(int j=0; j<10; j++){
          if(dot >= j/10.0){
            pP(j)++;
          }
          else if(dot <= -j/10.0){
            pN(j)++;
          } 
        }
      }      
      Linalg<TPrecision>::Scale(pP, 1.0/P.N(), pP);
      Linalg<TPrecision>::Scale(pN, 1.0/P.N(), pN);

      

        

      for(int k=0; k<nPerms; k++){
        for(int i=0; i<10; i++){
          pValsP(i) +=  ( pP0(i, k) >= pP(i) );
          pValsN(i) +=  ( pN0(i, k) >= pN(i) );
        }
      }
        
 
      Linalg<TPrecision>::Scale(pValsP, 1.0/nPerms, pValsP);
      Linalg<TPrecision>::Scale(pValsN, 1.0/nPerms, pValsN);

      pP.deallocate();
      pN.deallocate();
        
    };


/*
    void orthonormalize(DenseVector<TPrecision> v, int i){
       TPrecision d = Linalg<TPrecision>::DotColumn(V, i, v);
       Linalg<TPrecision>::AddScale(v, -d, V, i);
       Linalg<TPrecision>::Normalize(v);
    };


    void orthonormalize(int i1, int i2){
       TPrecision d = Linalg<TPrecision>::DotColumnColumn(V, i1, V, i2);
       Linalg<TPrecision>::ColumnAddScale(V, i2, -d, V, i1);
       Linalg<TPrecision>::NormalizeColumn(V, i2);
    };
*/



    std::string rowname(int i){
      return rownames[i];
    };

    std::string colname(int i){
      return colnames[i];
    };

    void color(int i, float c[3]){
      if(colormapping){
          return cmap.getColor(coloring(i), c);
      }else{
        c[0] = labelR(i);
        c[1] = labelG(i);
        c[2] = labelB(i);
      } 
    }

    float labelR(int index){
      return colLabelR(labels(index));
    };
    float labelG(int index){
      return colLabelG(labels(index));
    };

    float labelB(int index){
      return colLabelB(labels(index));
    };

    float colLabelR(int l){
    static float r[9] = {0.8941176, 0.2156863, 0.3019608, 0.5960784, 1.0000000, 1.0000000, 0.6509804, 0.9686275, 0.6000000};
      return r[l];
    };

    float colLabelG(int l){
      static float g[9] = {0.1019608, 0.4941176, 0.6862745, 0.3058824, 0.4980392, 1.0000000, 0.3372549, 0.5058824, 0.6000000};
      return g[l];
    };

    float colLabelB(int l){
      static float b[9] = {0.1098039, 0.7215686, 0.2901961, 0.6392157, 0.0000000, 0.2000000, 0.1568627, 0.7490196, 0.6000000};
      return b[l];

    };


 
    void getCorIndicies(int n, int *indices, bool cor){
      using namespace FortranLinalg;
          DenseVector<int> knn(n);
          DenseVector<TPrecision> knnd(n);
          SquaredEuclideanMetric<TPrecision> metric;
          DenseVector<TPrecision> p(1);
          p(0) = cor?1:-1;
          DenseMatrix<TPrecision> row(1, P.N());
          Linalg<TPrecision>::SetRow(row, 0, P, 0);
          Geometry<TPrecision>::computeKNN(row, p, knn, knnd, metric);
          for(int i=0; i<n; i++){
            indices[i] = knn(i);
          }
          knn.deallocate();
          knnd.deallocate();
          p.deallocate();
    }; 

/*
    void writeSelected(){
       int index = getSelected();
       if(index >= 0){
          DenseVector<int> knn(centered.N());
          DenseVector<TPrecision> knnd(knn.N());
          SquaredEuclideanMetric<TPrecision> metric;
          Geometry<TPrecision>::computeKNN(centered, index, knn, knnd, metric);
          std::vector<std::string> cor(20);
          std::vector<std::string> acor(20);
          for(int i=0; i<20; i++){
            int i1 = knn(i);
            std::stringstream s1;
            s1 << rownames[i1] << " " << 
            Linalg<TPrecision>::DotColumnColumn(centered, index, centered, i1);
            cor[i] = s1.str();

            int i2 = knn(knn.N()-1-i);
            std::stringstream s2;
            s2 << rownames[i1] << " " << 
            Linalg<TPrecision>::DotColumnColumn(centered, index, centered, i2);
            acor[i] = s2.str();
          }
	  IO<TPrecision>::writeStringList("cor.txt", cor);
          IO<TPrecision>::writeStringList("acor.txt", acor);
          knn.deallocate();
          knnd.deallocate();
       }
    };

   

 
    void readSearched(){
      std::list<std::string> names = IO<TPrecision>::readStringList(SEARCHFILE);
      for(std::list<std::string>::iterator it=names.begin(); it!=names.end(); ++it){
         std::string name = *it;
         for(int i=0; i<rownames.size(); i++){
             if(name.compare(rownames[i]) == 0 ){
               showRowname.push_back(i);
               break;
             }
          }  
       }
    };
*/
    
    void addShowed(std::string name){
      for(int i=0; i<rownames.size(); i++){
        if(name.compare(rownames[i]) == 0 ){
          showRowname.push_back(i);
          showProfile.push_back(i);
          break;
        }
      }  
    };

private:


  void computeNull(){
      using namespace FortranLinalg;
      if(nPerms == 0) return;

      int nSS = raw.N();
      //if(nSS > 10000){ nSS = 10000;}

      //compute null distribution
      DenseMatrix<TPrecision> tmp(raw.M(), nSS);


      Linalg<TPrecision>::Zero(pP0);
      Linalg<TPrecision>::Zero(pN0);
      for(int np = 0; np < nPerms; np++){
        for(int k=0; k < nSS; k++){
          int index = k;

          if(nSS < raw.N()){
            index = ( rand() / (double) RAND_MAX ) * raw.N();
            if(index == raw.N() ){ index = raw.N()-1; }
          }
     
          //create random permutation
          tmp(0, k) = raw(0, index);
          for(int j=1; j < tmp.M(); j++){
            int el = (  rand() / (double) RAND_MAX ) * (j+1);
            if(el == j+1){ el =j; }

            tmp(j, k) = tmp(el, k);
            tmp(el, k) = raw(j, index);
          }
        }
  
        if(useSpearman){
          //rank each row
          std::vector< Pair<TPrecision> > pairs(tmp.M());

          for(int i=0; i<tmp.N(); i++){

            for(int j=0; j<tmp.M(); j++){
              pairs[j].index = j;
              pairs[j].value = tmp(j, i);
            }
 
            std::sort(pairs.begin(), pairs.end());
            for(int j=0; j < tmp.M(); j++){
              tmp(pairs[j].index, i) = j;
            }

          }
        }

        DenseVector<TPrecision> means = Linalg<TPrecision>::SumRows(tmp);
        Linalg<TPrecision>::Scale(means, 1.0f/raw.M(), means);
        Linalg<TPrecision>::SubtractRowwise(tmp, means, tmp);
        means.deallocate();
  
 
        for(int i=0; i<tmp.N(); i++){
          TPrecision l = Linalg<TPrecision>::LengthColumn(tmp, i);
          Linalg<TPrecision>::ScaleColumn(tmp, i, 1.0f/l);
        } 
        
        PCA<TPrecision> tmpPCA(tmp, 1, false);
        DenseMatrix<TPrecision> proj = tmpPCA.project(tmp);
        tmpPCA.cleanup();

        for(int i=0; i<tmp.N(); i++){
          double dot = proj(0, i);
          for(int k=0; k<10; k++){
            if(dot > k/10.0){
              pP0(k, np)++;
            }
            else if(dot < -k/10.0){
              pN0(k, np)++;
            } 
          }
        }
        proj.deallocate();
      }
      Linalg<TPrecision>::Scale(pP0, 1.0/nSS, pP0);
      Linalg<TPrecision>::Scale(pN0, 1.0/nSS, pN0);
      tmp.deallocate();

  };
};



#endif
