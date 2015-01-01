#ifndef DATA_H
#define DATA_H

#include "DenseMatrix.h"
#include "Linalg.h"
#include "LinalgIO.h"
#include "PCA.h"

#include <vector>
#include <map>
#include <set>

template<typename TPrecision>
class Data{



  public:
 
    DenseMatrix<TPrecision> units; 
    DenseMatrix<TPrecision> raw;
    DenseMatrix<TPrecision> centered;
    DenseVector<int> labels;
    std::vector<std::string> names;
   
     
    PCA<TPrecision> *pca; 
    std::map<int, PCA<TPrecision>* > lPCA;
 




    Data( DenseMatrix<TPrecision> &r, DenseVector<int> l,
          std::vector<std::string> n ) : raw(r),labels(l), names(n){
     
      units = DenseMatrix<TPrecision>(raw.M(), raw.M());
      double ul = sqrt(1.0/raw.M());
      Linalg<TPrecision>::Set(units, -ul);
      for(int i=0; i<units.N(); i++){
        units(i,i) = 1.0 - ul;
      };      

 
      centered = Linalg<TPrecision>::Copy(raw);
      DenseVector<TPrecision> means = Linalg<TPrecision>::SumRows(centered);
      Linalg<TPrecision>::Scale(means, 1.0f/raw.M(), means);
     
      Linalg<TPrecision>::SubtractRowwise(centered, means, centered);
   
      for(int i=0; i<centered.N(); i++){
        TPrecision l = Linalg<TPrecision>::LengthColumn(centered, i);
        Linalg<TPrecision>::ScaleColumn(centered, i, 1.0f/l);
      } 
      pca = new PCA<TPrecision>(centered, centered.M()-1, false);


      //per label pca
      std::set<int> nlabels;
      for(int i=0; i<labels.N(); i++){
        nlabels.insert(labels(i));
      };
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
        tmp.deallocate();
      }
    };




    void shave(std::vector<int> keep){
      DenseMatrix<TPrecision> c2(raw.M(), keep.size());
      DenseVector<int> l2(keep.size());
      std::vector<std::string> n2;

      for(int i=0; i<keep.size(); i++){
        int index = keep[i];
        Linalg<TPrecision>::SetColumn(c2, i, centered, index);
        l2(i) = labels(index);
        n2.push_back(names[index]);
      }
      names = n2;
      centered.deallocate();
      centered = c2;
      labels.deallocate();
      labels = l2;
      pca->cleanup();
      pca = new PCA<TPrecision>(centered, centered.M()-1, false);
    };




    std::string name(int i){
      return names[i];
    };

    double colorLabelR(int i){
      return colorR(labels(i));
    }

    double colorLabelG(int i){
      return colorG(labels(i));
    };

    double colorLabelB(int i){
      return colorB(labels(i));
    };

    double colorR(int l){
    static double r[9] = {0.8941176, 0.2156863, 0.3019608, 0.5960784, 1.0000000, 1.0000000, 0.6509804, 0.9686275, 0.6000000};
      return r[l];
    };

    double colorG(int l){
      static double g[9] = {0.1019608, 0.4941176, 0.6862745, 0.3058824, 0.4980392, 1.0000000, 0.3372549, 0.5058824, 0.6000000};
      return g[l];
    };

    double colorB(int l){
      static double b[9] = {0.1098039, 0.7215686, 0.2901961, 0.6392157, 0.0000000, 0.2000000, 0.1568627, 0.7490196, 0.6000000};
      return b[l];

    };
   
};



#endif
