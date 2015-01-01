#ifndef FONT_H
#define FONT_H


class Font{
protected:
  double size;
public:
  Font(){
    size = 12;
  };

  virtual void renderString(std::string text, int x, int y, int z=0, bool vertical=false) = 0;

  void setSize(double s){
    size = s;
  };
};

#endif
