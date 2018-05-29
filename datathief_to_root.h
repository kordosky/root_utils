#ifndef DT2ROOT
#define DT2ROOT
#include <vector>
//class TGraph;
class TGraphErrors;

struct dt_data{
  public:
  double x,y,euy,ely;
  
};

class dt_to_root {
 public:
  bool read(const char* file);
  void clear();
  TGraphErrors* make_graph(const char* name, const char* title);
 private:
  int nmax;
  std::vector<dt_data> v;

};


#endif
