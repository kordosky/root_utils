#include "datathief_to_root.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include "TGraphErrors.h"

bool dt_to_root::read(const char* file){

  std::ifstream f(file);
  if(!f) return false;
  std::string b;
  std::getline(f,b);
  std::cout<<b<<std::endl;
  int nmax=0;
  while(std::getline(f,b)){    
    std::istringstream s(b);
    dt_data p;
    double tmp;
    int nvar=0;
    s>>tmp; if(s.good()) { nvar++; p.x=tmp;}
    s>>tmp; if(s.good()) { nvar++; p.y=tmp;}
    s>>tmp; if(s.good()) { nvar++; p.ely=tmp;}
    s>>tmp; if(s.good()) { nvar++; p.euy=tmp;}
    if(nvar>nmax) nmax=nvar;
    if(nvar==3) p.euy = p.ely;
    if(nvar==2) p.euy = p.ely = 0.0;
    //    std::cout<<p.x<<" "<<p.y<<" "<<p.ely<<" "<<p.euy<<std::endl;
    v.push_back(p);
  }
  return true;
}

void dt_to_root::clear(){
  v.clear();
}

TGraphErrors* dt_to_root::make_graph(const char* name, const char* title){
  TGraphErrors* g = new TGraphErrors();
  g->SetNameTitle(name,title);
  for(uint i=0; i<v.size(); i++){
    g->SetPoint(i,v[i].x,v[i].y);
    double big = v[i].ely;
    if(v[i].euy>big) big=v[i].euy;
    g->SetPointError(i,0.0, big);    
  }
  return g;
}
