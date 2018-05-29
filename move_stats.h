#ifndef MOVESTATS_H
#define MOVESTATS_H
#include <vector>
#include <map>
#include <string>
#include <cassert>


#include "Gtypes.h"
#include "Rtypes.h"
#include "TVirtualFitter.h"

class TObject;
class TH1;
class TH1F;
class TH1D;
class TProfile;
class TH2;
class TH3;
class TGraph;
class TGraphErrors;
class TMultiGraph;
class TTree;
class TFile;
class TDirectory;
class TLegend;
class TAxis;
class TPad;
class TCanvas;
class TAttLine;
class TAttFill;
class TAttMarker;
class TChain;
class THStack;
class TGraphAsymmErrors;

#include "TParameter.h"

#ifdef __MAKECINT__

#pragma link C++ class vector<TH1*>;
#pragma link C++ class vector<TH2*>;
#pragma link C++ class vector<TTree*>;
#pragma link C++ class vector<TAttLine*>;
#pragma link C++ class vector<TAttFill*>;
#pragma link C++ class vector<TAttMarker*>;
#pragma link C++ class vector<float>;
#pragma link C++ class vector<int>;
#pragma link C++ class vector<TPad*>;
#pragma link C++ class vector<TGraph*>;
#pragma link C++ class vector<TCanvas*>;
#pragma link C++ class vector< vector<TH1*> >;
#pragma link C++ class vector< TParameter<double> >;
#pragma link C++ class map<string,string>;
#pragma link C++ class vector<TGraphErrors*>;

#endif


///////////////////////////////////////////////////////////////////////////////
// functions in this file
///////////////////////////////////////////////////////////////////////////////

// set line style for all hists and graphs in memory in gDirectory
void set_ls(Int_t v, Int_t lw=-1);
// set line style for all hists and graphs in gDirectory (opens files on disk)
void set_ls_file(Int_t v, Int_t lw=-1);
// set line color for all hists and graphs in memory in gDirectory
void set_lc(Int_t v);
// set line color for all hists and graphs in gDirectory (opens files on disk)
void set_lc_file(Int_t v);
// set fill style for all hists and graphs in memory in gDirectory
void set_fs(Int_t v);
// set fill style for all hists and graphs in gDirectory (opens files on disk)
void set_fs_file(Int_t v);
// set fill color for all hists and graphs in memory in gDirectory
void set_fc(Int_t v);
// set fill color for all hists and graphs in gDirectory (opens files on disk)
void set_fc_file(Int_t v);

// set marker style for all hists and graphs in memory in gDirectory
void set_ms(Int_t v);
// set marker style for all hists and graphs in gDirectory(opens files on disk)
void set_ms_file(Int_t v);

// set marker color for all hists and graphs in memory in gDirectory
void set_mc(Int_t v);
// set marker color for all hists and graphs in gDirectory(opens files on disk)
void set_mc_file(Int_t v);

// define a greyscale palette
void grey_palette(int opt=0, Color_t=kGray); //[low->high]:0,default=[white->black] 1=inverse


// set colors and styles for a vector of stuff
//template<class T> void set_lc(std::vector<T*>& v, Color_t x);

void set_lc(std::vector<TH1*>& v, Color_t x);
void set_fc(std::vector<TH1*>& v, Color_t x);
void set_mc(std::vector<TH1*>& v, Color_t x);
void set_ls(std::vector<TH1*>& v, Style_t x);
void set_fs(std::vector<TH1*>& v, Style_t x);
void set_ms(std::vector<TH1*>& v, Style_t x);

void set_range(std::vector<TH1*>& v, 
	       std::vector<float>& r1, std::vector<float>& r2);
void rebin(std::vector<TH1*>& v, std::vector<int>& r);

// make vectors from arrays
std::vector<int> make_vector(const int* p1, const int* p2);
std::vector<float> make_vector(const float* p1, const float* p2);


// rearrange and color multiple stats boxes in gPad
// boxes are colored according to the histogram or graph line color
// opt = "v" arrange boxes vertically
//       "h" arrange boxes horizontally
// x,y = NDC coordinates of block of stats boxes
//       allows you to move them around left and right, up and down
//       if x,y=-1  use gStyle->GetStatsX,Y() for x and y
void move_stats(const char* opt="v", float x=-1, float y=-1);

// move stats on all pads of all canvases
void move_stats_all_cans(const char* opt="v", float x=-1, float y=-1);

// fix figure maximum on gPad so all histograms fit vertically
void fix_maximum();

// fix maximum on all pads of all canvases
void fix_maximum_all_cans();

// add a histogram to a stack and set colors at the same time
void add_to_stack(THStack*, TH1*, Color_t fill);

// call print on all open canvases
void write_all_cans(const char* prefix="", const char* type="eps");

// center axis titles
void center_titles(const char* opt="XY");
void center_titles_all_cans(const char* opt="XY");

// draw a legend, fixing the box color and border size
// one should be able to set these via gStyle
void draw_leg(TLegend* leg);

// find histograms, stick them in a vector
int get_histos(std::vector<TH1*>& v, const char* names);

// get things from a file, close the file (get_t is an exception)
// histograms are normally owned by the file but this function
// removes that ownership (so the file may be closed)
// get_t gets a tree but can't close the file
TH1* get_h1(const char* hname, const char* fname);
TH1* get_h1_from_can(const char* hname, const char* fname, const char* can);
TH2* get_h2(const char* hname, const char* fname);
TH3* get_h3(const char* hname, const char* fname);
TGraph* get_g(const char* hname, const char* fname);
TMultiGraph* get_mg(const char* hname, const char* fname);
TTree* get_t(const char* hname, const char* fname);
TChain* get_c(const char* hname, const char* fname);

TH1* get_h2_proj(const char* hname, const char* fname,const char* opt,
		 int b1=-1, int b2=-1, 
		 const char* name="", float scale_factor=1);

TH1* make_h2_proj(TH2* h2, const char* opt, int b1=-1, int b2=-1, 
		  const char* name="", float scale_factor=1);

void put(TObject* obj, const char* name, const char* file);

std::vector<TH1*> get_projections(TH2* h, const char* opt, 
				   std::vector<int>& b1, std::vector<int>& b2,
				   const char* name_form,
				   float scale_factor);

std::vector<TH1*> get_projections(const char* hname, const char* fname, 
				  const char* opt, 
				  std::vector<int>& b1, std::vector<int>& b2,
				  const char* name_form,
				  float scale_factor);

std::vector<TH1*> get_projections(TH2* h2, 
				  const char* opt, 
				  std::vector<float>& b1, 
				  std::vector<float>& b2,
				  const char* name_form,
				  float scale_factor);

std::vector<TH1*> get_projections(const char* hname, const char* fname, 
				  const char* opt, 
				  std::vector<float>& b1, 
				  std::vector<float>& b2,
				  const char* name_form,
				  float scale_factor);


std::vector<TH1*> get_projections(TH2* h2,
				  const char* opt, 
				  float low, float high, int n,
				  const char* name_form,
				  float scale_factor);



std::vector<TH1*> get_projections(const char* hname, const char* fname, 
				  const char* opt, 
				  float low, float high, int n,
				  const char* name_form,
				  float scale_factor);

void agg_axis(TAxis* a, float llim, float hlim, unsigned int n, 
	 std::vector<int>& low, std::vector<int>& high);

// choose_h3_axes used by make_h3_proj
const char* choose_h3_axes(TH3* h3, const char* opt, 
			   TAxis* & a_axis, TAxis* & b_axis, TAxis* & p_axis);

TH1* make_h3_proj(TH3* h3, const char* opt, int a1=-1, int a2=-1,
		  int b1 = -1, int b2 = -1,
		  const char* name="", float scale_factor=1);

std::vector<TH1*> get_h3_projections(const char* hname, const char* fname,
				     const char* opt,
				     std::vector<int>& a1,
				     std::vector<int>& a2,
				     std::vector<int>& b1,
				     std::vector<int>& b2,
				     const char* name_form,
				     float scale_factor);

std::vector<TH1*> get_h3_projections(TH3* h3, const char* opt,
				     std::vector<int>& a1,
				     std::vector<int>& a2,
				     std::vector<int>& b1,
				     std::vector<int>& b2,
				     const char* name_form,
				     float scale_factor);


std::vector<TH1*> get_h3_projections(TH3* h3, const char* opt,
				    TH2* h2,
				    const char* name_form,
				    float scale_factor);


// plot vectors of histograms
void plot_many(std::vector<TH1*>& v, const std::vector<TPad*>& pads, 
	       const char* dopt="");

void plot_many(std::vector<TH1*>& v, int nh, int nv, 
	       const char* xtit, const char* ytit, const char* opt,
	       std::vector<TCanvas*>& cans, std::vector<TPad*>& pads,
	       const char* cantit="can",
	       float xmargin=0.07, float ymargin=0.07,
	       float xsubmargin=0.01, float ysubmargin=0.01);

void print_canvases(std::vector<TCanvas*>& v, const char* name="can_%i", const char* extensions="eps,gif,root");
void print_canvas(TCanvas* can,const char* name,const char* extensions="eps,gif,root");
void print_canvas(const char* cname ,const char* name,const char* extensions="eps,gif,root");

//  Histogram with logarithmic x-axis scale
// sets bin sizes such that when gPad->Logx,y() is called
// the bin widths look uniform.
// This is similar to filling a histogram with log10(quantity)
// but the axis labels will look like 10^N rather than N
TH1* log10_h1(const char* hname, const char* tit, 
	    int nbins, float low, float high);
TProfile* log10_prof(const char* hname, const char* tit, 
	    int nbins, float low, float high);

// natural logs, unfortunately the bin sizes will
// not look uniform when gPad->Logx,y() is called since
// that's a log10 scale. You can't do a Ln scale in root. 
// They ought to add it.
TH1* log_h1(const char* hname, const char* tit, 
	    int nbins, float low, float high);

TProfile* log_prof(const char* hname, const char* tit, 
	    int nbins, float low, float high);

// couldn't make this templated version work with root 4.02.00
// due to a bug in CINT (found by Brett actually)
//template<class T> T* log_h1(const char* hname, const char* tit, 
//			    int nbins, float low, float high);


// histogram operations
// vanilla root doesn't include methods two divide two histograms
// without getting a pointer to one of them... this is annoying
// if "new" is included in the option then a new histogram will be returned
// otherwise the histogram referenced by n1 will modified and returned
//
// divide two histograms, opt is passed to the TH1::Divide() call
//TH1* dhist(const char* n1, const char* n2, const char* opt=0);
//TH1* dhist(const char* n1, const char* n2, const char* opt, 
//	   Double_t c1=1, Double_t c2=1,);
// multiply two histograms, opt is passed to the TH1::Multiply() call
//TH1* mhist(const char* n1, const char* n2, const char* opt);
// add two histograms
///////////////////////////////////////////////////////////////////////////////

///////////////// integrate a histogram to the right /////////////////////////
TH1* sum_right(TH1* h);

//////////////  set axes ranges //////////////////////////////////
void set_axis_range(TH1* h, float nsig, float quant);
void same_axis_range(TH1* h1, TH1* h2);
void set_axis_range2(TH1* h, float nsigup, float nsigdn, float quant);

///////////// fitting utilities //////////////////////////////////

/// The following functions use gMinuit and assume you've just done a fit

// perform a chi2 scan of the recently fit function in par[par_num]
// scan over the range [ (par-epar*nsig_down) , (par+epar*nsig_up) ]
// in ndiv steps, both up and down from the best fit point
// the graph shows delta chi2 vs. the parameter value
TGraph* chisquare_scan(Int_t par_num, Int_t nsig_up=3, Int_t nsig_down=3, Int_t ndiv=10);
// do the scan and also draw (convenience function)
TGraph* draw_chisquare_scan(Int_t par_num, Int_t nsig_up=3, Int_t nsig_down=3, Int_t ndiv=10);
// get ncont  two-dimensional fit contours,
// starting at sig_start  (in std. deviations), stepping by sig_step
/// p1,p2 are function parameter numbers --> map to x,y on graphs
// this get_fit_contours and draw_fit_contours assume gMinuit was used
// to do the fit
// scan_fit_2d is more general, and maybe more robust but probably slower

TMultiGraph* get_fit_contours(Int_t p1, Int_t p2, Int_t ncont=3, Float_t sig_start=1.0, Float_t sig_step=1.0, Int_t npts=30);

// get the contours and draw them (convenience function)
TMultiGraph* draw_fit_contours(Int_t p1, Int_t p2, Int_t ncont=3, Float_t sig_start=1.0, Float_t sig_step=1.0, Int_t npts=30);


//TMultiGraph* get_fit_contours_scan(TH2* h, Int_t p1, Int_t p2, Int_t ncont=3, Float_t sig_start=1.0, Float_t sig_step=1.0);
// scan p1,p2 in bins of input histogram h, minimizing other pars in each bin.
//void scan_fit_2d(TVirtualFitter* vfit, TH2* h,Int_t p1,Int_t p2);
void contours_on_surface(TH2* h, const char* dchi2="2.3 4.61");
std::vector<TGraph*> get_hist_contour(int level);


// compute chi2
double chi2_h1(TH1* h1, TH1* h2, const char* opt="");
double chi2_h1(TH1* h1, TH1* h2, double& chi2, int& N, const char* opt="");

// these functions are used by the ones above
// they should work fine when called alone but the author
// hasn't used them that way
TGraph* best_fit_point(Int_t p1, Int_t p2);
TGraph* two_dim_contour(Int_t p1, Int_t p2, Double_t sig=1.0, Int_t npts=30);


/// chebyshev polynomials

double chebyshev(int order, double x);

double cheby10(double* xx, double* pp);

/*
///// templated stuff
// couldn't make this work with 4.02.00... see above
template<class T> T* log_h1(const char* hname, const char* tit, int nbins, float low, float high){
  // xbins construction taken from rootalk : P.Sizun
  
  Double_t logxmin = TMath::Log(xmin);
  Double_t logxmax = TMath::Log(xmax);
  Double_t binwidth = (logxmax-logxmin)/nbins;
  Double_t xbins[nbins+1];
  xbins[0] = xmin;
  for (Int_t i=1;i<=nbins;i++) {
    //    xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);
    xbins[i] = xmin + exp(logxmin+i*binwidth);
  }
    
  T* h = new T(hname, tit, nbins,xbins);  
  return h;
}
#pragma link C++ function log_h1<TProfile>(const char* , const char*, int, float, float);
#pragma link C++ function log_h1<TH1F>(const char* , const char*, int, float, float);
#pragma link C++ function log_h1<TH1D>(const char* , const char*, int, float, float);
*/

float get_pot(const char* file, int beam, const char* name="h_tortgt");
float get_pot_mc(const char* file, float scale=2.42e13, const char* name="h_snarls");

// filter TTree t according to selection and write out passing events to fout
// only write variables indicated in vars, which is a comma or space seperated list
TTree* filter_tree(TTree* t, const char* selection, const char* vars, TFile* fout);
// interface taking ntuple name and "files" for TChain
// as well as output file name
void filter_tree(const char* ntuple, const char* files, const char* selection, const char* vars, const char* fout);

/*
template<class T> void set_lc(std::vector<T*>& v, Color_t x){
  for(unsigned int i=0; i<v.size(); i++) v[i]->SetLineColor(x);
}

#pragma link C++ function set_lc<TH1>(std::vector<TH1*>&,Color_t);
void set_lc(std::vector<TH1*>&, Color_t);

//#pragma link C++ function set_lc<TH1>;
*/

// rebin input likethis
// handy with make_h1() to convert a uniformly binned histogram
// into one with asymmetric bins
TH1* rebin_like(TH1* input, TH1* likethis, const char* newname);
TH2* rebin_like(TH2* input, TH2* likethis, const char* newname);

// make a histogram, possibly with asymmetric bins, in a convenient way
// make a TH1F according to binning
// here is how it works:
// opt="W": low_limit, nbins binw, nbins binw, etc.
// opt="R": low_limit, nbins high_limit, nbins high_limit, etc.
TH1* make_h1(const char* name, const char* title, const char* binning, const char* opt);

TH2* make_h2(const char* name, const char* title, const char* binning_x, const char* opt_x, const char* binning_y, const char* opt_y);


// utlitity function used by make_h1
Int_t make_binning(const char* binning, const char* opt, std::vector<Double_t>& low_edges, std::vector<Double_t>& high_edges);

// divide each bin by its width, and multiply by m
void divide_binw(TH1* h, float m=1.0);

// divide each bin by its width, and multiply by m
void divide_binw_2d(TH2* h, float m=1.0);

// shift/translate/move a TH1
TH1* translate_h1(TH1* hin, double x, const char* hname="%s_tr");


// utilities for merging histograms
/*
// removed due to a compilation problem in root 5.34
// see the .C file for details
void myhadd(const char* outfile, const char* infiles);

void myhadd(TDirectory* target, const char* infiles);
*/

void copy_hists(TDirectory* target, TDirectory* from);
void merge_hists(TDirectory* target, TDirectory* from);

void print_exposure(TDirectory* dir);

//////////////// utilities for I/O //////////////////////////////

//removed due to a compilation problem in root 5.34
// see the .C file for details
//void glob_files(const char* infiles, std::vector<std::string>& v);

void write_parameter(const char* name, double val, TDirectory* dir=0);
//void read_parameter(const char* name, double& val, TDirectory* dir=0);

// read from an ascii file and construct TParameters as
// <begin key>
// par1 val1
// par2 val2
// ...
// <end key>
void pars_from_ascii(const char* filename, std::vector<TParameter<double> >& v,
		const char* begin_key="#pars_begin#",
		const char* end_key="#pars_end#");
void pars_from_dir(TDirectory* d, std::vector<TParameter<double> >& v);

// use input TParameter vector v to define a TTree with name,title in dir
// v should not have duplicate names
TTree* define_tree_from_pars(TDirectory* dir, const char* name, const char* title, std::vector<TParameter<double> >& v);

void fill_tree_from_pars(TTree* t, std::vector<TParameter<double> >& v);



/////////////// plot collection ////////////////////////////////
//
// a collection of plots with the same attributes
//
////////////////////////////////////////////////////////////////

class plot_collection {
 public:
  plot_collection(TH3* h, float scale_factor, const char* projopt, const char* projtit, int* alow, int* ahigh, int na, int* blow, int* bhigh, int nb);

  plot_collection(const char* hname, const char* fname, float scale_factor, const char* projopt, const char* projtit, int* alow, int* ahigh, int na, int* blow, int* bhigh, int nb);


  plot_collection(TH2* h, float scale_factor, const char* projopt, const char* projtit, int* alow, int* ahigh, int na);

  plot_collection(const char* hname, const char* fname, float scale_factor, const char* projopt, const char* projtit, int* alow, int* ahigh, int na);


  //  plot_collection(TH2* h, float scale_factor, const char* projopt, int* alow, int* ahigh, int na);
  
  std::vector<TH1*>& get_histos() {return v;}
  void set_leg(const char* label, const char* opt){leg_label=label; leg_opt=opt;}
  void set_attline(Int_t c, Int_t s=-2, Int_t w=-2);//as in TAttLine()
  void set_attfill(Int_t c, Int_t s=-2);// as in TAttFill()  
  void set_attmarker(Int_t c, Int_t s=-2, Int_t z=-2); // as in TAttMarker()
  void set_range(float* d, float* u, unsigned int n);// d[i]>u[i] does nothing
  void set_rebin(int* r, unsigned int n);// 0 does nothing


  std::vector<TH1*> v;
  std::string leg_label;
  std::string leg_opt;
  std::vector<float> range_down;
  std::vector<float> range_up;

};

#ifdef __MAKECINT__
#pragma link C++ class plot_collection;
#pragma link C++ class vector<plot_collection*>;
#endif

////////////////////////// pc_plotter ///////////////////////////
//
// a class to draw plot_collections
//
/////////////////////////////////////////////////////////////////

class pc_plotter {
 public:
  pc_plotter(int nh, int nv, const char* xtit, const char* ytit, const char* cantit = "can", float xmargin=0.07, float ymargin=0.07);
  void add_pc(plot_collection* pc, const char* dopt="");
  void draw_plots();
  void set_legend(Double_t x1=0.6, Double_t y1=0.6, Double_t x2=0.85, Double_t y2=0.85, const char* header="");
  void draw_legend(unsigned int ipad);
  

  std::vector<plot_collection*> pcs;
  std::vector<TCanvas*> cans;
  std::vector<TPad*> pads;
  std::vector<std::string> opts;
  int nh;
  int nv;
  std::string xt;
  std::string yt;
  std::string ct;
  float xm;
  float ym;
  TLegend* leg; 
};


void make_minos_style();

TGraphAsymmErrors* setPoissonErrors (TGraphAsymmErrors* g1, TH1D* h1);


#endif
