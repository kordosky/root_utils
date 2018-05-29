#include "move_stats.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
//#include <glob.h>

#include "TObject.h"
#include "TROOT.h"
#include "TString.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TIterator.h"
#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TPaveStats.h"
#include "TDirectory.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttMarker.h"
#include "TColor.h"
#include "TKey.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TProfile.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TIterator.h"
#include "TLatex.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TEventList.h"
#include "TChain.h"
#include "TF1.h"
#include "TTreeFormula.h"
#include "TRandom.h"
//#include "TFitterMinuit.h"
#include "TFitter.h"
//#include "TFitterFumili.h"
//#include "TFumili.h"
#include "TLinearFitter.h"
#include "TFractionFitter.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"

using namespace std;

void write_all_cans(const char* prefix, const char* type){
   if(!gPad) return;
//   TVirtualPad* save = gPad;
   
   // loop on all canvases
   TIter canit(gROOT->GetListOfCanvases());
   TObject* obj;
   while((obj=canit.Next())){
      if(obj->InheritsFrom("TCanvas") ){
	 TCanvas* can = static_cast<TCanvas*>(obj);
	 can->cd();
	 TString name=prefix; name+=can->GetName();
	 //	 name+="."; name+=type;
	 //	 can->Print(name.Data(), type);
	 print_canvas(can,name.Data(),type);
      }
   }
}

void fix_maximum_all_cans(){
   if(!gPad) return;
   TVirtualPad* save = gPad;
   
   // loop on all canvases
   TIter canit(gROOT->GetListOfCanvases());
   TObject* obj;
   while((obj=canit.Next())){
      if(obj->InheritsFrom("TCanvas") ){
	 TCanvas* can = static_cast<TCanvas*>(obj);
	 can->cd();
	 TIter it(can->GetListOfPrimitives());
	 while((obj=it.Next())){
	    if(obj->InheritsFrom("TPad") ){
	       TPad* pad = static_cast<TPad*>( obj);
	       pad->cd();
	       fix_maximum();
	    }
	 }
      }
   }
   save->cd();
}

void move_stats_all_cans(const char* opt, float x, float y){
   if(!gPad) return;
   TVirtualPad* save = gPad;
   
   // loop on all canvases
   TIter canit(gROOT->GetListOfCanvases());
   TObject* obj;
   while((obj=canit.Next())){
      if(obj->InheritsFrom("TCanvas") ){
	 TCanvas* can = static_cast<TCanvas*>(obj);
	 can->cd();
	 TIter it(can->GetListOfPrimitives());
	 while((obj=it.Next())){
	    if(obj->InheritsFrom("TPad") ){
	       TPad* pad = static_cast<TPad*>( obj);
	       pad->cd();
	       move_stats(opt,x,y);
	    }
	 }
      }
   }
   save->cd();
}
void fix_maximum(){

   // sets the maximum in the current pad so that all histograms are on scale
     gPad->Modified();
     gPad->Update();

     int nps=0;
     float max=0.0;
     // iterate over primitives looking for histograms
     TH1* hvec[10000];
     TIter it(gPad->GetListOfPrimitives());
     TObject* obj;
     while((obj=it.Next())){
	  if(obj->InheritsFrom("TH1") ){
//	     && !obj->InheritsFrom("TH2") 
//	     && !obj->InheritsFrom("TH3")){
	       // is a 1d histogram
	       // try to get the stats box
	       TH1* h= static_cast<TH1*>(obj);
	       hvec[nps]=h;

	       if(h->GetMaximum()>max) max=h->GetMaximum();

	       nps++;
	  }
	  else if(obj->InheritsFrom("TPad")){
	     // make it recursive
	     TVirtualPad* spad = gPad;
	     TPad* npad = static_cast<TPad*>(obj);
	     npad->cd();
	     fix_maximum();
	     spad->cd();
	  }
     }
     // include an additional bit of whitespace at the top of the figure
     if(gPad->GetLogy()) max*=2.3;
     else max*=1.1;
     for(int i=0; i<nps; i++){
	TH1* h = hvec[i];
	h->SetMaximum(max);
     }
     gPad->Modified();
     gPad->Update();
}

// add a histogram to a stack and set colors at the same time
void add_to_stack(THStack* hstack, TH1* h, Color_t fill){
    hstack->Add(h); 
    h->SetFillColor(fill);
    h->SetLineColor(fill);
}


//void move_stats(){
void move_stats(const char* opt, float x, float y){

    // Moves and aligns the stats boxes for histograms in the current pad
    // 
    // author: mike kordosky
    // kordosky@fnal.gov
    //
 
     gPad->Modified();
     gPad->Update();
     
     std::string opts = opt;
     
     TPaveStats* ps_array[50]={0}; // Limitation: No more than 50 histograms in pad
     int nps=0;
     // iterate over primitives looking for histograms
     TIter it(gPad->GetListOfPrimitives());
     TObject* obj;
     int first=1;
     int blow=0;
     int bhigh=1;
     while((obj=it.Next())){
	  if(obj->InheritsFrom("TH1") ){
	       // is a 1d histogram
	       // try to get the stats box
	       TH1* h= static_cast<TH1*>(obj);


	       // this sordid busines resets the axis limits based on
	       // the first histogram it finds in the pad
	       // the reason is that this will make the stats box display
	       // consistent
	       if(first){
		  blow=h->GetXaxis()->GetFirst();
		  bhigh=h->GetXaxis()->GetLast();
		  first=0;
	       }
	       else{
		  h->GetXaxis()->SetRange(blow,bhigh);
	       }

	       TObject* ps_obj = h->FindObject("stats");
	       if(ps_obj&&ps_obj->InheritsFrom("TPaveStats")){
		    TPaveStats* ps = static_cast<TPaveStats*>(ps_obj);
		    // rename stats box
		    TString s = h->GetName();
		    s+="_stats";

		    ps->SetName(s.Data());
		    ps->SetTextColor(h->GetLineColor());
		    // store in array
		    ps_array[nps]=ps;
		    nps++;
	       }	       
	  }
	  else if(obj->InheritsFrom("TGraph") ){
	       // is a 1d histogram
	       // try to get the stats box
	       TGraph* gr= static_cast<TGraph*>(obj);


	       // this sordid busines resets the axis limits based on
	       // the first histogram it finds in the pad
	       // the reason is that this will make the stats box display
	       // consistent
	       if(first){
		  blow=gr->GetXaxis()->GetFirst();
		  bhigh=gr->GetXaxis()->GetLast();
		  first=0;
	       }
	       else{
		  gr->GetXaxis()->SetRange(blow,bhigh);
	       }

	       TObject* ps_obj = gr->FindObject("stats");
	       if(ps_obj&&ps_obj->InheritsFrom("TPaveStats")){
		    TPaveStats* ps = static_cast<TPaveStats*>(ps_obj);
		    // rename stats box
		    TString s = gr->GetName();
		    s+="_stats";

		    ps->SetName(s.Data());
		    ps->SetTextColor(gr->GetLineColor());
		    // store in array
		    ps_array[nps]=ps;
		    nps++;
	       }	       
	  }
	  else if(obj->InheritsFrom("TPad") ){
	     // make it recursive
	     TVirtualPad* spad = gPad;
	     TPad* npad = static_cast<TPad*>(obj);
	     npad->cd();
	     move_stats();
	     spad->cd();
	  }
     }
     if(nps<=0) return;
     // configuration parameters
     //     double tot_h = 0.6;
     double tot_h = gPad->GetTopMargin()+gPad->GetBottomMargin() 
       + gStyle->GetLabelSize("X")+gStyle->GetLabelOffset("X")
	 // this line from TGAxis::PaintAxis(): title offset is a scale factor
       + gStyle->GetTitleOffset("X")*1.3*gStyle->GetTitleSize("X");
     tot_h=1-tot_h;
     double stats_h = tot_h/nps;
     double start_h=gStyle->GetStatX()-gStyle->GetStatW();
     double start_v=gStyle->GetStatY();
     bool vertical=true;
     if((opts.find("h")!=std::string::npos)
	|| (opts.find("H")!=std::string::npos)){
       vertical=false;
       start_h=gStyle->GetStatX();
       start_v=gStyle->GetStatY()-gStyle->GetStatH();
       double tot_h = gPad->GetLeftMargin()+gPad->GetRightMargin()
	 + gStyle->GetLabelSize("Y") + gStyle->GetLabelOffset("Y")
	 // this line from TGAxis::PaintAxis(): title offset is a scale factor
	 + gStyle->GetTitleOffset("Y")*1.6*gStyle->GetTitleSize("Y");
       tot_h=1-tot_h;
     }
     
     // a fudge factor accounting for space taken up by axis

     if(y>=0){
       start_v=y;
     }
     if(x>=0){
       start_h=x;
     }

     for(int i=0; i<nps; i++){
	  TPaveStats* p = ps_array[i];
//	  p->Print();

//	  p->SetY2NDC(start_h-i*(stats_h+0.01));
//	  p->SetY1NDC(p->GetY2NDC()-stats_h);
	  if(vertical){
	    p->SetY2NDC(start_v-i*(stats_h+0.01));
	    p->SetY1NDC(p->GetY2NDC()-stats_h);
	    p->SetX2NDC(start_h+gStyle->GetStatW());
	    p->SetX1NDC(start_h);
	  }
	  else {
	    p->SetX2NDC(start_h-i*(stats_h+0.01));
	    p->SetX1NDC(p->GetX2NDC()-stats_h);
	    p->SetY2NDC(start_v+gStyle->GetStatH());
	    p->SetY1NDC(start_v);
	  }
     }

     gPad->Modified();
     gPad->Update();
}

void set_ls(Int_t v, Int_t lw){
   TIter it(gDirectory->GetList()); // list of objects in memory
   TObject* obj=0;
   while((obj=it.Next())){
      if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	 TAttLine* al = dynamic_cast<TAttLine*>(obj);
	 al->SetLineStyle(v);
	 if(lw>=0) al->SetLineWidth(lw);
      }
   }
}

void set_ls_file(Int_t v, Int_t lw){
   TIter it(gDirectory->GetListOfKeys());
   TObject* obj=0;
   TKey* key=0;
   while((key=static_cast<TKey*>(it.Next()))){
      obj=key->ReadObj();
      if(!obj) continue;
      if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	 TAttLine* al = dynamic_cast<TAttLine*>(obj);
	 al->SetLineStyle(v);
	 if(lw>=0) al->SetLineWidth(lw);	   
      }
   }
}

void set_lc(Int_t v){
   TIter it(gDirectory->GetList()); // list of objects in memory
   TObject* obj=0;
   while((obj=it.Next())){
      if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	 TAttLine* al = dynamic_cast<TAttLine*>(obj);
	 al->SetLineColor(v);
	}
   }
}

void set_lc_file(Int_t v){
     TIter it(gDirectory->GetListOfKeys());
     TObject* obj=0;
     TKey* key=0;
     while((key=static_cast<TKey*>(it.Next()))){
	obj=key->ReadObj();
	if(!obj) continue;
	if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	   TAttLine* al = dynamic_cast<TAttLine*>(obj);
	   al->SetLineColor(v);
	}
     }
}
/*
TH1* dhist(const char* n1, const char* n2, const char* opt){
   TH1* h1=0;
   h1=dynamic_cast<TH1*>(gDirectory->Get(n1));
   if(!h1){
      gROOT->FindObject(n1);
      
   }

   if(!h1){
      std::cerr<<"Cannot find '"<<n1<<"'"<<std::endl;
      return 0;
   }

   TH1* h2=0;
   h2= dynamic_cast<TH1*>(gDirectory->Get(n2));
   if(!h2){
      gROOT->FindObject(n2);
   }
   if(!h2){
      std::cerr<<"Cannot find '"<<n2<<"'"<<std::endl;
      return 0;
   }

   TString sopt(opt);
   bool nhist=false;
   if(sopt.Contains("new")) nhist=true;

}
*/


void center_titles(const char* opt){
     // iterate over primitives looking for histograms
     TIter it(gPad->GetListOfPrimitives());
     TObject* obj;
     TString sopt(opt);
     sopt.ToLower();

     while((obj=it.Next())){

	if(obj->InheritsFrom("TH1") ){
	   TH1* h= static_cast<TH1*>(obj);
	   if(sopt.Contains("x")){
	      h->GetXaxis()->CenterTitle();
	   }
	   if(sopt.Contains("y")){
	      h->GetYaxis()->CenterTitle();
	   }
	}
	else if(obj->InheritsFrom("TGraph")){
	   TGraph* g = static_cast<TGraph*>(obj);
	   TH1* h= g->GetHistogram();
	   if(sopt.Contains("x")){
	      h->GetXaxis()->CenterTitle();
	   }
	   if(sopt.Contains("y")){
	      h->GetYaxis()->CenterTitle();
	   }
	}
	else if(obj->InheritsFrom("TMultiGraph")){
	   TMultiGraph* g = static_cast<TMultiGraph*>(obj);
	   TH1* h= g->GetHistogram();
	   if(sopt.Contains("x")){
	     TAxis* ax = h->GetXaxis();
	     if(ax) ax->CenterTitle();
	   }
	   if(sopt.Contains("y")){
	     TAxis* ax = h->GetYaxis();
	     if(ax) ax->CenterTitle();
	   }
	}
	else if(obj->InheritsFrom("TF1")){
	   TF1* h = static_cast<TF1*>(obj);
	   if(sopt.Contains("x")){
	      h->GetXaxis()->CenterTitle();
	   }
	   if(sopt.Contains("y")){
	      h->GetYaxis()->CenterTitle();
	   }
	}
	else if(obj->InheritsFrom("TPad")){
	  // recursion!
	   TPad* pp = static_cast<TPad*>(obj);
	   TVirtualPad* sp = gPad;
	   pp->cd();
	   center_titles();
	   sp->cd();
	}


     }
     gPad->Modified();
     gPad->Update();
}

void center_titles_all_cans(const char* opt){
   if(!gPad) return;
   TVirtualPad* save = gPad;
   
   // loop on all canvases
   TIter canit(gROOT->GetListOfCanvases());
   TObject* obj;
   while((obj=canit.Next())){

      if(obj->InheritsFrom("TCanvas") ){
	 TCanvas* can = static_cast<TCanvas*>(obj);
	 can->cd();
	 TIter it(can->GetListOfPrimitives());
	 while((obj=it.Next())){
	    if(obj->InheritsFrom("TPad") ){
	       TPad* pad = static_cast<TPad*>( obj);
	       pad->cd();
	       center_titles(opt);
	    }
	 }
	 center_titles(opt);
      }

   }
   save->cd();
}

TH1* sum_right(TH1* h){
   if(!h) return 0;
   TString snew=h->GetName();
   snew+="_sumr";
   TH1* hnew;
   if(h->InheritsFrom("TProfile")){
      TProfile* hprof = dynamic_cast<TProfile*>(h);
      hnew = hprof->ProjectionX("_sumr");
   }
   else hnew = (TH1*) h->Clone(snew.Data());
   TAttLine* hal = dynamic_cast<TAttLine*>(h);
   hal->Copy(dynamic_cast<TAttLine&>(*hnew));
   TAttFill* haf = dynamic_cast<TAttFill*>(h);
   haf->Copy(dynamic_cast<TAttFill&>(*hnew));
   TAttMarker* ham = dynamic_cast<TAttMarker*>(h);
   ham->Copy(dynamic_cast<TAttMarker&>(*hnew));

   for(int i=2; i<=hnew->GetNbinsX(); i++){
     //      std::cout<<" bin["<<i-1<<"]="<<hnew->GetBinContent(i-1)
     //	       <<"  bin["<<i<<"]="<<hnew->GetBinContent(i)<<std::endl;
      double bc=hnew->GetBinContent(i-1) + hnew->GetBinContent(i);
      //      double be=sqrt(hnew->GetBinError(i-1)*hnew->GetBinError(i-1)
      //		     + hnew->GetBinError(i)*hnew->GetBinError(i));
      //      std::cout<<bc<<"  "<<be<<std::endl;
      hnew->SetBinContent(i,bc);
//      hnew->SetBinError(i,be);
   }
   return hnew;
}

void set_fs(Int_t v){
   TIter it(gDirectory->GetList()); // list of objects in memory
   TObject* obj=0;
   while((obj=it.Next())){
      if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	 TAttFill* al = dynamic_cast<TAttFill*>(obj);
	 al->SetFillStyle(v);
	}
   }
}

void set_fs_file(Int_t v){
     TIter it(gDirectory->GetListOfKeys());
     TObject* obj=0;
     TKey* key=0;
     while((key=static_cast<TKey*>(it.Next()))){
	obj=key->ReadObj();
	if(!obj) continue;
	if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	   TAttFill* al = dynamic_cast<TAttFill*>(obj);
	   al->SetFillStyle(v);
	}
     }
}


void set_fc(Int_t v){
   TIter it(gDirectory->GetList()); // list of objects in memory
   TObject* obj=0;
   while((obj=it.Next())){
      if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	 TAttFill* al = dynamic_cast<TAttFill*>(obj);
	 al->SetFillColor(v);
	}
   }
}

void set_fc_file(Int_t v){
     TIter it(gDirectory->GetListOfKeys());
     TObject* obj=0;
     TKey* key=0;
     while((key=static_cast<TKey*>(it.Next()))){
	obj=key->ReadObj();
	if(!obj) continue;
	if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	   TAttFill* al = dynamic_cast<TAttFill*>(obj);
	   al->SetFillColor(v);
	}
     }
}


void set_ms(Int_t v){
   TIter it(gDirectory->GetList()); // list of objects in memory
   TObject* obj=0;
   while((obj=it.Next())){
      if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	 TAttMarker* al = dynamic_cast<TAttMarker*>(obj);
	 al->SetMarkerStyle(v);
	}
   }
}

void set_ms_file(Int_t v){
     TIter it(gDirectory->GetListOfKeys());
     TObject* obj=0;
     TKey* key=0;
     while((key=static_cast<TKey*>(it.Next()))){
	obj=key->ReadObj();
	if(!obj) continue;
	if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	   TAttMarker* al = dynamic_cast<TAttMarker*>(obj);
	   al->SetMarkerStyle(v);
	}
     }
}


void set_mc(Int_t v){
   TIter it(gDirectory->GetList()); // list of objects in memory
   TObject* obj=0;
   while((obj=it.Next())){
      if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	 TAttMarker* al = dynamic_cast<TAttMarker*>(obj);
	 al->SetMarkerColor(v);
	}
   }
}

void set_mc_file(Int_t v){
     TIter it(gDirectory->GetListOfKeys());
     TObject* obj=0;
     TKey* key=0;
     while((key=static_cast<TKey*>(it.Next()))){
	obj=key->ReadObj();
	if(!obj) continue;
	if((obj->InheritsFrom("TH1")) || (obj->InheritsFrom("TGraph")) ){
	   TAttMarker* al = dynamic_cast<TAttMarker*>(obj);
	   al->SetMarkerColor(v);
	}
     }
}

void set_range(std::vector<TH1*>& v, 
	       std::vector<float>& r1, std::vector<float>& r2){
  for(unsigned int i=0; i<r1.size(); i++){
    v[i]->GetXaxis()->SetRangeUser(r1[i],r2[i]);
  }
  
}


void grey_palette(int opt, Color_t color){
  // [low->high]:  0,default=[white->black] 1=inverse
  TColor* tc = gROOT->GetColor(color);
  
  const int ncol=10;
  double red[ncol];
  double green[ncol];
  double blue[ncol];
  double stops[ncol];
  double dcol=-1/double(ncol);
  float R,G,B; tc->GetRGB(R,G,B);
  double gray=1.0;
  if(opt==1){
    dcol = 1/double(ncol);
    gray = 0.0;
  }

  for (int j = 0; j < ncol; j++) {
    //   ...... Define color with RGB equal to : gray, gray, gray .......
    stops[j]=double(j)/double(ncol-1);
    red[j]=gray*R;
    blue[j]=gray*B;
    green[j]=gray*G;
    //   cout<<gray<<endl;
    gray += dcol;
    cout<<red[j]<<" "<<blue[j]<<" "<<green[j]<<endl;
  }
  const int totcol=50;
#if ROOT_VERSION_CODE < ROOT_VERSION(5,16,0)
  gStyle->CreateGradientColorTable(ncol, stops, red,green,blue,totcol);
#else
  TColor::CreateGradientColorTable(ncol, stops, red,green,blue,totcol);
#endif
}

void rebin(std::vector<TH1*>& v, std::vector<int>& r){
  for(unsigned int i=0; i<r.size(); i++){
    if(r[i]==0 || r[i]==1) continue;
    else v[i]->Rebin(r[i]);
  }
}

std::vector<int> make_vector(const int* p1, const int* p2){
  return std::vector<int>(p1,p2);
}

std::vector<float> make_vector(const float* p1, const float* p2){
  return std::vector<float>(p1,p2);
}



int get_histos(std::vector<TH1*>& v, const char* names){
   std::istringstream iss(names);
   int cntr=0;
   while(!iss.eof()){
      std::string hn;
      iss>>hn;
      std::cout<<hn<<std::endl;
      TH1* h = dynamic_cast<TH1*>(gDirectory->Get(hn.c_str()));
      if(h){
	 v.push_back(h);
	 cntr++;
      }
      else{
	 std::cout<<"Could not locate "<<hn<<" or it's not a TH1."<<std::endl;
      }
   }

   return cntr;
}

void set_axis_range(TH1* h, float nsig, float quant){
     float mean=h->GetMean();
     float rms=h->GetRMS();
     float upper=mean+rms*nsig;
     upper= upper - fmod(upper,quant) + quant;
     h->GetXaxis()->SetRangeUser(0,upper);
}

void set_axis_range2(TH1* h, float nsigup, float nsigdn, float quant){
     float mean=h->GetMean();
     float rms=h->GetRMS();
     float upper=mean+rms*nsigup;
     float lower=mean-rms*nsigdn;
     upper= upper - fmod(upper,quant) + quant;
     lower= lower - fmod(lower,quant) + quant;

     h->GetXaxis()->SetRangeUser(lower,upper);
}

void same_axis_range(TH1* h1, TH1* h2){
     h2->GetXaxis()->SetRange(h1->GetXaxis()->GetFirst(), h1->GetXaxis()->GetLast());
}

TH1* get_h1(const char* hname, const char* fname){
  TFile f(fname);
  if(!f.IsOpen()) return 0;
  TH1* h = dynamic_cast<TH1*>(f.Get(hname));
  if(!h) return 0;
  else h->SetDirectory(0);
  return h;
}
TH1* get_h1_from_can(const char* hname, const char* fname, const char* can){
  TFile f(fname);
  if(!f.IsOpen()) return 0;
  TCanvas* c = dynamic_cast<TCanvas*>(f.Get(can));
  if(!c) return 0;
  
  TH1* h = dynamic_cast<TH1*>(c->GetListOfPrimitives()->FindObject(hname));
  if(!h) return 0;
  else h->SetDirectory(0);
  return h;
}

TH2* get_h2(const char* hname, const char* fname){
   TFile f(fname);
   if(!f.IsOpen()) return 0;
   TH2* h = dynamic_cast<TH2*>(f.Get(hname));
   if(!h) return 0;
   else h->SetDirectory(0);
   return h;
}

TH3* get_h3(const char* hname, const char* fname){
   TFile f(fname);
   if(!f.IsOpen()) return 0;
   TH3* h = dynamic_cast<TH3*>(f.Get(hname));
   if(!h) return 0;
   else h->SetDirectory(0);
   return h;
}


TGraph* get_g(const char* hname, const char* fname){
   TFile f(fname);
   if(!f.IsOpen()) return 0;
   TGraph* h = dynamic_cast<TGraph*>(f.Get(hname));
   if(!h) return 0;
   return h;
}

TMultiGraph* get_mg(const char* hname, const char* fname){
   TFile f(fname);
   if(!f.IsOpen()) return 0;
   TMultiGraph* h = dynamic_cast<TMultiGraph*>(f.Get(hname));
   if(!h) return 0;
   return h;
}

TTree* get_t(const char* hname, const char* fname){
   TFile* f = new TFile(fname);
   if(!f->IsOpen()) return 0;
   TTree* h = dynamic_cast<TTree*>(f->Get(hname));
   if(!h) return 0;
   return h;
}

TChain* get_c(const char* hname, const char* fname){
  TChain* c = new TChain(hname,"");
  c->Add(fname);
  if(c->GetNtrees()<=0){
    std::cerr<<"get_c: could not build a chain of "<<hname<<" trees from "<<fname<<" files!"<<std::endl;
    delete c; c=0;
  }
  return c;
}

void put(TObject* obj, const char* ss, const char* file){
  if(!obj) return;
  std::string s(ss);
  std::string::size_type isep = s.find_last_of('/');
  bool make_dir=true;
  std::string name="";
  std::string dir="";
  if(isep == string::npos) {
    make_dir = false;    
    name=s;
  }
  else{
    dir.assign(s,0,isep);
    if(isep == s.size()-1) name=obj->GetName();
    else name.assign(s,isep+1,s.size()-isep);
  }
  //  TFile f(file,"update");
  //  cout<<dir<<" : "<<name<<endl;
  
  // remove slashes


  TFile f(file,"update");
  if(make_dir){
    // strip out slashes
    for(std::string::iterator j=dir.begin(); j!=dir.end(); j++){
      if(*j == '/') *j =' ';      
    }
    std::istringstream ss(dir);
    std::string tmp;
    TDirectory* wd = &f;
    while(ss.good()) {
      ss>>tmp;      
      if(! wd->cd(tmp.c_str())) { 
	std::cout<<"creating directory : *"<<tmp<<"*"<<std::endl;
	TDirectory* d =wd->mkdir(tmp.c_str()); 
	d->cd(); 
	wd=d;
      }
      else wd = gDirectory;
    }
  }
  obj->Write(name.c_str());
}



TGraph* best_fit_point(Int_t p1, Int_t p2){
  if(!gMinuit) return 0;
  TGraph* g = new TGraph(1);
  g->SetMarkerStyle(kMultiply);
  Double_t x, ex,y,ey;
  gMinuit->GetParameter(p1,x,ex);
  gMinuit->GetParameter(p2,y,ey);
  g->SetPoint(0,x,y);
  return g;
}

TGraph* two_dim_contour(Int_t p1, Int_t p2, Double_t sig, Int_t npts){
  if(!gMinuit) return 0;
  gMinuit->SetErrorDef(sig*sig);
  TGraph* g = static_cast<TGraph*>(gMinuit->Contour(npts,p1,p2));
  if(g){
    g->SetLineWidth(gStyle->GetHistLineWidth()); // TGraph() uses gStyle->GetLineWidth
  
    int N = g->GetN();
    g->SetPoint(N, g->GetX()[0],g->GetY()[0]); // make closed contour
  }
  else{
    std::cout<<"error in two_dim_contour: failure to find a contour"
	     <<std::endl;
  }
    return g;
}

TMultiGraph* get_fit_contours(Int_t p1, Int_t p2, Int_t ncont, Float_t sig_start, Float_t sig_step, Int_t npts){

  if(!gMinuit) return 0;
  TMultiGraph* a = new TMultiGraph();
  TGraph* bf = best_fit_point(p1,p2);
  a->Add(bf);
  for(int i=0; i<ncont; i++){
    double sig = sig_start + i*sig_step;
    TGraph* g = two_dim_contour(p1,p2,sig,npts);
    if(g){
      g->SetLineStyle(i+1);
      a->Add(g);
    }
    else{
      std::cout<<"error from get_fit_contours: no contour #"<<i<<std::endl;
    }
  }
  return a;
}

TMultiGraph* draw_fit_contours(Int_t p1, Int_t p2, Int_t ncont, Float_t sig_start, Float_t sig_step, Int_t npts){

  TMultiGraph* g = get_fit_contours(p1,p2,ncont,sig_start,sig_step,npts);
  g->Draw("apc");
  return g;
}

TGraph* chisquare_scan(Int_t par_num, Int_t nsig_up, Int_t nsig_down, Int_t ndiv){
  if(!gMinuit) return 0;
  if(par_num<0) return 0;
  Int_t npar = gMinuit->GetNumPars();
  if(npar==0 || par_num>(npar-1)) return 0;
  
  Double_t* pars = new Double_t[npar];
  Double_t* epars = new Double_t[npar];
  
  Double_t p,ep;
  for(int i=0; i<npar; i++){
    gMinuit->GetParameter(i,p,ep);
    pars[i]=p;
    epars[i]=ep;
  }
  
  TGraph* g = new TGraph(ndiv*2+1);
  g->SetLineWidth(gStyle->GetHistLineWidth());
  
  const Double_t x=pars[par_num];
  const Double_t ex=epars[par_num];
  
  // scan parameter from minimum to best fit
  Double_t w=ex*nsig_down/( (double)ndiv);
  Double_t y=0;
  Int_t N=0;
  for(int i=0; i<ndiv; i++){
    pars[par_num] = x-nsig_down*ex + i*w;
    gMinuit->Eval(npar,0,y,pars,3);
    g->SetPoint(N,pars[par_num],y); N++;
  }
  // best fit
  pars[par_num]=x;  
  gMinuit->Eval(npar,0,y,pars,3);
  g->SetPoint(N,pars[par_num],y); N++;
  Double_t minchi2=y;
  // scan parameter from best fit to maximum
  w=ex*nsig_up/( (double)ndiv);
  y=0;
  for(int i=0; i<ndiv; i++){
    pars[par_num] = x + (i+1)*w;
    gMinuit->Eval(npar,0,y,pars,3);
    g->SetPoint(N,pars[par_num],y); N++;
  }
  // subtract value of minimum chi2 from each point
  for(int i=0; i<N; i++){
    Double_t yy=g->GetY()[i]-minchi2;
    Double_t xx=g->GetX()[i];
    g->SetPoint(i,xx,yy);
  }

  return g;

}

TGraph* draw_chisquare_scan(Int_t par_num, Int_t nsig_up, Int_t nsig_down, Int_t ndiv){
  TGraph* g = chisquare_scan(par_num,nsig_up,nsig_down,ndiv);
  g->Draw("apc");
  return g;
}

/*
void scan_fit_2d(TVirtualFitter* vfit, TH2* h, Int_t p1, Int_t p2){
  if(!vfit || !h) return;
  
  // save current parameters 
  std::vector<double> p;
  std::vector< std::vector<double> > cov;
  for(int i=0; i<vfit->GetNumberTotalParameters(); i++){
    p.push_back(vfit->GetParameter(i));
    cov.push_back(std::vector<double>(vfit->GetNumberTotalParameters()));
    for(int j=0; j<vfit->GetNumberTotalParameters(); j++){
      cov[i][j]=vfit->GetCovarianceMatrixElement(i,j);      
    }
  }

  // loop over bins of h, fix p1 as the x center, and p2 as the y center
  // reminimize w.r.t. other parameters and fill with chi2
  
  // get some parameter information, as set by user
  double p1_low,p1_high,p2_low,p2_high, p1v,p1ev,p2v,p2ev;
  char* p1_name = new char[256];
  char* p2_name = new char[256];
  double* fit_pars= new double[vfit->GetNumberTotalParameters()];

  vfit->GetParameter(p1,p1_name,p1v,p1ev,p1_low,p1_high);
  vfit->GetParameter(p2,p2_name,p2v,p2ev,p2_low,p2_high);

  //  for(int ix=1; ix<=h->GetNbinsX(); ix++){
  for(int ix=h->GetNbinsX(); ix>=1; ix--){
    double x = h->GetXaxis()->GetBinCenter(ix);
    double ex= sqrt(fabs(cov[p1][p1]));
    vfit->SetParameter(p1,p1_name,x,ex, p1_low, p1_high);
    vfit->FixParameter(p1);
    //    for(int iy=1; iy<=h->GetNbinsY(); iy++){
    for(int iy=h->GetNbinsY(); iy>=1; iy--){
      double y = h->GetYaxis()->GetBinCenter(iy);
      double ey= sqrt(fabs(cov[p2][p2]));
      vfit->SetParameter(p2,p2_name,y,ey, p2_low, p2_high);
      vfit->FixParameter(p2);
      // minimization step
      // Oh, I just can't believe it. There is no TVirtualFitter::Fit() ?!
      // AND, TFractionFitter need its own special care...
      Int_t fit_result=-1;
      Double_t chi2=0;
      if(dynamic_cast<TFractionFitter*>(vfit->GetObjectFit())!=0){
	TFractionFitter* fit 
	  = dynamic_cast<TFractionFitter*>(vfit->GetObjectFit());
	int cntr=0;
	while(cntr<3){
	  if(fit->Fit()==0){
	    fit_result=0;
	    Double_t a,b,c; Int_t d,e;
	    vfit->GetStats(a,b,c,d,e);
	    chi2=a;
	    //	    chi2=fit->GetChisquare();
	    std::cout<<"x y chi2 = "<<x<<" "<<y<<" "<<chi2<<std::endl;
	    break;
	  }
	  else{
	    Double_t a,b,c; Int_t d,e;
	    vfit->GetStats(a,b,c,d,e);
	    chi2=a;
	    //	    chi2=fit->GetChisquare();
	    fit_result=1;
	    cntr++;
	  }
	}
      }
      else{
	if(dynamic_cast<TFitterFumili*>(vfit) !=0){
	  TFitterFumili* fit = dynamic_cast<TFitterFumili*>(vfit);
	  fit_result=fit->Minimize();
	  
	}
	else if(dynamic_cast<TFitterMinuit*>(vfit) !=0 ){
	  TFitterMinuit* fit = dynamic_cast<TFitterMinuit*>(vfit);
	  fit_result=fit->Minimize();
	}
	else if(dynamic_cast<TFitter*>(vfit) !=0){
	  TFitter* fit = dynamic_cast<TFitter*>(vfit);
	  fit_result=fit->ExecuteCommand("MINI",0,0);	
	}
	else if(dynamic_cast<TFumili*>(vfit) !=0 ){
	  TFumili* fit = dynamic_cast<TFumili*>(vfit);
	  fit_result=fit->ExecuteCommand("MINI",0,0);
	}
	else if(dynamic_cast<TLinearFitter*>(vfit) !=0){
	  TLinearFitter* fit = dynamic_cast<TLinearFitter*>(vfit);
	  fit->Eval();
	  // guess we just have to assume it converged... what shit.
	  fit_result=0;
	}
	else{
	  std::cout<<"scan_fit_2d: ERROR! Unknown derived class!"<<std::endl;
	  delete fit_pars; delete p1_name; delete p2_name; return;
	}
	
	for(int k=0; k<vfit->GetNumberTotalParameters(); k++){
	  fit_pars[k]=vfit->GetParameter(k);
	}
	chi2=vfit->Chisquare(vfit->GetNumberTotalParameters(),fit_pars);
      }
      if(fit_result!=0){
	std::cout<<"possible bad fit! for p1 = "
		 <<x<<" and p2 = "<<y<<std::endl;
	//	chi2=0;
      }
      //      h->Fill(vfit->Chisquare(vfit->GetNumberTotalParameters(),0));

      h->SetBinContent(ix,iy,chi2);
    }
  }
  delete fit_pars;
  delete p1_name;
  delete p2_name;
  return;
}
*/


void contours_on_surface(TH2* h, const char* dchi2){
  if(!h) return;
  std::istringstream in(dchi2);
  double x=0;
  int ncont=0;
  std::vector<double> v;
  while(!in.eof()){
    in>>x;
    v.push_back(h->GetMinimum()+x); ncont++;
    if(in.eof()) break;
  }
  if(ncont>0) {
    h->SetContour(ncont,&v[0]);
  }  
}

std::vector<TGraph*> get_hist_contour(int level){
  // call contours_on_surface first, draw with "CONT Z LIST" option
  // after drawing the countours for a histogram get one of them
  // saving it as a vector of graphs 
  // vector is needed since a contour may have more than one island
  
  std::vector<TGraph*> gv;
  TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

  if(level<0) return gv;
  if(level>conts->GetSize()) return gv;
  TList* contour = (TList*) conts->At(level-1);
  
  //  contour->Print();
  for(int ig=0; ig<contour->GetSize(); ig++){
    TGraph* g = (TGraph*) contour->At(ig);    
    if(g){
      //      cout<<"got graph"<<endl;
      TGraph* gcopy=(TGraph*) g->Clone();
      gv.push_back(gcopy);
    }
    else{
      //      cout<<"no graph"<<endl;
    }
  }
  
  return gv;
}


double chebyshev(int order, double x){

  const double t0=1;
  const double t1=x;
  double t=0.0;
  if(order<0) t=0;
  else if(order==0) t=t0;
  else if(order==1) t=t1;
  else{
    double tn=t1;
    double tn1=t0;
    for(int n=2; n<=order; n++){
      t=2*x*tn -tn1;
      tn1=tn;
      tn=t;
    }
  }
  return t;
}

double cheby10(double* xx, double* pp){
  const double low=pp[0];
  const double high=pp[1];
  double x=2*((*xx)-low)/(high-low) -1 ; // define x on [-1,1]
  double f=0.0;
  for(int i=0; i<10; i++){
    f+=pp[i+2]*chebyshev(i,x);
  }
  return f;

}

TH1* log_h1(const char* hname, const char* tit, int nbins, float xmin, float xmax){
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
    
  TH1* h = new TH1F(hname, tit, nbins,xbins);  
  return h;
}

TProfile* log_prof(const char* hname, const char* tit, int nbins, float xmin, float xmax){
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
    
  TProfile* h = new TProfile(hname, tit, nbins,xbins);  
  return h;
}



TH1* log10_h1(const char* hname, const char* tit, int nbins, float xmin, float xmax){
  // xbins construction taken from rootalk : P.Sizun
  
  Double_t logxmin = TMath::Log10(xmin);
  Double_t logxmax = TMath::Log10(xmax);
  Double_t binwidth = (logxmax-logxmin)/nbins;
  Double_t xbins[nbins+1];
  xbins[0] = xmin;
  for (Int_t i=1;i<=nbins;i++) {
    xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);

  }
    
  TH1* h = new TH1F(hname, tit, nbins,xbins);  
  return h;
}

TProfile* log10_prof(const char* hname, const char* tit, int nbins, float xmin, float xmax){
  // xbins construction taken from rootalk : P.Sizun
  
  Double_t logxmin = TMath::Log10(xmin);
  Double_t logxmax = TMath::Log10(xmax);
  Double_t binwidth = (logxmax-logxmin)/nbins;
  Double_t xbins[nbins+1];
  xbins[0] = xmin;
  for (Int_t i=1;i<=nbins;i++) {
    xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);

  }
    
  TProfile* h = new TProfile(hname, tit, nbins,xbins);  
  return h;
}


// draw a legend, fixing the box color and border size
// one should be able to set these via gStyle
void draw_leg(TLegend* leg){
  if(!leg) return;
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextAlign(12);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->Draw();



}

float get_pot(const char* file, int beam, const char* name){
  TH1* h = get_h1(name,file);
  if(!h){
    cout<<"get_pot: couldn't find histogram "<<name<<" in "<<file<<endl;
    return 0;
  }
  cout<<"=============== Exposure  ================="<<endl;
  cout<<"LE / LE-10 : "<<h->GetBinContent(2)<<" 10^17 POT"<<endl;
  cout<<"pME        : "<<h->GetBinContent(5)<<" 10^17 POT"<<endl;
  cout<<"pHE        : "<<h->GetBinContent(6)<<" 10^17 POT"<<endl;
  cout<<"Unknown    : "<<h->GetBinContent(1)<<" 10^17 POT"<<endl;
  cout<<"==========================================="<<endl;
  float pot = h->GetBinContent(h->FindBin(beam));
  delete h;
  return pot;
  
}

float get_pot_mc(const char* file, float scale, const char* name){
  TH1* h = get_h1(name,file);
  if(!h){
    cout<<"get_pot: couldn't find histogram "<<name<<" in "<<file<<endl;
    return 0;
  }

  cout<<"=============== Exposure  =========================="<<endl;
  cout<<"Unknown/MC : "<<h->GetBinContent(1)*(scale/1e17)<<" 10^17 POT"
      <<" ( "<<h->GetBinContent(1)<<" snarls ) "<<endl;
  cout<<"INACCURATE WHEN RUN w/ PRESELECTION ON MC TRUTH!!!!!"<<endl;
  /*
  cout<<"LE / LE-10 : "<<h->GetBinContent(2)<<" 10^17 POT"<<endl;
  cout<<"pME        : "<<h->GetBinContent(5)<<" 10^17 POT"<<endl;
  cout<<"pHE        : "<<h->GetBinContent(6)<<" 10^17 POT"<<endl;
  cout<<"Unknown    : "<<h->GetBinContent(1)<<" 10^17 POT"<<endl;
  */
  cout<<"===================================================="<<endl;
  float pot = h->GetBinContent(1)*(scale/1e17);
  delete h;
  return pot;
  
}


std::vector<TH1*> get_projections(TH2* h2, const char* opt, 
				   std::vector<int>& b1, std::vector<int>& b2,
				   const char* name_form,
				   float scale_factor){
  std::vector<TH1*> v;
  if(b1.size() != b2.size()) {
    cout<<"Error in get_projections : b1 != b2 "<<endl;
    return v;
  }
  if(h2==0){
    cout<<"Error in get_projections : h2 == 0"<<endl;
    return v;
  }
  std::string formstr=name_form;
  // deal with title in name_form
  // name_form = name ; title
  // split them into two strings
  std::string::size_type isplit = formstr.find(';');
  std::string::size_type name_end=isplit;
  //  if(isplit!=std::string::npos) name_end-=1;//avoid ;
  std::string namestr( formstr,0,name_end);
  std::string titstr = "";
  if(isplit!=std::string::npos) {
    titstr.assign(formstr,isplit+1,formstr.size()-isplit);
  }

  if(namestr=="") {namestr=h2->GetName(); namestr+="_%i";}
  cout<<"titstr: "<<titstr<<endl;
  TAxis* ax = h2->GetYaxis();
  std::string optstr=opt;
  if(optstr=="Y") ax = h2->GetXaxis();

  for(unsigned int i=0; i<b1.size(); i++){

    TH1* h=make_h2_proj(h2,opt,b1[i],b2[i],
			::Form(namestr.c_str(),gRandom->Integer(999999)),scale_factor);
    v.push_back(h);
    std::string::size_type nkey=std::string::npos;
    std::string ttit=titstr;
    
    // replace all occurences of #{binl} with value of lowbin
    while( (nkey=ttit.find("#{binl}"))!=std::string::npos ){
      std::ostringstream os; os<<b1[i];
      ttit.replace(nkey,7,os.str());
    }
    // replace all occurences of #{binh} with value of highbin
    while( (nkey=ttit.find("#{binh}"))!=std::string::npos ){
      std::ostringstream os; os<<b2[i];
      ttit.replace(nkey,7,os.str());
    }
    // replace all occurences of #{edgel} with low edge of low bin
    while( (nkey=ttit.find("#{edgel}"))!=std::string::npos ){
      std::ostringstream os; os<<ax->GetBinLowEdge(b1[i]);
      ttit.replace(nkey,8,os.str());
    }
    // replace all occurences of #{edgeh} with up endge of highest bin
    while( (nkey=ttit.find("#{edgeh}"))!=std::string::npos ){
      std::ostringstream os; os<<ax->GetBinUpEdge(b2[i]);
      ttit.replace(nkey,8,os.str());
    }
    h->SetTitle(ttit.c_str());
  }

  return v;
}


TH1* make_h2_proj(TH2* h2, const char* opt, int b1, int b2, 
		  const char* name, float scale_factor){

  TString optstr=opt;
  TH1* h=0;
  if(optstr.Contains("X",TString::kIgnoreCase)){
    h = h2->ProjectionX(name,b1,b2,"e");
  }
  else if(optstr.Contains("Y",TString::kIgnoreCase)){
    h = h2->ProjectionY(name,b1,b2,"e");
  }
  h->Scale(scale_factor);
  //  h->SetDirectory(0);
  return h;

}

TH1* get_h2_proj(const char* hname, const char* fname,const char* opt, 
		 int b1, int b2, 
		 const char* name, float scale_factor){
  TH2* h2 = get_h2(hname,fname);
  if(!h2) return 0;
  TH1* h = make_h2_proj(h2,opt,b1,b2,name,scale_factor);
  return h;
}


std::vector<TH1*> get_projections(const char* hname, const char* fname, 
				  const char* opt, 
				  std::vector<int>& b1, std::vector<int>& b2,
				  const char* name_form,
				  float scale_factor){
  TH2* h2 = get_h2(hname,fname);
  std::vector<TH1*> v;
  if(!h2) return v;
  v = get_projections(h2,opt,b1,b2,name_form,scale_factor);
  return v;

}


std::vector<TH1*> get_projections(const char* hname, const char* fname, 
				  const char* opt, 
				  std::vector<float>& b1, 
				  std::vector<float>& b2,
				  const char* name_form,
				  float scale_factor){

  TH2* h2 = get_h2(hname,fname);
  std::vector<TH1*> v;
  if(!h2) return v;
  v=get_projections(h2,opt, b1, b2, name_form, scale_factor);
  return v;
}

std::vector<TH1*> get_projections(TH2* h2, 
				  const char* opt, 
				  std::vector<float>& b1, 
				  std::vector<float>& b2,
				  const char* name_form,
				  float scale_factor){
  //  TH2* h2 = get_h2(hname,fname);
  std::vector<TH1*> v;
  if(!h2) return v;
  std::vector<int> i1;
  std::vector<int> i2;

  TString optstr=opt;
  bool isx=true;
  if(optstr.Contains("X",TString::kIgnoreCase)){
    isx=true;
  }
  else if(optstr.Contains("Y",TString::kIgnoreCase)){
    isx=false;
  }

  TAxis* ax = 0;
  if(isx) ax = h2->GetXaxis();
  else ax = h2->GetYaxis();
  for(unsigned int j = 0; j<i1.size(); j++){
    i1.push_back(ax->FindBin(b1[j]));
    i2.push_back(ax->FindBin(b2[j]));

  }
  v = get_projections(h2,opt,i1,i2,name_form,scale_factor);
  return v;

}


std::vector<TH1*> get_projections(TH2* h2,
				  const char* opt, 
				  float low, float high, int n,
				  const char* name_form,
				  float scale_factor){

  TString optstr=opt;
  bool isx=true;
  if(optstr.Contains("X",TString::kIgnoreCase)){
    isx=true;
  }
  else if(optstr.Contains("Y",TString::kIgnoreCase)){
    isx=false;
  }
  TAxis* ax = 0;
  // opt refers to projection direction
  // we want the perpendicular direction here
  if(isx) ax = h2->GetYaxis(); 
  else ax = h2->GetXaxis();
  
  std::vector<int> b1;
  std::vector<int> b2;
  
  agg_axis(ax,low,high,n,b1,b2);
  std::vector<TH1*> v=get_projections(h2,opt,b1,b2,name_form,scale_factor);
  return v;
}



std::vector<TH1*> get_projections(const char* hname, const char* fname, 
				  const char* opt, 
				  float low, float high, int n,
				  const char* name_form,
				  float scale_factor){

  TH2* h2 = get_h2(hname,fname);
  std::vector<TH1*> v;
  if(!h2) return v;
  v = get_projections(h2,opt,low,high,n,name_form,scale_factor);

  return v;

}



void agg_axis(TAxis* a, float llim, float hlim, unsigned int n,
	      std::vector<int>& low,
	      std::vector<int>& high){
  // aggregate the axis up into n mostly-equal slices
  // mostly: slop is put into largest bin.
  a->Print();
  unsigned int nbin = a->FindBin(hlim) - a->FindBin(llim) +1;
  cout<<"nbin: "<<nbin<<endl;
  unsigned int binw=1;
  if(n>0) binw=static_cast<Int_t>(floor(nbin/double(n)));
  cout<<"binw: "<<binw<<endl;
  unsigned int ll=a->FindBin(llim);
  cout<<"ll: "<<ll<<endl;
  for(unsigned int i=0; i<n; i++){
    unsigned int hl=ll+(binw-1);
    if(i+1 == n){
      // care for slop in the last bin
      if(hl!=static_cast<unsigned int>(a->FindBin(hlim))) hl=a->FindBin(hlim);
    }
    low.push_back(ll);
    high.push_back(hl);
    cout<<i<<" low bin: "<<ll<<"  high bin "<<hl
	<<" : ["<<a->GetBinLowEdge(ll)<<" , "<<a->GetBinUpEdge(hl)<<"]"<<endl;
    ll=hl+1;
  }

}


const char* choose_h3_axes(TH3* h3, const char* opt, 
			   TAxis* & a_axis, TAxis* & b_axis, TAxis* & p_axis){
  a_axis=0;
  b_axis=0;
  p_axis=0;

  TString optstr=opt;
  // encode axes XYZ
  // 0 = a_axis, 1=b_axis, -1=projection_axis
  int axis_xyz[3]={-1,-1,-1};
  // search for axis options
  int i=0;
  for(; i<optstr.Length(); i++){
    const char c = optstr(i);
    if(c=='x' || c=='X') {axis_xyz[0]=0;i++;break;}
    else if(c=='y' || c=='Y') {axis_xyz[1]=0;i++;break;}
    else if(c=='z' || c=='Z') {axis_xyz[2]=0;i++;break;}
  }
  for(; i<optstr.Length(); i++){
    const char c = optstr(i);
    if(c=='x' || c=='X') {axis_xyz[0]=1;i++;break;}
    else if(c=='y' || c=='Y') {axis_xyz[1]=1;i++;break;}
    else if(c=='z' || c=='Z') {axis_xyz[2]=1;i++;break;}
  }
  // assign the axes
  if(axis_xyz[0]==-1) p_axis=h3->GetXaxis();
  else if(axis_xyz[0]==0) a_axis=h3->GetXaxis();
  else if(axis_xyz[0]==1) b_axis=h3->GetXaxis();

  if(axis_xyz[1]==-1) p_axis=h3->GetYaxis();
  else if(axis_xyz[1]==0) a_axis=h3->GetYaxis();
  else if(axis_xyz[1]==1) b_axis=h3->GetYaxis();
  
  if(axis_xyz[2]==-1) p_axis=h3->GetZaxis();
  else if(axis_xyz[2]==0) a_axis=h3->GetZaxis();
  else if(axis_xyz[2]==1) b_axis=h3->GetZaxis();

  // find the option for TH3::Project3D

  const char* popt="";
  if(axis_xyz[0]==-1) popt="x";
  else if(axis_xyz[1]==-1) popt="y";
  else if(axis_xyz[2]==-1) popt="z";
 
  return popt;
}

// make a 1d projection of this 3d histogram
// a1,a2 and b1,b2 are bin limits
// if opt="XY" then X-->a1,a2, Y-->b1,b2 and projection shows Z
// if opt="ZX" then Z-->a1,a2, X-->b1,b2 and projection shows Y
// etc.

TH1* make_h3_proj(TH3* h3, const char* opt, int a1, int a2,
		  int b1, int b2,
		  const char* name, float scale_factor){

  TAxis* a_axis=0;
  TAxis* b_axis=0;
  TAxis* p_axis=0;
  const char* popt=choose_h3_axes(h3,opt,a_axis,b_axis,p_axis);
  
  // prepare to set axis limits in advance of projection
  // start by saving old axis limits
  const int xl_save=h3->GetXaxis()->GetFirst();
  const int xh_save=h3->GetXaxis()->GetLast();

  const int yl_save=h3->GetYaxis()->GetFirst();
  const int yh_save=h3->GetYaxis()->GetLast();  

  const int zl_save=h3->GetZaxis()->GetFirst();
  const int zh_save=h3->GetZaxis()->GetLast();

  // set axis limits, project and rename histogram
  a_axis->SetRange(a1,a2);
  b_axis->SetRange(b1,b2);
  p_axis->SetRange(-1,p_axis->GetNbins()+1);
  TH1* h = h3->Project3D(Form("e%s",popt)); // make sure errors get set  
  h->SetName(name);
  h->Scale(scale_factor);
  // reset limits
  h3->GetXaxis()->SetRange(xl_save,xh_save);
  h3->GetYaxis()->SetRange(yl_save,yh_save);
  h3->GetZaxis()->SetRange(zl_save,zh_save);

  return h;
}

std::vector<TH1*> get_h3_projections(const char* hname, const char* fname,
				     const char* opt,
				     std::vector<int>& a1,
				     std::vector<int>& a2,
				     std::vector<int>& b1,
				     std::vector<int>& b2,
				     const char* name_form,
				     float scale_factor){
  TH3* h3 = get_h3(hname,fname);

  return get_h3_projections(h3,opt,a1,a2,b1,b2,name_form,scale_factor);
  
}

std::vector<TH1*> get_h3_projections(TH3* h3, const char* opt,
				     std::vector<int>& a1,
				     std::vector<int>& a2,
				     std::vector<int>& b1,
				     std::vector<int>& b2,
				     const char* name_form,
				     float scale_factor){
  std::vector<TH1*> v;
  if(a1.size() != a2.size()) {
    cout<<"Error in get_h3_projections : a1 != a2 "<<endl;
    return v;
  }

  if(b1.size() != b2.size()) {
    cout<<"Error in get_h3_projections : b1 != b2 "<<endl;
    return v;
  }
  if(h3==0){
    cout<<"Error in get_h3_projections : h3 == 0"<<endl;
    return v;
  }
  std::string formstr=name_form;
  // deal with title in name_form
  // name_form = name ; title
  // split them into two strings
  std::string::size_type isplit = formstr.find(';');
  std::string::size_type name_end=isplit;
  //  if(isplit!=std::string::npos) name_end-=1;//avoid ;
  std::string namestr( formstr,0,name_end);
  std::string titstr = "";
  if(isplit!=std::string::npos) {
    titstr.assign(formstr,isplit+1,formstr.size()-isplit);
  }

  if(namestr=="") {namestr=h3->GetName(); namestr+="_%i";}
  cout<<"titstr: "<<titstr<<endl;
  /*
  TAxis* ax = h2->GetYaxis();
  std::string optstr=opt;
  if(optstr=="Y") ax = h2->GetXaxis();
  */
  TAxis* a_axis;
  TAxis* b_axis;
  TAxis* p_axis;
  choose_h3_axes(h3,opt,a_axis,b_axis,p_axis);
  
  for(unsigned int i=0; i<b1.size(); i++){

    TH1* h=make_h3_proj(h3,opt,a1[i],a2[i],b1[i],b2[i],
			::Form(namestr.c_str(),gRandom->Integer(999999)),scale_factor);
    v.push_back(h);
    
    // format title
    std::string::size_type nkey=std::string::npos;
    std::string ttit=titstr;
    
    // replace all occurences of #{abinl} with value of a_axis lowbin
    while( (nkey=ttit.find("#{abinl}"))!=std::string::npos ){
      std::ostringstream os; os<<a1[i];
      ttit.replace(nkey,8,os.str());
    }
    // replace all occurences of #{abinh} with value of a_axis highbin
    while( (nkey=ttit.find("#{abinh}"))!=std::string::npos ){
      std::ostringstream os; os<<a2[i];
      ttit.replace(nkey,8,os.str());
    }

    // replace all occurences of #{binl} with value of lowbin
    while( (nkey=ttit.find("#{bbinl}"))!=std::string::npos ){
      std::ostringstream os; os<<b1[i];
      ttit.replace(nkey,8,os.str());
    }
    // replace all occurences of #{binh} with value of highbin
    while( (nkey=ttit.find("#{bbinh}"))!=std::string::npos ){
      std::ostringstream os; os<<b2[i];
      ttit.replace(nkey,8,os.str());
    }

    // replace all occurences of #{aedgel} with low edge of low bin
    while( (nkey=ttit.find("#{aedgel}"))!=std::string::npos ){
      std::ostringstream os; os<<a_axis->GetBinLowEdge(a1[i]);
      ttit.replace(nkey,9,os.str());
    }
    // replace all occurences of #{adgeh} with up endge of highest bin
    while( (nkey=ttit.find("#{aedgeh}"))!=std::string::npos ){
      std::ostringstream os; os<<a_axis->GetBinUpEdge(a2[i]);
      ttit.replace(nkey,9,os.str());
    }

    // replace all occurences of #{bedgel} with low edge of low bin
    while( (nkey=ttit.find("#{bedgel}"))!=std::string::npos ){
      std::ostringstream os; os<<b_axis->GetBinLowEdge(b1[i]);
      ttit.replace(nkey,9,os.str());
    }
    // replace all occurences of #{bdgeh} with up endge of highest bin
    while( (nkey=ttit.find("#{bedgeh}"))!=std::string::npos ){
      std::ostringstream os; os<<b_axis->GetBinUpEdge(b2[i]);
      ttit.replace(nkey,9,os.str());
    }
    h->SetTitle(ttit.c_str());
  }

  return v;
  
}

std::vector<TH1*> get_h3_projections(TH3* h3, const char* opt,
				     TH2* h2,
				     const char* name_form,
				     float scale_factor){
  std::vector<TH1*> v;
  //
  // project the 3rd dimension of h3, using bins of h2

  TAxis* a_axis;
  TAxis* b_axis;
  TAxis* p_axis;
  choose_h3_axes(h3,opt,a_axis,b_axis,p_axis);

  std::cout<<"a_axis: ";
  a_axis->Print();
  std::cout<<std::endl;
  std::cout<<"b_axis: ";
  b_axis->Print();
  std::cout<<std::endl;

  const Axis_t epsilon=1e-9; // not well programmed...
  for(int ix=1; ix<=h2->GetXaxis()->GetNbins(); ix++){
      const Axis_t xl= h2->GetXaxis()->GetBinLowEdge(ix);
      const Axis_t xh= h2->GetXaxis()->GetBinUpEdge(ix);
      const Int_t al=a_axis->FindBin(xl+epsilon);
      const Int_t ah=a_axis->FindBin(xh-epsilon);
    for(int iy=1; iy<=h2->GetYaxis()->GetNbins(); iy++){
      const Axis_t yl= h2->GetYaxis()->GetBinLowEdge(iy);
      const Axis_t yh= h2->GetYaxis()->GetBinUpEdge(iy);
      const Int_t bl=b_axis->FindBin(yl+epsilon);
      const Int_t bh=b_axis->FindBin(yh-epsilon);
      
      //      a_axis->SetRange(al,ah);
      //      b_axis->SetRange(bl,bh);


      std::string formstr=name_form;
      // deal with title in name_form
      // name_form = name ; title
      // split them into two strings
      std::string::size_type isplit = formstr.find(';');
      std::string::size_type name_end=isplit;
      //  if(isplit!=std::string::npos) name_end-=1;//avoid ;
      std::string namestr( formstr,0,name_end);
      std::string titstr = "";
      if(isplit!=std::string::npos) {
	titstr.assign(formstr,isplit+1,formstr.size()-isplit);
      }
      
      if(namestr=="") {namestr=h3->GetName(); namestr+="_%i";}
      //      cout<<"titstr: "<<titstr<<endl;
      // format title
      std::string::size_type nkey=std::string::npos;
      std::string ttit=titstr;

      TH1* h1 = make_h3_proj(h3,opt,al,ah,bl,bh,
			     ::Form(namestr.c_str(),gRandom->Integer(999999)),
			     scale_factor);
      
      v.push_back(h1);

      
      // replace all occurences of #{abinl} with value of a_axis lowbin
      while( (nkey=ttit.find("#{abinl}"))!=std::string::npos ){
	std::ostringstream os; os<<ah;
	ttit.replace(nkey,8,os.str());
      }
      // replace all occurences of #{abinh} with value of a_axis highbin
      while( (nkey=ttit.find("#{abinh}"))!=std::string::npos ){
	std::ostringstream os; os<<al;
	ttit.replace(nkey,8,os.str());
      }
      
      // replace all occurences of #{binl} with value of lowbin
      while( (nkey=ttit.find("#{bbinl}"))!=std::string::npos ){
	std::ostringstream os; os<<bl;
	ttit.replace(nkey,8,os.str());
      }
      // replace all occurences of #{binh} with value of highbin
      while( (nkey=ttit.find("#{bbinh}"))!=std::string::npos ){
	std::ostringstream os; os<<bh;
	ttit.replace(nkey,8,os.str());
      }
      
      // replace all occurences of #{aedgel} with low edge of low bin
      while( (nkey=ttit.find("#{aedgel}"))!=std::string::npos ){
	std::ostringstream os; os<<a_axis->GetBinLowEdge(al);
	ttit.replace(nkey,9,os.str());
      }
      // replace all occurences of #{adgeh} with up endge of highest bin
      while( (nkey=ttit.find("#{aedgeh}"))!=std::string::npos ){
	std::ostringstream os; os<<a_axis->GetBinUpEdge(ah);
	ttit.replace(nkey,9,os.str());
      }
      
      // replace all occurences of #{bedgel} with low edge of low bin
      while( (nkey=ttit.find("#{bedgel}"))!=std::string::npos ){
	std::ostringstream os; os<<b_axis->GetBinLowEdge(bl);
	ttit.replace(nkey,9,os.str());
      }
      // replace all occurences of #{bdgeh} with up endge of highest bin
      while( (nkey=ttit.find("#{bedgeh}"))!=std::string::npos ){
	std::ostringstream os; os<<b_axis->GetBinUpEdge(bh);
	ttit.replace(nkey,9,os.str());
      }
      h1->SetTitle(ttit.c_str());      

      std::cout<<ix<<","<<iy<<" : "<<ttit<<std::endl;

    }
  }
  return v;
}




void plot_many(std::vector<TH1*>& v, const std::vector<TPad*>& pads, 
	       const char* dopt){
  unsigned int npad=pads.size();
  for(unsigned int i=0; i<v.size(); i++){
    if(npad<=i) continue;
    TH1* h = v[i];
    pads[i]->cd();
    if(h){
      h->Draw(dopt);
      TIter fiter(h->GetListOfFunctions());
      TF1* func = 0;
      while( (func = dynamic_cast<TF1*>(fiter.Next()) ) )  {
	func->Draw("lsame");
      }
    }
  }
}

void plot_many(std::vector<TH1*>& v, int nh, int nv, 
	       const char* xtit, const char* ytit, const char* opt,
	       std::vector<TCanvas*>& cans, std::vector<TPad*>& pads,
	       const char* cantit,
	       float xmargin, float ymargin,
	       float xsubmargin,float ysubmargin){

  const int npad=nh*nv;
  int can_num=0;
  TCanvas* c=0;
  TPad* p1 = 0;
  for(unsigned int i=0; i<v.size(); i++){
    if(i%npad == 0){// new canvas
      TString s=cantit; s+="_"; s+=can_num;
      // control canvas size via gStyle settings
      c = new TCanvas(s.Data(),s.Data());
      cans.push_back(c);
      TString s2=cantit; s2+="_pad";
      p1 = new TPad(s2.Data(),s2.Data(),xmargin,ymargin,1.0,1.0,0,0,0);
      p1->Draw();
      p1->Divide(nh,nv,xsubmargin,ysubmargin);
      // loop over subpads and add to array of pads
      TIter itr(p1->GetListOfPrimitives());
      while(TObject* obj = itr.Next()){
	TPad* p = dynamic_cast<TPad*>(obj);
	if(p) pads.push_back(p);
      }
      c->cd(0);
      TLatex lx;
      lx.SetNDC();
      lx.SetTextAlign(22);
      lx.SetTextFont(gStyle->GetTitleFont());
      lx.SetTextSize(ymargin*(0.05/0.07));
      if(xtit!=""){
	lx.DrawLatex(0.5,ymargin/2.0,xtit);
      }
      if(ytit!=""){
	lx.SetTextAngle(90);
	lx.DrawLatex(xmargin/2.0,0.5,ytit);
      }
      can_num++;
    }
    int curpad=(i%npad)+1;
    //    cout<<"curpad: "<<curpad<<endl;
    p1->cd(curpad);
    TH1* h = v[i];
    if(h) {
      h->Draw(opt);
      // look for functions to draw
      TIter fiter(h->GetListOfFunctions());
      TF1* func = 0;
      while( (func = dynamic_cast<TF1*>(fiter.Next()) ) )  {
	func->Draw("lsame");
      }      
    }
  }

}

TTree* filter_tree(TTree* t, const char* selection, const char* vars, TFile* fout)
{

  // ttf implements the event selection
  TTreeFormula ttf("filter_tree",selection,t);

  // turn off all branches in input
  t->SetBranchStatus("*",0);

  // tokenize array of variable names
  TString s(vars);
  TObjArray* vlist = s.Tokenize(", ");
  // loop over variable names, turning on associated branch in input
  std::cout<<"filter_tree: writing out the following variables:"<<std::endl;
  for(int i=0; i<vlist->GetEntries(); i++){
    TObjString* obj = dynamic_cast<TObjString*>(vlist->At(i));
    TString v=obj->GetString();
    std::cout<<v<<std::endl;
    t->SetBranchStatus(v.Data(),1);
  }
  std::cout<<"end of variables"<<std::endl;

  fout->cd();
  TTree* tc = t->CloneTree(0);

  // go through and turn on branches needed by selection  
  std::cout<<"Turning on branches"<<std::endl;
  for(Int_t j=0; j< ttf.GetNcodes(); j++){
    std::cout<<"Enabling branch : "<<ttf.GetLeaf(j)->GetBranch()->GetName()<<std::endl;
    t->SetBranchStatus(ttf.GetLeaf(j)->GetBranch()->GetName(),1);
  }
  ttf.UpdateFormulaLeaves();
  //  const float tick_frac=0.02;
  //  float done_frac=0;
  Long64_t nbr=0;
  std::cout.precision(3);
  Long64_t i=0;
  Long64_t npass=0;
  bool keep_going=true;
  Int_t current_tree=-1;

  while(true) {
    Long64_t local_entry=t->LoadTree(i);
    if(local_entry<0) {
      std::cout<<"filter_tree: LoadTree()<0: finishing"<<std::endl;
      break;
    }
    if(t->GetTreeNumber()!=current_tree) {
      current_tree=t->GetTreeNumber(); 
      ttf.UpdateFormulaLeaves();
    }
    // load branches needed by selection
    
    for(Int_t j=0; j< ttf.GetNcodes(); j++){

      Long64_t nb = ttf.GetLeaf(j)->GetBranch()->GetEntry(local_entry);
      /*
      if(nb<=0){
	// done reading
	std::cout<<"Should be done (1)? : [i,j] = "<<i<<" , "<<j<<std::endl;
	std::cout<<"Branch giving 0: "<<ttf.GetLeaf(j)->GetBranch()->GetName()<<std::endl;
	std::cout<<"local_entry: "<<local_entry<<std::endl;
	keep_going=false;
	break;
      }
      nbr+=nb;
      */
    }
    if(!keep_going) break;
    
    // evaluate the selection, if true, fill output
    Int_t ndata = ttf.GetNdata();
    for( Int_t idata=0; idata<ndata; idata++){ // idata>1 in arrays
      Double_t val = ttf.EvalInstance(idata);
      if(val > 0.0){
	npass++;
	Long64_t nb=t->GetEntry(i,0);
	if(nb<=0){
	  // done reading
	  std::cout<<"Should be done (2)?"<<std::endl;
	  keep_going=false;
	  continue;
	}
	nbr+=nb;      
	tc->Fill();
	// if we are looping over an array only write one record
	// even if multiple array elements pass
	break; 
      }
    }
    // how's it going
    if(i%100000==0){
      std::cout<<i<<" entries  read : "<<npass<<" filled : "
      <<nbr/(1024.0*1024.0)<<" Mb read"<<std::endl;
    }
    
    i++; // don't forget
  }
  std::cout<<"Done                                              "<<std::endl;
  std::cout<<"filter_tree: ending on tree # "<<current_tree<<std::endl;
  std::cout<<"filter_tree: output has "
	   <<npass<<" entries."<<std::endl;  
  return tc;
}

/*
TTree* filter_tree(TTree* t, const char* selection, const char* vars, TFile* fout)
//void ttt(const char* vars);
//void ttt(const char* vars)
{
  // select events
  std::cout<<" filter_tree: selecting entries which satisfy\n( "<<selection<<" )\n"<<std::endl;
  TEventList* el = new TEventList("el", "selection");
  t->Draw(">>el",selection,"goff");
  t->SetEventList(el);
  const Long64_t nent=el->GetN();
  std::cout<<"filter_tree: "<<nent<<" entries passed selection."<<std::endl;
  // turn off all branches
  t->SetBranchStatus("*",0);

  // tokenize array of variable names
  TString s(vars);
  TObjArray* vlist = s.Tokenize(", ");
  // loop over variable names, turning on associated branch
  std::cout<<"filter_tree: writing out the following variables:"<<std::endl;
  for(int i=0; i<vlist->GetEntries(); i++){
    TObjString* obj = dynamic_cast<TObjString*>(vlist->At(i));
    TString v=obj->GetString();
    std::cout<<v<<std::endl;
    t->SetBranchStatus(v.Data(),1);
  }
  std::cout<<"end of variables"<<std::endl;
  //  TTree* tc = t->CopyTree("","");
  fout->cd();
  TTree* tc = t->CloneTree(0);
  //  tc->CopyEntries(t);

  const float tick_frac=0.02;
  float done_frac=0;
  //  const Long64_t tick_evts=static_cast<Long64_t>(tick_frac*nent);
  //  std::cout<<"tick_evts: "<<tick_evts<<std::endl;
  std::cout.precision(3);
  for (Int_t i=0;i<nent; i++) {
    t->GetEntry(el->GetEntry(i));
    tc->Fill();
    if(i/double(nent)>done_frac){
      done_frac+=tick_frac;
      std::cout<<"progress: "<<(i/double(nent))*100.0<<" %              \r"; 
      std::cout.flush();
    }
  }
  std::cout<<"Done                                              "<<std::endl;
  std::cout<<"filter_tree: output has "
	   <<tc->GetEntries()<<" entries."<<std::endl;  
  return tc;
}
*/

void filter_tree(const char* ntuple, const char* files, const char* selection, const char* vars, const char* fout){
  std::cout<<"filter_tree: building chain of '"
	   <<ntuple<<"' ntuples from "<<files<<std::endl;
  TChain c(ntuple,ntuple);
  c.Add(files);
  /*
  std::cout<<"filter_tree: using "<<c.GetNtrees()
	   <<" trees with a total of "
	   <<c.GetEntries()<<" entries."<<std::endl;
  */

  std::cout<<"filter_tree: using "<<c.GetNtrees()
	   <<" trees."<<std::endl;


  TFile f(fout,"recreate");
  filter_tree(&c,selection,vars,&f);
  
  f.Write();

}


void set_lc(std::vector<TH1*>& v, Color_t x){
  for(unsigned int i=0; i<v.size(); i++) v[i]->SetLineColor(x);
}
void set_fc(std::vector<TH1*>& v, Color_t x){
  for(unsigned int i=0; i<v.size(); i++) v[i]->SetFillColor(x);
}
void set_mc(std::vector<TH1*>& v, Color_t x){
  for(unsigned int i=0; i<v.size(); i++) v[i]->SetMarkerColor(x);
}
void set_ls(std::vector<TH1*>& v, Style_t x){
  for(unsigned int i=0; i<v.size(); i++) v[i]->SetLineStyle(x);
}
void set_fs(std::vector<TH1*>& v, Style_t x){
  for(unsigned int i=0; i<v.size(); i++) v[i]->SetFillStyle(x);
}
void set_ms(std::vector<TH1*>& v, Style_t x){
  for(unsigned int i=0; i<v.size(); i++) v[i]->SetMarkerStyle(x);
}

void print_canvases(std::vector<TCanvas*>& v, const char* fmt,
		const char* extensions){
  for(unsigned int i=0; i<v.size(); i++){ 
    std::string s=Form(fmt,i);
    print_canvas(v[i],s.c_str(),extensions);} 
  
}

void print_canvas(TCanvas* can,const char* name,const char* extensions){
  if(!can) return;
  TVirtualPad* save = gPad;
  can->cd(0);
  
  // tokenize extensions
  TString s(extensions);
  // they can be seperated by , or a space
  TObjArray* vlist = s.Tokenize(", ");
  // for each extension, go through and print the canvas
  for(int i=0; i<vlist->GetEntries(); i++){

    TObjString* obj = dynamic_cast<TObjString*>(vlist->At(i));
    // v is the extension
    TString v=obj->GetString();
    // remove leading '.' if one is there
    // prevents file..eps sort of names
    v.Remove(TString::kLeading,'.');
    // add name and '.' to extension
    TString v2(name);
    v2+="."; v2+=v;
    can->Print(v2.Data());
  }
  if(save) save->cd();
  return;
  
}

void print_canvas(const char* cname,const char* name,const char* extensions){
  TCanvas* c= 
    dynamic_cast<TCanvas*>(gROOT->GetListOfCanvases()->FindObject(cname));
  if(c) print_canvas(c,name,extensions);
}


TH1* rebin_like(TH1* input, TH1* likethis, const char* newname){
  if(!input) return 0;
  if(!likethis) return 0;
  // should check that input can actually be unambiguously rebinned likethis
  // do this by looping over the new histogram and checking that the bin edges match with those in the input
  const Double_t epsilon = 1e-10;
  for(int i=1; i<=likethis->GetNbinsX(); i++){
    Double_t ble_likethis = likethis->GetBinLowEdge(i);
    bool found = false;
    for(int j=1; j<=input->GetNbinsX(); j++){
      Double_t ble_input = input->GetBinLowEdge(j);
      if(fabs(ble_likethis-ble_input)<epsilon) {found=true; break;}
    }
    if(!found) {
      std::cout<<"rebin_like: bin edges of input histogram do not agree with bins of likethis histogram"<<std::endl;
      return 0;
    }    
  }
  // made it to here, now it's easy
  //  TH1* htemp = static_cast<TH1*>(gDirectory->Get(newname));
  //  if(htemp){htemp->SetDirectory(0); delete htemp;}

  TH1* h = static_cast<TH1*>(likethis->Clone(newname));
  h->Clear();
  h->Reset();

  for(int i=1; i<input->GetNbinsX();i++){
    Axis_t cen = input->GetBinCenter(i);
    Stat_t cont = input->GetBinContent(i);
    Stat_t err = input->GetBinError(i);
    Int_t bin = h->FindBin(cen);
    Stat_t obc = h->GetBinContent(bin);
    Stat_t obe = h->GetBinError(bin);
    h->SetBinContent(bin,cont+obc);
    h->SetBinError(bin,sqrt(obe*obe + err*err));
  }
 
  //  h->SetDirectory(0);
  return h;

}


TH2* rebin_like(TH2* input, TH2* likethis, const char* newname){
  if(!input) return 0;
  if(!likethis) return 0;
  // should check that input can actually be unambiguously rebinned likethis
  // do this by looping over the new histogram and checking that the bin edges match with those in the input
  const Double_t epsilon = 1e-10;
  for(int i=1; i<=likethis->GetNbinsX(); i++){
    Double_t ble_likethis = likethis->GetXaxis()->GetBinLowEdge(i);
    bool found = false;
    for(int j=1; j<=input->GetNbinsX(); j++){
      Double_t ble_input = input->GetXaxis()->GetBinLowEdge(j);
      if(fabs(ble_likethis-ble_input)<epsilon) {found=true; break;}
    }
    if(!found) {
      std::cout<<"rebin_like: X axis bin edges of input histogram do not agree with bins of likethis histogram"<<std::endl;
      return 0;
    }    
  }

  for(int i=1; i<=likethis->GetNbinsY(); i++){
    Double_t ble_likethis = likethis->GetYaxis()->GetBinLowEdge(i);
    bool found = false;
    for(int j=1; j<=input->GetNbinsY(); j++){
      Double_t ble_input = input->GetYaxis()->GetBinLowEdge(j);
      if(fabs(ble_likethis-ble_input)<epsilon) {found=true; break;}
    }
    if(!found) {
      std::cout<<"rebin_like: Y axis bin edges of input histogram do not agree with bins of likethis histogram"<<std::endl;
      return 0;
    }    
  }

  // made it to here, now it's easy
  //  TH1* htemp = static_cast<TH1*>(gDirectory->Get(newname));
  //  if(htemp){htemp->SetDirectory(0); delete htemp;}

  TH2* h = static_cast<TH2*>(likethis->Clone(newname));
  h->Clear();
  h->Reset();

  for(int i=1; i<input->GetNbinsX();i++){
    for(int j=1; j<input->GetNbinsY();j++){
      Axis_t xcen = input->GetXaxis()->GetBinCenter(i);
      Axis_t ycen = input->GetYaxis()->GetBinCenter(j);
      Stat_t cont = input->GetBinContent(i,j);
      Stat_t err = input->GetBinError(i,j);
      Int_t bin = h->FindBin(xcen,ycen);
      Stat_t obc = h->GetBinContent(bin);
      Stat_t obe = h->GetBinError(bin);
      h->SetBinContent(bin,cont+obc);
      h->SetBinError(bin,sqrt(obe*obe + err*err));
    }
  }
  //  h->SetDirectory(0);
  return h;

}

Int_t make_binning(const char* binning, const char* opt, std::vector<Double_t>& low_edges, std::vector<Double_t>& high_edges){
  
  // make a TH1F according to binning
  // here is how it works:
  // opt="W": low_limit, nbins binw, nbins binw, etc.
  // opt="R": low_limit, nbins high_limit, nbins high_limit, etc.

  // clear out bin edge vectors 
  low_edges.clear();
  high_edges.clear();

  // first, tokenise on commas
  TString s(binning);
  TObjArray* vlist = s.Tokenize(",");

  // must have at least two entries
  assert(vlist->GetEntries()>1);
  
  Double_t low_limit = atof(dynamic_cast<TObjString*>(vlist->At(0))->GetString().Data());
  
  std::vector<Int_t> vn;
  std::vector<Double_t> vx;
  // for each binning set, place the nbins and floating value in arrays
  for(int i=1; i<vlist->GetEntries(); i++){
    istringstream iss(dynamic_cast<TObjString*>(vlist->At(i))->GetString().Data());
    Int_t n; Double_t x;
    iss>>n>>x;
    vn.push_back(n);
    vx.push_back(x);
  }
  
  // figure out which option we have
  TString sopt(opt);
  Int_t iopt=0;
  if(sopt.Contains("W")) iopt=1;
  else if(sopt.Contains("R")) iopt=2;
  assert(iopt);

  // loop over the arrays to fill bin edges vectors
  for(unsigned int i=0; i<vn.size(); i++){
    // figure out the bin width for this binning
    Double_t binw=0.0;
    // first time through start on the lower limit
    Double_t start_edge = low_limit;
    // successive times, start at the end of the high_edges array
    if(i>0) start_edge = high_edges[ high_edges.size()-1 ];

    if(iopt==1){
      binw=vx[i];
    }
    else if(iopt==2){
      binw=(vx[i]-start_edge)/vn[i];
    }
    // now, calculate bin edges and push into edges arrays
    for(int j=0; j<vn[i]; j++){
      Double_t low = start_edge + binw*j;
      Double_t high = low+binw;
      low_edges.push_back(low);
      high_edges.push_back(high);
    }
  }
  for(unsigned int i=0; i<low_edges.size(); i++){
    //    std::cout<<i<<" [ "<<low_edges[i]<<" , "<<high_edges[i]<<" ] "<<std::endl;

  }

  Int_t result = low_edges.size();
  return result;
}


TH1* make_h1(const char* name, const char* title, const char* binning, const char* opt ){

  TH1* result=0;
  // make a TH1F according to binning
  // see make_binning() for a description of the binning string and opt
  std::vector<Double_t> vlow;
  std::vector<Double_t> vhigh;
  Int_t nbins = make_binning(binning,opt,vlow,vhigh);
  assert(nbins>0);
  vlow.push_back(vhigh[ vhigh.size() -1 ]);
  result = new TH1F(name,title,nbins,&vlow[0]);
  return result;
}

TH2* make_h2(const char* name, const char* title, const char* binning_x, const char* opt_x, const char* binning_y, const char* opt_y){
  
  TH2* result=0;
  // make a TH2F according to binning
  // see make_binning() for a description of the binning string and opt
  std::vector<Double_t> xlow;
  std::vector<Double_t> xhigh;
  Int_t nbins_x = make_binning(binning_x,opt_x,xlow,xhigh);
  assert(nbins_x>0);
  xlow.push_back(xhigh[ xhigh.size() -1 ]);

  std::vector<Double_t> ylow;
  std::vector<Double_t> yhigh;
  Int_t nbins_y = make_binning(binning_y,opt_y,ylow,yhigh);
  assert(nbins_y>0);
  ylow.push_back(yhigh[ yhigh.size() -1 ]);


  result = new TH2F(name,title,nbins_x,&xlow[0],nbins_y,&ylow[0]);
  return result;

}


void divide_binw(TH1* h, float m){
  // divide each bin by its width
  if(!h) return;
  for(int i=1; i<=h->GetNbinsX(); i++){
    float bw = h->GetBinWidth(i);
    h->SetBinContent(i,h->GetBinContent(i)*m/bw);
    h->SetBinError(i,h->GetBinError(i)*m/bw);
  }

}

void divide_binw_2d(TH2* h, float m){
  // divide each bin by its width
  if(!h) return;
  for(int i=1; i<=h->GetNbinsX(); i++){
    for(int j=1; j<=h->GetNbinsY(); j++){
      float bw = h->GetXaxis()->GetBinWidth(i)*h->GetYaxis()->GetBinWidth(j);
      h->SetBinContent(i,j,h->GetBinContent(i,j)*m/bw);
      h->SetBinError(i,j,h->GetBinError(i,j)*m/bw);
    }
  }

}

// shift/translate/move a TH1
TH1* translate_h1(TH1* h, double t, const char* hname){
  // from Roottalk forum: O. Couet
   Double_t xmin = h->GetXaxis()->GetXmin()+t;
   Double_t xmax = h->GetXaxis()->GetXmax()+t;
   Int_t n = h->GetNbinsX();
   TH1 *ht = new TH1F(Form(hname,h->GetName()), h->GetTitle(), n, xmin, xmax);
   for (int i=1; i<=n; i++) {
      ht->SetBinContent(i,h->GetBinContent(i));
      ht->SetBinError(i,h->GetBinError(i));
   }
   ht->SetEntries(h->GetEntries());
   return ht;


}



void merge_hists(TDirectory* target, TDirectory* from){
  // loop over all histograms in 'target'
  // look for equivalent in 'from'
  // TH1::Add()
  if(!target || !from) return;
  TIter it(target->GetList());
  
  TObject* obj =0;
  while( (obj = it())){
    TH1* htarg = dynamic_cast<TH1*>(obj);
    if(!htarg) continue;
    TH1* hfrom = dynamic_cast<TH1*>(from->Get(htarg->GetName()));
    if(!hfrom) continue;
    htarg->Add(hfrom,1.0);    
  }
  return;
}
void copy_hists(TDirectory* from, TDirectory* target){
  if(!target || !from) return;
  TIter it(from->GetListOfKeys());
  TKey* key=0;
  while( (key = static_cast<TKey*>(it())) ){
    TH1* h = dynamic_cast<TH1*>(key->ReadObj());
    if(h) h->SetDirectory(target);
  }
  return;
}

/*

// had to comment out myhadd for root 5.34, due to errors compiling move_stats
// example:
// Info in <TUnixSystem::ACLiC>: creating shared library /home/kordosky/root_build/root_5_34_28_libs//home/kordosky/base_macros/move_stats_C.so
//Error: Symbol __size_t gl_pathc is not defined in current scope  /usr/include/glob.h:106:


void myhadd(TDirectory* target, const char* infiles){
  if(!target) return;
  glob_t found_files;
  int flags = 0;
  int n = glob(infiles,flags,NULL,&found_files);
  if(n!=0) {
    std::cout<<"myhadd: glob() returns "
	     <<n<<". This is bad!"<<std::endl; 
    globfree(&found_files);
    return;
  }
  if(found_files.gl_pathc == 0){
    std::cout<<"myhadd: found no input files"<<std::cout;
    return;
  }
  // for first file, make a copy of histograms
  {
    TFile fin0(found_files.gl_pathv[0]);
    copy_hists(&fin0,target);
  }
  // for the other files, merge
  for(size_t i=1; i<found_files.gl_pathc; i++){
    std::cout<<"Adding: "<<found_files.gl_pathv[i]<<std::endl;
    TFile fin(found_files.gl_pathv[i]);
    merge_hists(target,&fin);
  }
  globfree(&found_files);
  
  return;
}


void myhadd(const char* outfile, const char* infiles){
  TFile f(outfile,"recreate");
  myhadd(&f,infiles);
  f.Write();
}
*/


void print_exposure(TDirectory* dir){
  if(!dir) return;
  TH1* h_tortgt = dynamic_cast<TH1*>(dir->Get("h_tortgt"));
  TH1* h_tormi = dynamic_cast<TH1*>(dir->Get("h_tormi"));
  TH1* h_snarls = dynamic_cast<TH1*>(dir->Get("h_snarls"));

  cout<<"=============== Exposure  TORTGT : TOR101 ================="<<endl;
  cout<<"LE / LE-10 : "<<( h_tortgt ? h_tortgt->GetBinContent(2) : 0)<<" : "<<( h_tormi ? h_tormi->GetBinContent(2) : 0)<<" 10^17 POT"<<endl;
  cout<<"pME        : "<<( h_tortgt ? h_tortgt->GetBinContent(5) : 0)<<" : "<<( h_tormi ? h_tormi->GetBinContent(5) : 0)<<" 10^17 POT"<<endl;
  cout<<"pHE        : "<<( h_tortgt ? h_tortgt->GetBinContent(6) : 0)<<" : "<<( h_tormi ? h_tormi->GetBinContent(6) : 0)<<" 10^17 POT"<<endl;
  cout<<"Unknown    : "<<( h_tortgt ? h_tortgt->GetBinContent(1) : 0)<<" : "<<( h_tormi ? h_tormi->GetBinContent(1) : 0)<<" 10^17 POT"<<endl;
  cout<<"==========================================================="<<endl;
  cout<<"Snarls passing pre-selection: "<<(h_snarls ? h_snarls->GetEntries() : 0)<<endl;
  cout<<"POT for 2.4e13 pot/snarl : "
      <<(h_snarls ? h_snarls->GetEntries()*2.42e-4 : 0)<<" 10^17 POT"<<endl;

  return;
}

double chi2_h1(TH1* h1, TH1* h2, const char* opt){
  double chi2; int N;
  return chi2_h1(h1,h2,chi2,N,opt);
}
double chi2_h1(TH1* h1, TH1* h2, double& chi2, int& N, const char* opt){
  // direct implementation of numerical recipes calculation
  //
  double sum1=0.0;
  double sum2=0.0;
  N=0;
  int first=h1->GetXaxis()->GetFirst();
  int last=h1->GetXaxis()->GetLast();
  for(int i=first; i<=last; i++){
    N++;
    sum1+=h1->GetBinContent(i);
    sum2+=h2->GetBinContent(i);
  }
  chi2=0.0;
  for(int i=first; i<=last; i++){
    chi2+= pow(sqrt(sum2/sum1)*h1->GetBinContent(i)-sqrt(sum1/sum2)*h2->GetBinContent(i),2.0)/(pow(h1->GetBinError(i),2.0)+pow(h2->GetBinError(i),2.0));
  }
  N-=1;
  return TMath::Prob(0.5*chi2, static_cast<int>(floor(0.5*N+0.001)));

}

//////////////// utilities for I/O //////////////////////////////

/*
// removed due to compilation problems with root 5.34
// see note in myhadd above
void glob_files(const char* infiles, std::vector<std::string>& v){

  glob_t found_files;
  int flags = 0;
  int n = glob(infiles,flags,NULL,&found_files);
  if(n!=0) {
    std::cerr<<"glob_files: glob() returns "
	     <<n<<". This is bad!"<<std::endl; 
    globfree(&found_files);
    return;
  }
  if(found_files.gl_pathc == 0){
    std::cerr<<"glob_vfiles: found no input files"<<std::cout;
    return;
  }
  // for the other files, merge
  for(size_t i=0; i<found_files.gl_pathc; i++){
    v.push_back(found_files.gl_pathv[i]);
  }
  globfree(&found_files);
  
  return;
}
*/

void write_parameter(const char* name, double val, TDirectory* dir){
  TDirectory* d = gDirectory;
  if(dir){ dir->cd(); }

  TParameter<double> p(name,val);
  p.Write(name);
  if(d) d->cd();
}

//void read_parameter(const char* name, double& val, TDirectory* dir=0)

void pars_from_ascii(const char* filename, std::vector<TParameter<double> >& v,
		     const char* begin_key, const char* end_key){

  std::ifstream f(filename);
  std::string s;
  const std::string bk(begin_key);
  const std::string ek(end_key);
  while(!f.eof()){
    std::getline(f,s);
    std::cout<<s<<std::endl;
    if(s==bk && !f.eof()){// found begin key
      std::getline(f,s);
      while(s!=ek && !f.eof()){
	// inside of data area
	std::istringstream iss(s);
	std::string name;
	double val;
	iss>>name>>val;
	TParameter<double> p(name.c_str(),val);
	std::cout<<name<<" : "<<val<<std::endl;
	v.push_back(p);
	std::getline(f,s);
      }
    
    }

  }
  
}

void pars_from_dir(TDirectory* d, std::vector<TParameter<double> >& v){
  
  if(!d) return;
  TIter it(d->GetListOfKeys());
  //  TObject* obj=0;
  TKey* key=0;
  while((key=static_cast<TKey*>(it.Next()))){
    const std::string name=key->GetClassName();
    if(name=="TParameter<double>") {
      TParameter<double>* p =dynamic_cast<TParameter<double>* >(key->ReadObj());
      if(p) v.push_back(*p);
    }
    
  }


}

TTree* define_tree_from_pars(TDirectory* dir, const char* name, const char* title, std::vector<TParameter<double> >& v){
  TDirectory* d = gDirectory;
  if(dir) dir->cd();
  TTree* t = new TTree(name,title);
  for(unsigned int i=0; i<v.size(); i++){
    void* X = static_cast<void*>(const_cast<double*>(&(v[i].GetVal())));
    t->Branch(v[i].GetName(),X,Form("%s/D",v[i].GetName()));
  }
    
  if(d) d->cd();
  return t;
}

void fill_tree_from_pars(TTree* t, std::vector<TParameter<double> >& v){
  if(!t) return;
  for(unsigned int i=0; i<v.size(); i++){
    TBranch* b = t->GetBranch(v[i].GetName());
    if(!b){
      std::cerr<<"fill_tree_from_pars can't find branch : "
	       <<v[i].GetName()<<std::endl;
      continue;
    }
    else {
      void* X = static_cast<void*>(const_cast<double*>(&(v[i].GetVal())));
      b->SetAddress(X);
    }
  }
  t->Fill();
  
}


/////////////// plot collection ////////////////////////////////
//
// a collection of plots with the same attributes
//
////////////////////////////////////////////////////////////////

plot_collection::plot_collection(TH3* h, float scale_factor, const char* projopt, const char* projtit, int* alow, int* ahigh, int na, int* blow, int* bhigh, int nb){
  
  std::vector<int> a1 = make_vector(alow,alow+na);
  std::vector<int> a2 = make_vector(ahigh,ahigh+na);
  std::vector<int> b1 = make_vector(blow,blow+nb);
  std::vector<int> b2 = make_vector(bhigh,bhigh+nb);
  
  v=get_h3_projections(h,projopt,a1,a2,b1,b2,projtit,scale_factor);

}

plot_collection::plot_collection(const char* hname, const char* fname, float scale_factor, const char* projopt, const char* projtit, int* alow, int* ahigh, int na, int* blow, int* bhigh, int nb){
  TH3* h3 = get_h3(hname,fname);

  std::vector<int> a1 = make_vector(alow,alow+na);
  std::vector<int> a2 = make_vector(ahigh,ahigh+na);
  std::vector<int> b1 = make_vector(blow,blow+nb);
  std::vector<int> b2 = make_vector(bhigh,bhigh+nb);
  
  v=get_h3_projections(h3,projopt,a1,a2,b1,b2,projtit,scale_factor);
  
}

plot_collection::plot_collection(TH2* h, float scale_factor, const char* projopt, const char* projtit, int* alow, int* ahigh, int na){
  
  std::vector<int> a1 = make_vector(alow,alow+na);
  std::vector<int> a2 = make_vector(ahigh,ahigh+na);
  
  v=get_projections(h,projopt,a1,a2,projtit,scale_factor);

}

plot_collection::plot_collection(const char* hname, const char* fname, float scale_factor, const char* projopt, const char* projtit, int* alow, int* ahigh, int na){
  TH2* h2 = get_h2(hname,fname);

  std::vector<int> a1 = make_vector(alow,alow+na);
  std::vector<int> a2 = make_vector(ahigh,ahigh+na);
   
  v=get_projections(h2,projopt,a1,a2,projtit,scale_factor);
  
}


//  plot_collection(TH2* h, float scale_factor, const char* projopt, int* alow, int* ahigh, int na);

void plot_collection::set_attline(Int_t c, Int_t s, Int_t w){
  
  for(unsigned int i=0; i<v.size(); i++){
    if(c!=-2) v[i]->SetLineColor(c);
    if(s!=-2) v[i]->SetLineStyle(s);
    if(w!=-2) v[i]->SetLineWidth(w);
  }
}

void plot_collection::set_attfill(Int_t c, Int_t s){
  
  for(unsigned int i=0; i<v.size(); i++){
    if(c!=-2) v[i]->SetFillColor(c);
    if(s!=-2) v[i]->SetFillStyle(s);
  }

}


void plot_collection::set_attmarker(Int_t c, Int_t s, Int_t z){
  
  for(unsigned int i=0; i<v.size(); i++){
    if(c!=-2) v[i]->SetMarkerColor(c);
    if(s!=-2) v[i]->SetMarkerStyle(s);
    if(z!=-2) v[i]->SetMarkerSize(z);
  }
}

void plot_collection::set_range(float* d, float* u, unsigned int n){
  
  for(unsigned int i=0; (i<n && i<v.size()); i++){
    if(d[i]>u[i]) continue;
    v[i]->GetXaxis()->SetRangeUser(d[i],u[i]);
  }

}

void plot_collection::set_rebin(int* r, unsigned int n){

  for(unsigned int i=0; (i<n && i<v.size()); i++){
    if(r[i]<2) continue;
    v[i]->Rebin(r[i]);
  }

}



pc_plotter::pc_plotter(int h, int v, const char* xtit, const char* ytit, const char* cantit, float xmargin, float ymargin)
  : nh(h), nv(v), xt(xtit), yt(ytit), ct(cantit),xm(xmargin),ym(ymargin),leg(0){}

void pc_plotter::add_pc(plot_collection* pc, const char* dopt){
  
  pcs.push_back(pc);
  if(pcs.size()>1) {
    std::string s=dopt; s+=" same";
    opts.push_back(s);
  }
  else opts.push_back(dopt);  
  
}

void pc_plotter::draw_plots(){
  if(pcs.size()<1) return;
  
  plot_many(pcs[0]->v,nh,nv,xt.c_str(),yt.c_str(),opts[0].c_str(),cans,pads,ct.c_str(),xm,ym);
  
  for(unsigned int i=1; i<pcs.size(); i++){
    plot_many(pcs[i]->v, pads, opts[i].c_str());
  }

}

void pc_plotter::set_legend(Double_t x1, Double_t y1, Double_t x2, Double_t y2, const char* header){
  
  leg = new TLegend(x1,y1,x2,y2,header);

  for(unsigned int i=0; i<pcs.size(); i++){
    leg->AddEntry(pcs[i]->v[0],pcs[i]->leg_label.c_str(),pcs[i]->leg_opt.c_str());
  }

}

void pc_plotter::draw_legend(unsigned int ipad){
  if(ipad<pads.size()){
    TVirtualPad* save = gPad;
    pads[ipad]->cd();
    draw_leg(leg);
    save->cd();
  }
}


void make_minos_style(){

  TStyle* minosStyle = new  TStyle("minosStyle", "MINOS Style");
  //set the background color to white
  minosStyle->SetFillColor(10);
  minosStyle->SetFrameFillColor(10);
  minosStyle->SetCanvasColor(10);
  minosStyle->SetPadColor(10);
  minosStyle->SetTitleFillColor(0);
  minosStyle->SetStatColor(10);
  //dont put a colored frame around the plots
  minosStyle->SetFrameBorderMode(0);
  minosStyle->SetCanvasBorderMode(0);
  minosStyle->SetPadBorderMode(0); 
  //use the primary color palette
  minosStyle->SetPalette(1,0);
  //set the default line color for a histogram to be black
  minosStyle->SetHistLineColor(kBlack);
  //set the default line color for a fit function to be red
  minosStyle->SetFuncColor(kRed);
  //make the axis labels black
  minosStyle->SetLabelColor(kBlack,"xyz");
  //set the default title color to be black
  minosStyle->SetTitleColor(kBlack);
  //set the margins
  minosStyle->SetPadBottomMargin(0.2);
  minosStyle->SetPadTopMargin(0.075);
  minosStyle->SetPadLeftMargin(0.15); 
  //set axis label and title text sizes
  minosStyle->SetLabelSize(0.07,"xyz");
  minosStyle->SetTitleSize(0.08,"xyz");
  minosStyle->SetTitleOffset(0.9,"x");
  minosStyle->SetTitleOffset(0.8,"yz");
  minosStyle->SetStatFontSize(0.07);
  minosStyle->SetTextSize(0.08);
  minosStyle->SetTitleBorderSize(0);
  //set line widths
  minosStyle->SetHistLineWidth(2);
  minosStyle->SetFrameLineWidth(2);
  minosStyle->SetFuncWidth(2);
  //set the number of divisions to show
  minosStyle->SetNdivisions(506, "xy");
  //turn off xy grids
  minosStyle->SetPadGridX(0);
  minosStyle->SetPadGridY(0);
  //set the tick mark style
  minosStyle->SetPadTickX(1);
  minosStyle->SetPadTickY(1);
  //show the fit parameters in a box
  minosStyle->SetOptFit(1111);
  //turn off all other stats
  minosStyle->SetOptStat(0000000);
  //marker settings
  minosStyle->SetMarkerStyle(8);
  minosStyle->SetMarkerSize(0.8);
  // Fonts
  minosStyle->SetStatFont(42);
  minosStyle->SetLabelFont(42,"xyz");
  minosStyle->SetTitleFont(42,"xyz");
  minosStyle->SetTextFont(42);

}


TGraphAsymmErrors* setPoissonErrors (TGraphAsymmErrors* g1, TH1D* h1){
  Double_t jimeu[51];
  Double_t jimed[51];
  jimeu[0]=  1.841;  jimed[0]= 0.000 ;
  jimeu[1]=  3.3;    jimed[1]= 0.173 ;
  jimeu[2]=  4.638;  jimed[2]= 0.708 ;
  jimeu[3]=  5.918;  jimed[3]= 1.367 ;
  jimeu[4]=  7.163;  jimed[4]= 2.086 ;
  jimeu[5]=  8.382;  jimed[5]= 2.84  ;
  jimeu[6]=  9.584;  jimed[6]= 3.62  ;
  jimeu[7]=  10.77;  jimed[7]= 4.419 ;
  jimeu[8]=  11.95;  jimed[8]= 5.232 ;
  jimeu[9]=  13.11;  jimed[9]= 6.057 ;
  jimeu[10]= 14.27;  jimed[10]=6.891 ;
  jimeu[11]= 15.42;  jimed[11]=7.734 ;
  jimeu[12]= 16.56;  jimed[12]=8.585 ;
  jimeu[13]= 17.7;   jimed[13]=9.441;
  jimeu[14]= 18.83;  jimed[14]=10.3  ;
  jimeu[15]= 19.96;  jimed[15]=11.17 ;
  jimeu[16]= 21.08;  jimed[16]=12.04 ;
  jimeu[17]= 22.2 ;  jimed[17]=12.92 ;
  jimeu[18]= 23.32;  jimed[18]=13.8  ;
  jimeu[19]= 24.44;  jimed[19]=14.68 ;
  jimeu[20]= 25.55;  jimed[20]=15.57 ;
  jimeu[21]= 26.66;  jimed[21]=16.45 ;
  jimeu[22]= 27.76;  jimed[22]=17.35 ;
  jimeu[23]= 28.87;  jimed[23]=18.24 ;
  jimeu[24]= 29.97;  jimed[24]=19.14 ;
  jimeu[25]= 31.07;  jimed[25]=20.03 ;  
  jimeu[26]= 32.16;  jimed[26]=20.93 ;
  jimeu[27]= 33.26;  jimed[27]=21.84 ;
  jimeu[28]= 34.35;  jimed[28]=22.74 ;
  jimeu[29]= 35.45;  jimed[29]=23.65 ;
  jimeu[30]= 36.54;  jimed[30]=24.55 ;
  jimeu[31]= 37.63;  jimed[31]=25.46;
  jimeu[32]= 38.72;  jimed[32]=26.37;
  jimeu[33]= 39.80;  jimed[33]=27.28;
  jimeu[34]= 40.89;  jimed[34]=28.20;
  jimeu[35]= 41.97;  jimed[35]=29.11;
  jimeu[36]= 43.06;  jimed[36]=30.03;
  jimeu[37]= 44.14;  jimed[37]=30.94;
  jimeu[38]= 45.22;  jimed[38]=31.86;
  jimeu[39]= 46.30;  jimed[39]=32.78;
  jimeu[40]= 47.32;  jimed[40]=33.70;
  jimeu[41]= 48.36;  jimed[41]=34.62;
  jimeu[42]= 49.53;  jimed[42]=35.55;
  jimeu[43]= 50.61;  jimed[43]=36.47;
  jimeu[44]= 51.68;  jimed[44]=37.39;
  jimeu[45]= 52.76;  jimed[45]=38.32;
  jimeu[46]= 53.83;  jimed[46]=39.24;
  jimeu[47]= 54.90;  jimed[47]=40.17;
  jimeu[48]= 55.98;  jimed[48]=41.10;
  jimeu[49]= 57.05;  jimed[49]=42.02;
  jimeu[50]= 58.12;  jimed[50]=42.95;

  Int_t n = h1->GetNbinsX();
  cout << n << " <---test";
  //static const Int_t n;                                                                                                                                               
  cout<<"histogram bins= "<<n<<endl;                                                                                                                                  
  TArrayD yValues(n);
  TArrayD yErrorsLow(n);
  TArrayD yErrorsHigh(n);

  cout<<"defined the arrays"<<endl;                                                                                                                                   
  cout<<"Values\t"<<"Error Low\t"<<"Error High\n"<<endl;                                                                                                              
  for (Int_t i=1; i<=n; i++){
    Double_t x, y;
    g1->GetPoint(i, x, y);
    yValues[i-1]=y;
    cout <<"i: "<< i <<" (i -1): " << i-1<<endl;
    cout <<"x: "<<x<<"y: "<<y<<endl;
    if(y>0&&y<=50){
      yErrorsLow[i-1] = yValues[i-1] - jimed[(Int_t)yValues[i-1]];
      yErrorsHigh[i-1] = jimeu[(Int_t)yValues[i-1]] - yValues[i-1];

      cout<<"yvalue: "<<yValues[i-1]<<" jimed: "<<jimed[(Int_t)yValues[i-1]] <<" jimeu: "<<jimeu[(Int_t)yValues[i-1]]<<endl;
      //cout<<"high: "<<yErrorsHigh[i-1]<<"low: "<<yErrorsLow[i-1]<<endl;
      g1->SetPointEYlow(i, yErrorsLow[i-1]);
      g1->SetPointEYhigh(i, yErrorsHigh[i-1]);
      cout<<"Set Poisson errors for bin "<<i<<endl;                                 }
  }
  return g1;
}
