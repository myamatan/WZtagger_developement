#include<iostream>
#include<fstream>
#include<string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TEventList.h"
#include "TF1.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "sigList.hh"
#include "func.hh"

TString  boson = "W";
TString  chs = "1Lep";
int nbin = 45;//30
double mShift = 0;
const int npback = 5;
std::vector<float> ptRanges = {500,700};
TString memo = "500_700";
double x1=50;
double x2=200;
bool fitWZmass  = false;   // Switch on to float W mass. Defaulted to off to ensure closure linearity.
bool fitWZwidth = false;   // Switch on to float main Gaussian width. Default to off for above reason.
double secondGaussianRelativeWidth = 2.;   // Starting value only
double secondGaussianRate = 0.20;          // Starting value only
bool fitWithSecondGaussian = true; // Ratio of single to double Gaussian difference from 1 taken as syst
double fitWwidth;
double dfitWwidth;

struct fitresult {
  double sc;      // Signal scale parameter
  double dsc;     // Error on signal scale parameter
  double nev;     // Number of W+Z events from fit on data
  double dnev;    // Stat error on above
  double nevexp;  // Expected number from fit on MC
  double mu;      // Measured SF = nev/nevexp
  double dmu;     // Stat error on above
  double mwshift; // Peak position shift
  double dmwshift;// Error on mwshift
};

// Fit function
Double_t bkdfunc1D(Double_t *x, Double_t *par) {
  Double_t mass=x[0]-100;
  double scale = par[0];
  double toexp = 0;
  double massp = mass;
  for (int i=1;i<npback;++i) {
    toexp += massp*par[i];
    massp *= mass;
  }
  double y = scale*exp(toexp);
  return y;
}

Double_t sigfuncWZ1D(Double_t *x, Double_t *par,double massWZ) {
  // Fit a Double Gassian with tails centred  on  massWZ
  // The second Gaussian has a width which is a multiple of the main Gaussian width.
  Double_t mass=x[0];
  Double_t scale = par[0];
  Double_t offset = par[1];
  Double_t wid = par[2];
  Double_t tailWid = par[3]*par[2];
  Double_t tailScale = par[4];  

  if (tailWid == 0 || wid == 0) {
    return 0;
  }
  double dm=(mass-massWZ-offset)/wid;
  double dmt=(mass-massWZ-offset)/tailWid;
  
  Double_t z =  exp(-0.5*dm*dm) + tailScale*exp(-0.5*dmt*dmt);
  z *= scale;
  return z;
}
Double_t sigfuncW1D(Double_t *x, Double_t *par) {
  return sigfuncWZ1D(x,par,80.3);
}
Double_t sigfuncZ1D(Double_t *x, Double_t *par) {
  return sigfuncWZ1D(x,par,91.18);
}
Double_t totWjets(Double_t *x, Double_t *par) { // W+jets don't need qcd bkg anymore since we have truth matching
  double result = sigfuncWZ1D(x,par,80.3);
  result += bkdfunc1D(x,par+5);
  return result;
}
Double_t totZjets(Double_t *x, Double_t *par) {// Z+jets modell as double Gaussina, added bkdg func not necessary
  double result = sigfuncWZ1D(x,par,91.18);
  result += bkdfunc1D(x,par+5);
  return result;
}
Double_t datfunc1D(Double_t *x, Double_t *par) {
  // 0-2 control signal fit, 10-19 is w pars, 20-29 z pars and 30-50 bkd.
  double xscale[1];
  double parsW[50];
  double parsZ[50];
  double signalrate = par[0];
  double zfrac = par[1];
  xscale[0] = x[0];

  parsW[0] = signalrate*(1-zfrac);
  for (int i=1;i<5;++i) parsW[i] = par[10+i];  // W parameters

  parsZ[0] = signalrate*zfrac;
  for (int i=1;i<5;++i) parsZ[i] = par[20+i]; // Z parameters

  // Now make the Z peak offset be additive to the W one
  parsZ[1] += parsW[1];

  double bkd = bkdfunc1D(x,par+30);
  double sigW = sigfuncW1D(xscale,parsW);
  double sigZ = sigfuncZ1D(xscale,parsZ);
  double y = bkd+sigW+sigZ;
  return y;
}

double   getFittedWZ(TH1F * hist,TF1 * fn,int fitType,double x1,double x2) {
  // Fittype: 0 for W, 1 for Z 2 for data
  std::cout<<"getFittedWZ: 1D starting"<<std::endl;

  double total = 0;
  int nps=0;
  TF1 * f_nobk=0;

  switch (fitType) {
  case 0: // w
    nps = 5;
    f_nobk = new  TF1("f_nobk",totWjets,x1,x2,10);
    break;
  case 1: // z
    nps = 5;
    f_nobk = new  TF1("f_nobk",totZjets,x1,x2,10);
    break;
  case 2:
    nps = 30;
    f_nobk = new  TF1("f_nobk",datfunc1D,x1,x2,nps+npback);
    break;
  default:
    std::cout<<"getFittedWZ called for unknown fitType-"<<fitType<<std::endl;
    return 0;
  }    
  //npback==5
  // Copy the signal parameters
  for (int i=0;i<nps;++i) f_nobk->FixParameter(i,fn->GetParameter(i));
  // Zero the bkd component.
  for (int i=nps;i<nps+npback;++i) {
    //std::cout<<"getFittedWZ Zeroing parameter "<<i<<std::endl;
    f_nobk->FixParameter(i,0);
  }

  int nx = hist->GetXaxis()->GetNbins();
  
  for (int ix=0;ix<nx;++ix) {
    double x = hist->GetXaxis()->GetBinCenter(ix+1);
    if (x < x1 || x > x2) continue;
    double rate =  f_nobk->Eval(x);
    total += rate;
  }
  delete f_nobk;
  return total;  
}


void sigbkg_modeling(TString cutDef="veryloose", TString WP="MaxVVSemi_v5", int Nbin=30, double MShift=30){

	SetAtlasStyle();
	TF1 *fD2,*fML,*fMU;

	nbin=Nbin;
	mShift=MShift;

	TString wpName = WP;
	//TString wpName = "MaxVVJJ_3D";
	//TString wpName = "MaxVVJJ";

	TFile *f[3];
	TTree *tr[3];
	TH1F *h[3];
	h[0] = new TH1F("h" ,"h" , nbin,x1,x2);
	h[1] = new TH1F("hW","hW", nbin,x1,x2);
	h[2] = new TH1F("hZ","hZ", nbin,x1,x2);

	// W/Z tagger cut function
	std::string wpfile = std::string("../output/WPs/SmoothedContained"+boson+wps[wpName]);
	std::vector<TString> cutWPs = retrieveWPs( wpfile );
	fD2 = new TF1("D2"+boson, cutWPs.at(2), 200, 2500);
	fMU = new TF1("MU"+boson, cutWPs.at(1), 200, 2500);
	fML = new TF1("ML"+boson, cutWPs.at(0), 200, 2500);
	float D2cut=0.;

	// Pre-cut definition
	TCut precut;
	if( cutDef=="tight") precut = "fatjet_D2>0 && fatjet_m >=50 &&" + preCutfor1Lep;
	else if( cutDef=="veryloose") precut = "fatjet_D2>0 && fatjet_m >=50 &&" + preCutfor1Lep_veryloose;
	else precut = "fatjet_D2>0 && fatjet_m >=50 &&" + preCutfor1Lep_loose;
	
	// Background function
	TF1 *bkdFcn = new TF1("bkdFcn1D",bkdfunc1D,x1,x2,npback);
  	TString pname="sl";
  	bkdFcn->SetParName(0, pname+" scale");
  	for (int i=1;i<npback;++i) {
  	  TString parname = pname+" p";
  	  parname += i;
  	  parname += "   ";
  	  bkdFcn->SetParName(i,parname);
  	}
	bkdFcn->SetParameter(0,1);
  	if (npback > 1) bkdFcn->SetParameter(1, -0.01);
  	if (npback > 2) bkdFcn->SetParameter(2, 0.0001);
  	if (npback > 3) bkdFcn->SetParameter(3, 1.E-6);
  	if (npback > 4) bkdFcn->SetParameter(4, 0.);

	// Background only Fit
	nameOption = "Wjets";
	TString fname = "../input/merged_new/"+nameOption+"/"+nameOption+"." + chs + ".root";
	std::cout << "++++++ " << fname << " ++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	f[0] = new TFile(fname);
	tr[0] = (TTree*)f[0]->Get("Nominal");
	SetBranch(tr[0], chs);
			
	tr[0]->Draw(">>elist",precut);
	TEventList *elist = (TEventList*)gROOT->FindObject("elist");
	long nEvent = elist->GetN();
	for(long l=0; l<nEvent; l++){
		tr[0]->GetEntry(elist->GetEntry(l));
		D2cut = fD2->Eval(pt);
		if( pt > ptRanges.at(0) && pt <= ptRanges.at(1) ){
		    if( d2<D2cut ) h[0]->Fill( m, weight );
		}
		if( l%(int(nEvent/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEvent, float(float(l)/nEvent)*100) << std::endl;
	}
	std::cout << h[0]->Integral() << std::endl;

	h[0]->Fit(bkdFcn,"R0","",x1,x2);
	for (int i=0;i<5;++i) bkdFcn->ReleaseParameter(i);
	h[0]->Fit(bkdFcn,"R0E","",x1,x2);//for error

	std::cout << "Reduced chi2 : " << bkdFcn->GetChisquare()/bkdFcn->GetNDF() << std::endl;
	//std::cout << "GetChisquare : " << f1->GetChisquare() << std::endl;
	//std::cout << "GetNDF : " << f1->GetNDF() << std::endl;
	//std::cout << "GetNumberFitPoints : " << f1->GetNumberFitPoints() << std::endl;
	//std::cout << "GetNumberFreeParameters : " << f1->GetNumberFreeParameters() << std::endl;

	// Signal only Fit
	TF1 *f_sigWZ[2];
	double fittedWZ[2];
	double sw=0, sz=0;
	for(int i=0; i<2; i++){
		std::cout << i << std::endl;
		if( i==0 )nameOption = "signalW";
		else if( i==1 )nameOption = "signalZ";
		fname = "../input/merged/"+nameOption+"/"+nameOption+"." + chs + ".root";
		std::cout << "++++++ " << fname << " ++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		f[i+1]  = new TFile(fname);
		tr[i+1] = (TTree*)f[i+1]->Get("Nominal");
		SetBranch(tr[i+1], chs);

		tr[i+1]->Draw(">>elist",precut);
		elist = (TEventList*)gROOT->FindObject("elist");
		nEvent = elist->GetN();
		for(long l=0; l<nEvent; l++){
			tr[i+1]->GetEntry(elist->GetEntry(l));
			D2cut = fD2->Eval(pt);
			if( pt > ptRanges.at(0) && pt <= ptRanges.at(1) ){
			    if( d2<D2cut ) h[i+1]->Fill( m+mShift, weight );
			}
			if( l%(int(nEvent/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEvent, float(float(l)/nEvent)*100) << std::endl;
		}
		std::cout << "hWZ[i]->Integral() : " << h[i+1]->Integral() << std::endl;

		TString part = "W";
	        if (i==1) part="Z";
	    
	        if (i==0) f_sigWZ[i] = new TF1("f_sigW",totWjets,x1,x2,10);
	        else if( i==1 )f_sigWZ[i] = new TF1("f_sigZ",totZjets,x1,x2,10);
	    
	        f_sigWZ[i]->SetParName(0,"scale     ");
	        f_sigWZ[i]->SetParName(1,"offset    ");
	        f_sigWZ[i]->SetParName(2,"wid       ");
	        f_sigWZ[i]->SetParName(3,"tailRelWid");
	        f_sigWZ[i]->SetParName(4,"tailscale ");
	        f_sigWZ[i]->SetParName(5,"bkd scale ");

		for (int j=1;j<npback;j++) {
	  	      TString parname = "bkd p";
	  	      parname += j;
	  	      parname += "    ";
	  	      f_sigWZ[i]->SetParName(5+j,parname);
	  	      f_sigWZ[i]->FixParameter(5+j,bkdFcn->GetParameter(j));
	  	}
	  	
	  	f_sigWZ[i]->SetParameter(0,h[i+1]->GetMaximum());
	  	f_sigWZ[i]->SetParameter(1,4.);
	  	f_sigWZ[i]->SetParameter(2,8.);
	  	f_sigWZ[i]->SetParameter(3,secondGaussianRelativeWidth);
	  	f_sigWZ[i]->FixParameter(4,secondGaussianRate);    
	  	f_sigWZ[i]->FixParameter(5,f_sigWZ[i]->GetParameter(0)/10.);
	  	f_sigWZ[i]->SetParameter(5,0.); //background scale in signal modelling
	  	
	  	if (!fitWithSecondGaussian) {
	  	  std::cout<<"Horeay, I switched off the tail"<<std::endl;
	  	  f_sigWZ[i]->FixParameter(3,1.);
	  	  f_sigWZ[i]->FixParameter(4,0.00);
	  	}

		std::cout<<"About to fit WZ "<<i<<std::endl;
	        h[i+1]->Fit(f_sigWZ[i],"lR0","",x1,x2);
	    
	        if (fitWithSecondGaussian) {
	          f_sigWZ[i]->ReleaseParameter(3);
	          f_sigWZ[i]->ReleaseParameter(4);
	        }
	    
	        f_sigWZ[i]->ReleaseParameter(5); 
	        // Important! For VV analysis QCD shape under W/Z looks same as in general QCD - only need to float scale
	        // This is NOT TRUE in general - may need to fit all. But its denagerous...too flexibile.
	        // It would be so much better to have truth tagging in the W/Z signal.
	        //f_sigWZ[i]->ReleaseParameter(6);
	        //f_sigWZ[i]->ReleaseParameter(7);
	        //f_sigWZ[i]->ReleaseParameter(8);
	        //f_sigWZ[i]->ReleaseParameter(9);
	        h[i+1]->Fit(f_sigWZ[i],"R0","",x1,x2);
	    
	        switch (i) {
	        case 0: sw = f_sigWZ[0]->GetParameter(0); break;
	        case 1: sz = f_sigWZ[1]->GetParameter(0); break;
	        }
	        fittedWZ[i] = getFittedWZ(h[i+1],f_sigWZ[i],i,x1,x2); //(MC_hists, function, 0=W, 1=Z, xmin, xmax)
	}


	// Data fit (closure here --> signal + background fit on bakground )
	double calczscale = 0.5;
        double totz = h[1]->GetSumOfWeights();
        double totw = h[2]->GetSumOfWeights();
        std::cout<<" "<<std::endl;
        std::cout<<"Find "<<totz<<" z and "<<totw<<" w giving frac="<<(totz/(totw+totz))<<" of "<<(totz+totw)<<std::endl;
        double wfint =  f_sigWZ[0]->Integral(x1,x2);
        double zfint =  f_sigWZ[1]->Integral(x1,x2);
        std::cout<<"Integral W="<<wfint<<" integral z="<<zfint<<std::endl;
        std::cout<<"gives fraction of "<<(zfint/(wfint+zfint))<<std::endl;
        // find calczscale: height of the Z peak relative to the W
        // Will be used in the data function to link them together
        calczscale = sz/(sw+sz);
        std::cout<<"sw="<<sw<<" sz="<<sz<<" giving z scale="<<calczscale<<std::endl;
        TF1 *f_dat = new TF1("f_dat",datfunc1D,x1,x2,30+npback);
        std::cout<<"Next set the data pars"<<std::endl;
        f_dat->SetParName(0,"SignalRate ");
        f_dat->SetParName(1,"Z rel scale");
        std::cout<<"About to define pars 3-9"<<std::endl;
        for (int i=2;i<10;++i) {
          f_dat->SetParName(i,"dummy      ");
          f_dat->FixParameter(i,0);
        }
        std::cout<<"Just defined pars 2-9 as dummy"<<std::endl;
        // Next copy in the W and Z parameters
        for (int iwz=0;iwz<2;++iwz) {
          f_dat->SetParName(10*iwz+10,"dummy       ");
          f_dat->FixParameter(10*iwz+10,0);
          for (int i=1;i<5;++i) {
            f_dat->SetParName(iwz*10+10+i,f_sigWZ[iwz]->GetParName(i));
            //f_dat->FixParameter(iwz*10+10+i,0.);
            f_dat->FixParameter(iwz*10+10+i,f_sigWZ[iwz]->GetParameter(i));
          }
          for (int i=5;i<10;++i) {
            f_dat->SetParName(iwz*10+10+i,"dummy       ");
            f_dat->FixParameter(iwz*10+10+i,0);
          }
        }
        // Make the Z offset additive to the W one. 
        f_dat->FixParameter(21,f_dat->GetParameter(21)-f_dat->GetParameter(11));
        // Now copy the QCD parameers in
        for (int i=0;i<npback;++i) {
          f_dat->SetParName(30+i,bkdFcn->GetParName(i));
          f_dat->FixParameter(30+i,bkdFcn->GetParameter(i));
        }
      
        f_dat->ReleaseParameter(30);// Release QCD rate               
      
        // Make a first guess at signal parameters:
        f_dat->SetParameter(0,10000);
        f_dat->FixParameter(1,calczscale);
        f_dat->FixParameter(2,1.);
        h[0]->Fit(f_dat,"qR0"," ",x1,x2); //First reduced fit
        
        // Next release the QCD shape parameters
        for (int i=31;i<30+npback;++i) {
          f_dat->ReleaseParameter(i);                    
        }                                               
      
        // Fit W mass if aksed to. Z mass floats off it.
        double fittedW = f_dat->GetParameter(11); // mW
        if (fitWZmass) {
          f_dat->ReleaseParameter(11); // mW
          f_dat->SetParLimits(11,fittedW - 5,fittedW+5); //Allow 5 GeV range
          f_dat->ReleaseParameter(21); // mZ
        }
        double fittedWwidth = f_dat->GetParameter(12); // sigma_W
      
        if (fitWZwidth) {
          f_dat->ReleaseParameter(12); // GammaW
          f_dat->SetParLimits(12,0.9*fittedWwidth,1.1*fittedWwidth); // Allow 30% variation
        }
      
        //  f_dat->ReleaseParameter(1); // z scale
        // Do a full fit. with "R" to respect range.  
        h[0]->Fit(f_dat,"qR0"," ",x1,x2);
        // Repeat, calling MINOS for errors too.
        h[0]->Fit(f_dat,"qMR0"," ",x1,x2);
        // Fit AGAIN, returning fit results. 
        // This  triplicated  fit is partly because it is more lmore likely to converge.
        TFitResultPtr fitres = h[0]->Fit(f_dat,"ER0S"," ",x1,x2);
        double shiftW = 0;
        double dshiftW = 0;
      
        if (fitWZmass) {
          shiftW = f_dat->GetParameter(11) - fittedW;
          dshiftW = f_dat->GetParError(11);
          std::cout<<"Shift in mW is "<<shiftW<<"+-"<<dshiftW<<" for "<<std::endl;
        }
        if (fitWZwidth) {
          fitWwidth = f_dat->GetParameter(12);
          dfitWwidth = f_dat->GetParError(12);
          std::cout<<"Fitted w width is "<<fitWwidth<<"+-"<<dfitWwidth<<" c/f MC "<<fittedWwidth<<" for "<<std::endl;
        }
        //  std::cout<<"Passed w mass shift"<<std::endl;
      
        double totalWZpred = fittedWZ[0]+fittedWZ[1]; //Number of fitted W/Z from MC
        double dWpred = f_sigWZ[0]->GetParError(0)/f_sigWZ[0]->GetParameter(0)*fittedWZ[0];
        double dZpred = f_sigWZ[1]->GetParError(0)/f_sigWZ[1]->GetParameter(0)*fittedWZ[1];
        double dtotalPred = pow((dWpred*dWpred + dZpred*dZpred), .5);
        
      
        double totalWZfit = getFittedWZ(h[0],f_dat,2,x1,x2); //Number of fitted W/Z from Data
        std::cout<<"Passed getFitted, data totalZW="<<totalWZfit<<" W_jets samp="<<fittedWZ[0]
      	   <<" Z_jets samp="<<fittedWZ[1]<<std::endl;
        double dtotalWZfit = f_dat->GetParError(0)/f_dat->GetParameter(0)*totalWZfit;
        std::cout<<"Passed dtotalwz"<<std::endl;
        std::cout << "dWpred = " << dWpred << std::endl;
        std::cout << "dZpred = " << dZpred << std::endl;
      
        TString showName = "Data";
        std::cout<<"About to call showFit, totalWZpred="<<totalWZpred<<std::endl;

	fitresult theresult;
  	theresult.sc = f_dat->GetParameter(0); //signal scale
  	theresult.dsc = f_dat->GetParError(0); //signal scale error
  	theresult.nev = totalWZfit;
	theresult.dnev = totalWZfit*f_dat->GetParError(0)/f_dat->GetParameter(0);
  	theresult.nevexp = totalWZpred;
  	theresult.mwshift = shiftW;
  	theresult.dmwshift = dshiftW;
  	theresult.mu = theresult.nev / totalWZpred;
  	theresult.dmu = theresult.mu*(theresult.dsc/theresult.sc);



  	std::cout << "WZfitted = " << totalWZfit << " +/- " << theresult.dnev << "<<<<<<<<<<<<" <<std::endl;
  	//std::cout << "WZpredic = " << totalWZpred << "<<<<<<<<<<<<" <<std::endl;
  	std::cout << "mu       = " << theresult.mu << " +/- " << theresult.dmu << "<<<<<<<<<<<<" <<std::endl;
  	std::cout << ">>>>>>>>>>>>>>>>>>>>>>gamma_w = " << f_dat->GetParameter(12) << " <<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

	// Plot
	TCanvas    *cv[4];
	TLegend   *leg[4];
	TLatex    *lat[4];
	TLatex *latsub[4];
	TLatex *latStat[4];
	
	for(UInt_t i=0; i<3; i++){
		cv[i] = new TCanvas(Form("c1%d",i),"Canvas",200,10,500,480);
		leg[i]= new TLegend(0.66,0.32,0.88,0.45);
		leg[i]->SetBorderSize(0);


        	float k_Nbin  = float(h[i]->GetNbinsX()+1);
        	float k_range = float(h[i]->GetXaxis()->GetXmax() - h[i]->GetXaxis()->GetXmin());
        	float k_xmin  = float(h[i]->GetXaxis()->GetXmin()) + k_range/k_Nbin*float(h[i]->FindFirstBinAbove());
        	float k_xmax  = float(h[i]->GetXaxis()->GetXmin()) + k_range/k_Nbin*float(h[i]->FindLastBinAbove());
        	float XrangeMax = k_xmax;
        	float XrangeMin = k_xmin;
        	float k_ymax = h[i]->GetBinContent( h[i]->GetMaximumBin() );
        	float YrangeMax = k_ymax * 1.1;//40.1
        	float YrangeMin = k_ymax * 0;//1000.

		float chi2_ndf;
		if(i==0) chi2_ndf=bkdFcn->GetChisquare()/bkdFcn->GetNDF();
		else chi2_ndf=f_sigWZ[i-1]->GetChisquare()/f_sigWZ[i-1]->GetNDF();
        	lat[i] = new TLatex(XrangeMin+(XrangeMax-XrangeMin)*0.35,
        	                YrangeMin+(YrangeMax-YrangeMin)*0.55,
        	                "#splitline{ #it{#bf{ATLAS}} #it{#scale[0.8]{Simulation Internal}} }{ #scale[0.8]{#sqrt{s} = 13 TeV, "+boson+" tagger} }");
        	latsub[i] = new TLatex(XrangeMin+(XrangeMax-XrangeMin)*0.42,
        	                YrangeMin+(YrangeMax-YrangeMin)*0.40,
        	                Form("#scale[0.8]{#chi^{2}/ndf=%0.3f}",chi2_ndf));

		cv[i]->cd();
		h[i]->Draw();
		h[i]->SetStats(0);
		h[i]->GetYaxis()->SetRangeUser(YrangeMin, YrangeMax);
		h[i]->GetYaxis()->SetTitle("Events");
		h[i]->GetXaxis()->SetTitle("m_{J} [GeV]");
		h[i]->SetMarkerSize(1);
		h[i]->SetMarkerStyle(8);
		h[i]->SetLineWidth(2);

		if( i==0 ){
			leg[i]->AddEntry(h[i], "W(lv)+jets", "lep");
			leg[i]->AddEntry(bkdFcn, "Fit", "l");
			bkdFcn->SetLineColor(kGreen-3);
			bkdFcn->Draw("same");
			//leg->AddEntry(h, Form("%d < p_{T}(J) < %d GeV",int(ptRanges.at(i).at(0)),int(ptRanges.at(i).at(1))), "lep");
		}else if(i==1){
			leg[i]->AddEntry(h[i], "W signal", "lep");
			leg[i]->AddEntry(f_sigWZ[i-1], "Fit", "l");
			f_sigWZ[i-1]->SetLineColor(kGreen-3);
			f_sigWZ[i-1]->Draw("same");
			//leg->AddEntry(h, Form("%d < p_{T}(J) < %d GeV",int(ptRanges.at(i).at(0)),int(ptRanges.at(i).at(1))), "lep");
		}else if(i==2){
			leg[i]->AddEntry(h[i], "Z signal", "lep");
			leg[i]->AddEntry(f_sigWZ[i-1], "Fit", "l");
			f_sigWZ[i-1]->SetLineColor(kGreen-3);
			f_sigWZ[i-1]->Draw("same");
			//leg->AddEntry(h, Form("%d < p_{T}(J) < %d GeV",int(ptRanges.at(i).at(0)),int(ptRanges.at(i).at(1))), "lep");
		}

		leg[i]->Draw();
		lat[i]->Draw();
		latsub[i]->Draw();
		if( i==0 )nameOption = "Wjets";
		else if( i==1 )nameOption = "signalW";
		else if( i==2 )nameOption = "signalZ";

		TString str_mShift = Form("mShift%d",int(mShift));	
		cv[i]->Print("../output/sigbkg_modeling." + cutDef + "." + boson + "." + nameOption + "." + wpName + "." + Form("nb%d",nbin) + "." + memo + "." + str_mShift + ".eps");
	}

	// Closure test plot
	int c_index=3;
	cv[c_index] = new TCanvas("closure","Canvas",200,10,500,480);
	leg[c_index]= new TLegend(0.66,0.32,0.88,0.45);
	leg[c_index]->SetBorderSize(0);

	float k_Nbin  = float(h[0]->GetNbinsX()+1);
	float k_range = float(h[0]->GetXaxis()->GetXmax() - h[0]->GetXaxis()->GetXmin());
	float k_xmin  = float(h[0]->GetXaxis()->GetXmin()) + k_range/k_Nbin*float(h[0]->FindFirstBinAbove());
	float k_xmax  = float(h[0]->GetXaxis()->GetXmin()) + k_range/k_Nbin*float(h[0]->FindLastBinAbove());
	float XrangeMax = k_xmax;
	float XrangeMin = k_xmin;
	float k_ymax = h[0]->GetBinContent( h[0]->GetMaximumBin() );
	float YrangeMax = k_ymax * 1.1;//40.1
	float YrangeMin = k_ymax * 0;//1000.
	float chi2_ndf=f_dat->GetChisquare()/f_dat->GetNDF();
	lat[c_index] = new TLatex(XrangeMin+(XrangeMax-XrangeMin)*0.35,
	                YrangeMin+(YrangeMax-YrangeMin)*0.55,
	                "#splitline{ #it{#bf{ATLAS}} #it{#scale[0.8]{Simulation Internal}} }{ #scale[0.8]{#sqrt{s} = 13 TeV, "+boson+" tagger} }");
	latsub[c_index] = new TLatex(XrangeMin+(XrangeMax-XrangeMin)*0.42,
	                YrangeMin+(YrangeMax-YrangeMin)*0.40,
	                Form("#scale[0.8]{#chi^{2}/ndf=%0.3f}",chi2_ndf));
        latStat[c_index] = new TLatex(XrangeMin+(XrangeMax-XrangeMin)*0.08,
                           YrangeMin+(YrangeMax-YrangeMin)*0.14,
                           //Form("#scale[0.8]{ #splitline{ Fitted WZ: %3.2f #pm %3.2f}{ #mu=%1.3f #pm %1.3f} }",totalWZfit,theresult.dnev,theresult.mu,theresult.dmu));
                           Form("#scale[0.8]{ #splitline{ Fitted WZ: %3.2f #pm %3.2f}{} }",totalWZfit,theresult.dnev));

	cv[c_index]->cd();
	h[0]->Draw();
	h[0]->SetStats(0);
	h[0]->GetYaxis()->SetRangeUser(YrangeMin, YrangeMax);
	h[0]->GetYaxis()->SetTitle("Events");
	h[0]->GetXaxis()->SetTitle("m_{J} [GeV]");
	h[0]->SetMarkerSize(1);
	h[0]->SetMarkerStyle(8);
	h[0]->SetLineWidth(2);

	leg[c_index]->AddEntry(h[0], "W(lv)+jets", "lep");
	leg[c_index]->AddEntry(f_dat, "Fit", "l");
	leg[c_index]->AddEntry(bkdFcn, "Fit bkd", "l");
	f_dat->SetLineColor(kGreen-3);
	f_dat->Draw("same");
	bkdFcn->SetLineColor(kGreen-3);
	bkdFcn->SetLineStyle(2);
	bkdFcn->Draw("same");

	leg[c_index]->Draw();
	lat[c_index]->Draw();
	latsub[c_index]->Draw();
	latStat[c_index]->Draw();
	nameOption="Wjets";
	TString str_mShift = Form("mShift%d",int(mShift));
	cv[c_index]->Print("../output/sigbkg_modeling." + cutDef + "." + boson + "." + nameOption + "." + wpName + "." + Form("nb%d",nbin) + "." + memo + "." + str_mShift + ".closure.eps");

}
