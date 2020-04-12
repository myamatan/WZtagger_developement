#include<iostream>
#include<fstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TCut.h"
#include "TEventList.h"
#include "TF1.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

TString filenameSig = "../input/merged/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged/allMC/allMC."+ch+".root";

void fsigni(){

	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");

	// Cut function
	TF1 *fD2Max;
	TF1 *fMLMax;
	TF1 *fMUMax;
	TF1 *fD2Fix;
	TF1 *fMLFix;
	TF1 *fMUFix;
	TF1 *fD2;
	TF1 *fML;
	TF1 *fMU;
	float d2cut = 0;
	float mcutu = 0;
	float mcutl = 0;
	// MaxSignificance
	std::string wpfile = std::string("../output/WPs/SmoothedContained"+tagger+"Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_MC16a_20171201.dat");
	std::vector<TString> cutWPs = retrieveWPs( wpfile );
	fD2Max = new TF1("D2Max", cutWPs.at(2), 200, 2500);
	fMUMax = new TF1("MUMax", cutWPs.at(1), 200, 2500);
	fMLMax = new TF1("MLMax", cutWPs.at(0), 200, 2500);
	cutWPs.clear();
	// FixedSignalEfficiency50
	wpfile = std::string("../output/WPs/SmoothedContained"+tagger+"Tagger_AntiKt10TrackCaloClusterTrimmed_FixedSignalEfficiency50_MC16a_20171201.dat");
	cutWPs = retrieveWPs( wpfile );
	fD2Fix = new TF1("D2Fix", cutWPs.at(2), 200, 2500);
	fMUFix = new TF1("MUFix", cutWPs.at(1), 200, 2500);
	fMLFix = new TF1("MLFix", cutWPs.at(0), 200, 2500);
	cutWPs.clear();

	// My result
	wpfile = std::string("../output/WPs/SmoothedContained"+tagger+"Tagger_AntiKt10TrackCaloClusterTrimmed_"+ch+"_"+sigType+"_"+nameOption+".txt");
	//wpfile = std::string("../output/WPs/SmoothedContained"+tagger+"Tagger_AntiKt10TrackCaloClusterTrimmed_Max_VV_combinedFit.txt");
	cutWPs = retrieveWPs( wpfile );
	fD2 = new TF1("D2", cutWPs.at(2), 200, 2500);
	fMU = new TF1("MU", cutWPs.at(1), 200, 2500);
	fML = new TF1("ML", cutWPs.at(0), 200, 2500);


	// Histograms for signal & background
	TH1F *hb= new TH1F("hb", "hb", nbins, rebin21);
	TH1F *hbMax= new TH1F("hbMax", "hbMax", nbins, rebin21);
	TH1F *hbFix= new TH1F("hbFix", "hbFix", nbins, rebin21);
	TH1F *hs[NmassPoints];
	TH1F *hsMax[NmassPoints];
	TH1F *hsFix[NmassPoints];
	for(UInt_t l=0; l<NmassPoints; l++){
		TString nameID = TString(std::to_string(l));
		hs[l]= new TH1F(nameID+"hs", nameID+"hs", nbins, rebin21);
		hsMax[l]= new TH1F(nameID+"hsMax", nameID+"hsMax", nbins, rebin21);
		hsFix[l]= new TH1F(nameID+"hsFix", nameID+"hsFix", nbins, rebin21);
	}

	// Background
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
	TCut precut = returnTCut(preSele);
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN(); //long nEventB = trB->GetEntries();
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l)); //trB->GetEntry(l);
		// Max D2 &&  Mass window
		d2cut = pt > 2500 ? fD2Max->Eval(2500) : fD2Max->Eval(pt);
		mcutu = pt > 2500 ? fMUMax->Eval(2500) : fMUMax->Eval(pt);
		mcutl = pt > 2500 ? fMLMax->Eval(2500) : fMLMax->Eval(pt);
		if( d2 < d2cut && m < mcutu && m > mcutl ){
			hbMax->Fill(vvm, weight*lumi/40.);
		}
		// Fix D2 &&  Mass window
		d2cut = pt > 2500 ? fD2Fix->Eval(2500) : fD2Fix->Eval(pt);
		mcutu = pt > 2500 ? fMUFix->Eval(2500) : fMUFix->Eval(pt);
		mcutl = pt > 2500 ? fMLFix->Eval(2500) : fMLFix->Eval(pt);
		if( d2 < d2cut && m < mcutu && m > mcutl ){
			hbFix->Fill(vvm, weight*lumi/40.);
		}
		// My D2 &&  Mass window
		d2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
		mcutu = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
		mcutl = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
		if( d2 < d2cut && m < mcutu && m > mcutl ){
			hb->Fill(vvm, weight*lumi/40.);
		}
		if( l%(int(nEventB/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventB, float(float(l)/nEventB)*100) << std::endl;
	}

	// Signal
	std::cout << "Signal filling ....." << std::endl;
	std::map<int, int> sigMap;
	UInt_t sigCount = 0;
	UInt_t sigArr[10];
	SetBranch(trS, ch);
	precut = returnTCut(preSele);
        trS->Draw(">>elist",precut);
        TEventList *elistS = (TEventList*)gROOT->FindObject("elist");
	long nEventS = elistS->GetN(); //long nEventS = trS->GetEntries();
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(mc) != sigWeight.end() ){
			weightS = weight * lumi * sigWeight[mc];
			//weightS = lumi * sigWeight[mc];
			//weightS = weight;
			//weightS = 1.0;
			if( sigMap.find(mc) == sigMap.end() ){
				sigMap[mc] = sigCount; 
				sigArr[sigCount] = mc;
				sigCount++;
			}
		}
		// Max D2 &&  Mass window
		d2cut = pt > 2500 ? fD2Max->Eval(2500) : fD2Max->Eval(pt);
		mcutu = pt > 2500 ? fMUMax->Eval(2500) : fMUMax->Eval(pt);
		mcutl = pt > 2500 ? fMLMax->Eval(2500) : fMLMax->Eval(pt);
		if( d2 < d2cut && m < mcutu && m > mcutl ){
			hsMax[sigMap[mc]]->Fill(vvm, weightS);
		}
		// Fix D2 &&  Mass window
		d2cut = pt > 2500 ? fD2Fix->Eval(2500) : fD2Fix->Eval(pt);
		mcutu = pt > 2500 ? fMUFix->Eval(2500) : fMUFix->Eval(pt);
		mcutl = pt > 2500 ? fMLFix->Eval(2500) : fMLFix->Eval(pt);
		if( d2 < d2cut && m < mcutu && m > mcutl ){
			hsFix[sigMap[mc]]->Fill(vvm, weightS);
		}
		// VV semi D2 & Mass windo
		d2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
		mcutu = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
		mcutl = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
		if( d2 < d2cut && m < mcutu && m > mcutl ){
			hs[sigMap[mc]]->Fill(vvm, weightS);
		}
		if( l%(int(nEventS/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}

	// Significance calculation
	std::cout << "Optimized cut searching  ....." << std::endl;
	float signi[NmassPoints]={0.0};
	float signiMax[NmassPoints]={0.0};
	float signiFix[NmassPoints]={0.0};
	float error[NmassPoints]={0.0};
	float errorMax[NmassPoints]={0.0};
	float errorFix[NmassPoints]={0.0};
	// D2 & Mass window
	for(UInt_t m=0; m<sigCount; m++){
		signi[m] = significance(hs[m], hs[m], hb, "mass");
		signiMax[m] = significance(hsMax[m], hsMax[m], hbMax, "mass");
		signiFix[m] = significance(hsFix[m], hsFix[m], hbFix, "mass");
		error[m] = significance(hs[m], hs[m], hb, "mass", "error");
		errorMax[m] = significance(hsMax[m], hsMax[m], hbMax, "mass", "error");
		errorFix[m] = significance(hsFix[m], hsFix[m], hbFix, "mass", "error");
	}

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	/*
	std::cout << "Signal	Max		Fix		My" << std::endl;
	for(UInt_t m=0; m<sigCount; m++){
		std::cout << Form(TString(sigName[sigArr[m]])+"	%0.2f +/- %0.2f	%0.2f +/- %0.2f	%0.2f +/- %0.2f", signiMax[m],errorMax[m],signiFix[m],errorFix[m],signi[m],error[m]) << std::endl;;
	}
	*/
	std::cout << "Value ++++++++++++++++++++++" << std::endl;
	std::cout << "Signal	Max	Fix	VV	MaxEr	FixEr	VVEr" << std::endl;
	for(UInt_t m=0; m<sigCount; m++){
		std::cout << Form(TString(sigName[sigArr[m]])+"	%0.3f	%0.3f	%0.3f	%0.3f	%0.3f	%0.3f", signiMax[m],signiFix[m],signi[m],errorMax[m],errorFix[m],error[m]) << std::endl;;
	}


	// Plot
	

}


