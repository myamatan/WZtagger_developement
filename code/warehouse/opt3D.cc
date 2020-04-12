#include<iostream>
#include<fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TCut.h"
#include "TEventList.h"
#include "TF1.h"
#include "TCanvas.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

TString filenameSig = "../input/merged_nv/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged_nv/allMC/allMC."+ch+".root";

float d2cut_max = 1.8;//1.8
float d2cut_min = 0.6;//0.6
float d2cut_d = 0.1;//0.1

float mcutu_max = 126;//126
float mcutu_min = 92;//92
float mcutu_d = 2;//2

float mcutl_max = 86;//86
float mcutl_min = 54;//54
float mcutl_d = 2;//2

float ktdrcut_max = 0.8;//0.8
float ktdrcut_min = 0.1;//0.1
float ktdrcut_d = 0.1;//0.1

void opt3D(){

	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");

	// Optiized value output
	ofstream outf(std::string("../output/cutValue.3D.") + ch.Data() + "." + tagger + "." + sigType + "." + nameOption + std::string(".txt"));
	//ofstream pyoutf(std::string("../output/pycutValue.") + ch.Data() + std::string(".txt"));

	// Cut Fcuntion
	TF1 *fD2;
	TF1 *fML;
	TF1 *fMU;
	std::string wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps["MaxVVSemi"]);
	std::vector<TString> cutWPs = retrieveWPs( wpfile );
	fD2 = new TF1(Form("D2"), cutWPs.at(2), 200, 2500);
	fMU = new TF1(Form("MU"), cutWPs.at(1), 200, 2500);
	fML = new TF1(Form("ML"), cutWPs.at(0), 200, 2500);
	
	UInt_t d2_range = UInt_t( (d2cut_max-d2cut_min)/d2cut_d );
	UInt_t mu_range = UInt_t( (mcutu_max-mcutu_min)/mcutu_d );
	UInt_t ml_range = UInt_t( (mcutl_max-mcutl_min)/mcutl_d );
	UInt_t ktdr_range = UInt_t( (ktdrcut_max-ktdrcut_min)/ktdrcut_d );

	// D2 & Mass window
	TH2F *hb[d2_range][mu_range][ml_range][ktdr_range];
	TH2F *hs[d2_range][mu_range][ml_range][ktdr_range][NmassPoints][2];
	TH2F *hsAll[d2_range][mu_range][ml_range][ktdr_range][2];


	// Pt flat
	TH1F *hpt[NmassPoints];
	for(UInt_t l=0; l<NmassPoints; l++){
		TString nameID = TString(std::to_string(l));
		hpt[l] = new TH1F(nameID+"h", nameID+"h", nPtbins, ptRange);
	}
	TH1F *hptAll= new TH1F("hptAll", "hptAll", nPtbins, ptRange);
	TH1F *hptb= new TH1F("hptb", "hptb", nPtbins, ptRange);

	// Histo initialization (x:VV mass, y:fatjet pt)
	for(UInt_t i=0; i<d2_range; i++){
		for(UInt_t j=0; j<mu_range; j++){
			for(UInt_t k=0; k<ml_range; k++){
				for(UInt_t s=0; s<ktdr_range; s++){
					TString nameID = TString(std::to_string(i) + "-" + std::to_string(j) + "-" + std::to_string(k)+ "-" + std::to_string(s));
					hb[i][j][k][s] = new TH2F(nameID+"hb", nameID+"hb",nbins, rebin21, nPtbins, ptRange);
					for(UInt_t l=0; l<NmassPoints; l++){
						nameID = TString(std::to_string(i) + "-" + std::to_string(j) + "-" + std::to_string(k) + "-" + std::to_string(s) + "-" + std::to_string(l));
						hs[i][j][k][s][l][0] = new TH2F(nameID+"hsa", nameID+"hsa", nbins, rebin21, nPtbins, ptRange);
						hs[i][j][k][s][l][1] = new TH2F(nameID+"hsb", nameID+"hsb", nbins, rebin21, nPtbins, ptRange);
					}
					hsAll[i][j][k][s][0] = new TH2F(nameID+"hsAllb", nameID+"hsAllb", nbins, rebin21, nPtbins, ptRange);
					hsAll[i][j][k][s][1] = new TH2F(nameID+"hsAlla", nameID+"hsAlla", nbins, rebin21, nPtbins, ptRange);
				}
			}
		}
	}

	// KtDR vs D2 significance map 
	TH2F *hmap[nPtbins][mu_range][ml_range];
	for(UInt_t t=0; t<nPtbins; t++){
		for(UInt_t j=0; j<mu_range; j++){
			for(UInt_t k=0; k<ml_range; k++){
				TString nameID = TString(std::to_string(t) + "-" + std::to_string(j) + "-" + std::to_string(k));
				hmap[t][j][k] = new TH2F(nameID+"hmap", nameID+"hmap", d2_range, d2cut_min, d2cut_max, ktdr_range, ktdrcut_min, ktdrcut_max); 
			}
		}
	}

	/////////////////////////////
	// Background
	/////////////////////////////
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
	// Filling
	TCut precut = returnTCut(preSele);
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN(); //long nEventB = trB->GetEntries();
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l)); //trB->GetEntry(l);
		// D2 &&  Mass window
		for(UInt_t i=0; i<d2_range; i++){
			for(UInt_t j=0; j<mu_range; j++){
				for(UInt_t k=0; k<ml_range; k++){
					for(UInt_t s=0; s<ktdr_range; s++){
						if( d2 < d2cut_min+d2cut_d*i &&
						    m < mcutu_max-mcutu_d*j &&
			 			    m > mcutl_min+mcutl_d*k &&
			 			    ktdr < ktdrcut_min+ktdrcut_d*s &&
						    //m < fMU->Eval(pt) &&
						    //m > fML->Eval(pt) &&
						    1
						){
							hb[i][j][k][s]->Fill(vvm, pt, weight*lumi/40.);
						}
					}
				}
			}
		}
		if( l%(int(nEventB/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventB, float(float(l)/nEventB)*100) << std::endl;
	}

	/////////////////////////////
	// Signal
	/////////////////////////////
	std::cout << "Signal filling ....." << std::endl;
	std::map<int, int> sigMap;
	UInt_t sigCount = 0;
	UInt_t sigArr[10];
	SetBranch(trS, ch);
	precut = returnTCut(preSele);
        trS->Draw(">>elist",precut);
        TEventList *elistS = (TEventList*)gROOT->FindObject("elist");
	long nEventS = elistS->GetN(); //long nEventS = trS->GetEntries();
	// For pt flat
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(mc) != sigWeight.end() ) weightS = weight * lumi * sigWeight[mc];
		if( sigMap.find(mc) == sigMap.end() ){
			sigMap[mc] = sigCount; 
			sigArr[sigCount] = mc;
			sigCount++;
		}
		hpt[sigMap[mc]]->Fill(pt, weightS);
		hptAll->Fill(pt, weightS);
	}
	// Filling
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(mc) != sigWeight.end() ) weightS = weight * lumi * sigWeight[mc];
		if( isPtFlat ) weightS *= ptFlat(hptAll, pt);
		// D2 & Mass window
		for(UInt_t i=0; i<d2_range; i++){
			for(UInt_t j=0; j<mu_range; j++){
				for(UInt_t k=0; k<ml_range; k++){
					for(UInt_t s=0; s<ktdr_range; s++){
						hs[i][j][k][s][sigMap[mc]][0]->Fill(vvm, pt, weightS);
						hsAll[i][j][k][s][0]->Fill(vvm, pt, weightS);
						if( d2 < d2cut_min+d2cut_d*i &&
						    m < mcutu_max-mcutu_d*j &&
			 			    m > mcutl_min+mcutl_d*k &&
			 			    ktdr < ktdrcut_min+ktdrcut_d*s &&
						    //m < fMU->Eval(pt) &&
						    //m > fML->Eval(pt) &&
						    1
						){  
							hs[i][j][k][s][sigMap[mc]][1]->Fill(vvm, pt, weightS);
							hsAll[i][j][k][s][1]->Fill(vvm, pt, weightS);
						}
					}
				}
			}
		}
		if( l%(int(nEventS/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}

	// Significance calculation
	std::cout << "Optimized cut searching  ....." << std::endl;
	float d2cut[nPtbins]={0};
	float mcutu[nPtbins]={0};
	float mcutl[nPtbins]={0};
	float ktdrcut[nPtbins]={0};
	float maxSigni=0.0;
	float signi=0.0;
	for(UInt_t i=0; i<nPtbins; i++){
		// D2 & Mass window
		for(UInt_t j=0; j<d2_range; j++){
			for(UInt_t k=0; k<mu_range; k++){
				for(UInt_t l=0; l<ml_range; l++){
					for(UInt_t s=0; s<ktdr_range; s++){
						signi = 0.;
				    		//signi = significance(hsAll[j][k][l][0], hsAll[j][k][l][1], hb[j][k][l], i+1, sigType);
		       		 		for(UInt_t m=0; m<sigCount; m++){
							//int Nptbin = ptSplit[sigName[sigArr[m]]];
							//if( Nptbin==int(i) ) signi += pow(significance(hs[j][k][l][m][0], hs[j][k][l][m][1], hb[j][k][l], i+1, sigType), 2);
							signi += pow(significance(hs[j][k][l][s][m][0], hs[j][k][l][s][m][1], hb[j][k][l][s], i+1, sigType), 2);
						}
						if( signi > maxSigni){
							d2cut[i] = d2cut_min+d2cut_d*j;
							mcutu[i] = mcutu_max-mcutu_d*k;
							mcutl[i] = mcutl_min+mcutl_d*l;
							ktdrcut[i] = ktdrcut_min+ktdrcut_d*s;
							maxSigni = signi;
						}
						hmap[i][k][l]->SetBinContent(s, j, signi);
					}
				}
			}
		}
		maxSigni=0.;
		signi=0.;
	}

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "Pt	:D2	:Upper Mass	: Lower Mass:	KtDR" << std::endl;	
	for(UInt_t i=0; i<nPtbins; i++){
		std::cout << ptRange[i] << "	" << d2cut[i] << "	" << mcutu[i] << "	" << mcutl[i] << "	" << ktdrcut[i] << std::endl;
		if( i==nPtbins-1 ) outf << int((ptRange[i]+ptRange[i+1])/2.) << "," << d2cut[i] << "," << mcutu[i] << "," << mcutl[i] << "," << ktdrcut[i];
		else outf << int((ptRange[i]+ptRange[i+1])/2.) << "," << d2cut[i] << "," << mcutu[i] << "," << mcutl[i] << "," << ktdrcut[i] << ",";
		//pyoutf << int((ptRange[i]+ptRange[i+1])/2.) << "," << d2cut[i] << "," << mcutu[i] << "," << mcutl[i] << std::endl;
	}
	outf.close();

	for(UInt_t i=0; i<nPtbins; i++){
		for(UInt_t k=0; k<mu_range; k++){
			for(UInt_t l=0; l<ml_range; l++){
				TCanvas *cv = new TCanvas("c1","c1",200,10,560,480);
				cv->cd();
				hmap[i][k][l]->Draw("colz");
				cv->Print("../output/sigificanceMap."+ch+"."+tagger+"."+sigType+"."+nameOption+"."+Form("%d-%d-%d",i,k,l)+".png");
			}
		}
	}

}
