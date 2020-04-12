TString ch = "1Lep";
TString model = "ggH"; // Old
//TString model = "WZ";   // New (2019.02.25)
TString tagger = "Z";
TString sigType = "logli"; // "mass", "eff", "logli"
TString preSele = "beforeBosonTag"; // "beforeTargetCut", "beforeBosonTag", "forPlot"
TString nameOption = "20190308_signal";// logliFit// 5 //ReferenceFix // TuneFit //_ConstantMassWindow
//TString nameOption = "20190130_140fb";// logliFit// 5 //ReferenceFix // TuneFit //_ConstantMassWindow
bool isPtFlat = false;
float lumi = 14E+1; //[fb-1]
bool isPlotMarker = true;

bool isMCSmooth = true;
float fitRangeMin =1000; //800 //1000
float smoothXmin = 1600;  //2000

// Drawing
float func_xMax = 2500; // 4000 for fullrange;

// Combination fit (e.g VVSemi + VVJJ)
bool isCombineFit = false;
float connectBin = 1100;

// Cut function connection
bool isSubset = true;

TString sigModel = "_" + model;
const UInt_t NmassPoints = 30;// 5 

float k_ymax;

//const int nbins=18;
//double rebin21[nbins+1]={0,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2500,6000};
//const int nbins=21;
//double rebin21[nbins+1]={0,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2500,3000,4000,5000,6000};

// W
const int nbins=16;
double rebin21[nbins+1]={0,400,600,800,1000,1200,1400,1600,1800,2000,2300,2600,3000,3500,4000,5000,6000};
// Z
//const int nbins=17;
//double rebin21[nbins+1]={0,400,600,800,1000,1200,1400,1600,1800,2000,2300,2600,3000,3400,3800,4200,5000,6000};

//const int nPtbins=16;
//double ptRange[nPtbins+1]={200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1500, 1700, 1900, 2100, 2500};

//const int nPtbins=8; // Final version
//double ptRange[nPtbins+1]={200,300,500,700,900,1300,1700,2100,2500};
const int nPtbins=7; // Efficiency (NOTE)
double ptRange[nPtbins+1]={200,300,450,700,1100,1550,1900,2500};
const int nPtbinsDraw=9; // Final version
double ptRangeDraw[nPtbinsDraw+1]={200,300,500,700,900,1300,1700,2100,2500,func_xMax};

// List of available working points
std::map<TString, TString> wps=
{
  {"MaxVVJJ",		  "Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_MC16a_20171201.dat"},
  {"MaxVVJJ_3D",	  "Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_MC16d_20181111.dat"},
  {"FixVVJJ",		  "Tagger_AntiKt10TrackCaloClusterTrimmed_FixedSignalEfficiency50_MC16a_20171201.dat"},
  {"Fix50VVJJ_LCT",	  "Tagger_AntiKt10LCTopoTrimmed_FixedSignalEfficiency50_MC15c_20161219.dat"},
  {"Fix80VVJJ_LCT",	  "Tagger_AntiKt10LCTopoTrimmed_FixedSignalEfficiency80_MC15c_20161219.dat"},
  {"MaxVVSemi",		  "Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_VVSemilep_MC16d_20181204.dat"},
  {"MaxVVSemi_v3",        "Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_VVSemilep_MC16d_20181217.txt"},
  {"MaxVVSemi_v4",        "Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_VVSemilep_MC16d_20190117.dat"},
  {"MaxVVSemi_v5",        "Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_VVSemilep_MC16d_20190130.dat"},
  {"MaxVVSemi_Flat05",    "Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_VVSemilep_pFlat05_MC16d_20181204.dat"},
  {"MaxVVJJSemiCombined", "Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_VVSemilep_VVJJ_combined_20180901.dat"},
  //{"MaxVVSemi_test",	  "Tagger_AntiKt10TrackCaloClusterTrimmed_"+ch+"_"+sigType+"_"+nameOption+".txt"},
  {"MaxVVSemi_test",	  "Tagger_AntiKt10TrackCaloClusterTrimmed_1Lep_"+sigType+"_"+nameOption+".txt"},
  {"Flat15",		  "Tagger_AntiKt10TrackCaloClusterTrimmed_Flat15.dat"},
  {"Flat10",		  "Tagger_AntiKt10TrackCaloClusterTrimmed_Flat10.dat"},
  {"Flat18",		  "Tagger_AntiKt10TrackCaloClusterTrimmed_Flat18.dat"},
  {"Flat08",		  "Tagger_AntiKt10TrackCaloClusterTrimmed_Flat08.dat"},
  {"preMaxVVSemi",	  "Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_VVSemileptonic_20181120.txt"},
  {"Fix50VVSemi",	  "Tagger_AntiKt10TrackCaloClusterTrimmed_FixedSignalEfficiency50_test_VVSemilep_MC16d_20181209.dat"},
  {"Fix80VVSemi",	  "Tagger_AntiKt10TrackCaloClusterTrimmed_FixedSignalEfficiency80_test_VVSemilep_MC16d_20181209.dat"}
};

/////////////////////////////////////////////////////////////////////////////////////////
// Combination of WPs ---> {"Strategy Name", {"D2 cut","UpperM cut","LowerM cut"}}
// HP
std::map<TString, std::vector<TString>> combs=
{ 
  //{"FixedSignalEfficiency50 VVJJ",                {"FixVVJJ","FixVVJJ","FixVVJJ"}}
  //{"VVJJ 3D HP(D2 passed) + LP(D2 failed)",                        {"MaxVVJJ_3D","MaxVVJJ_3D","MaxVVJJ_3D"}},
  //{"Default HP(D2 passed) + LP(D2 failed)",                   {"MaxVVJJSemiCombined","MaxVVJJSemiCombined","MaxVVJJSemiCombined"}},
  //{"MaxSignificance VVSemi-lep",                   {"MaxVVSemi","MaxVVSemi","MaxVVSemi"}},
  //{"Pre-test  HP(D2 passed) + LP(D2 failed)",                   {"preMaxVVSemi","preMaxVVSemi","preMaxVVSemi"}},
  //{"EfficiencyBased   HP(D2 passed) + LP(D2 failed)",                   {"MaxVVSemi","MaxVVSemi","MaxVVSemi"}},
  //{"D2:ver1, Mass:ver2",                   {"preMaxVVSemi","MaxVVSemi","MaxVVSemi"}},
  //{"D2:ver2, Mass:ver1",                   {"MaxVVSemi","preMaxVVSemi","preMaxVVSemi"}},
  //{"VVSemi optimization",                   	{"MaxVVSemi_test","MaxVVSemi_test","MaxVVSemi_test"}}
  //{"SignificanceBased HP(D2 passed) + LP(D2 failed)",                   {"MaxVVSemi_v5","MaxVVSemi_v5","MaxVVSemi_v5"}}
  {"VV Semileptonic tagger",                   {"MaxVVSemi_v5","MaxVVSemi_v5","MaxVVSemi_v5"}},
  {"MaxSignificance VVJJ",                        {"MaxVVJJ","MaxVVJJ","MaxVVJJ"}}
  //{"Harmonizaed HP(D2 passed) + LP(D2 failed)",                   {"MaxVVSemi_v4","MaxVVSemi_v4","MaxVVSemi_v4"}}
};

// LP
std::map<TString, std::vector<TString>> combs_LP=
{ 
  //{"FixedSignalEfficiency50 VVJJ",                {"FixVVJJ","FixVVJJ","FixVVJJ"}},
  //{"MaxSignificance VVJJ",                        {"MaxVVJJ","MaxVVJJ","MaxVVJJ"}},
  //{"MaxSignificance VVJJ 3D",                        {"Flat15","MaxVVJJ_3D","MaxVVJJ_3D"}},
  //{"pre",                       {"Flat15","MaxVVJJ_3D","MaxVVJJ_3D"}},
  //{"Default",                        {"Flat15","MaxVVJJ_3D","MaxVVJJ_3D"}},
  //{"MaxSignificance VVSemi-lep",                   {"MaxVVSemi","MaxVVSemi","MaxVVSemi"}},
  //{"HP+Lp VVSemi-lep ",                   {"Flat15","preMaxVVSemi","preMaxVVSemi"}},
  //{"HP+Lp VVSemi-lep 2",                   {"Flat15","MaxVVSemi","MaxVVSemi"}},
  //{"TCC HP(preFix50)+LP(preFix80)",                   	{"Fix80VVSemi","Fix80VVSemi","Fix80VVSemi"}}
  //{"TCC HP(preFix50)+LP(preFix80)-2",                   	{"Fix80VVSemi","Fix80VVSemi","Fix80VVSemi"}}
  {"VVSemi optimization",                   	{"Flat15","MaxVVSemi","MaxVVSemi"}},
  {"VaaaaVSemi optimization",                   	{"Flat15","MaxVVSemi","MaxVVSemi"}}
  //{"Harmoni",                   	{"Flat15","MaxVVSemi","MaxVVSemi"}}
};

float vvm, d2, pt, m, ktdr, mu;
float weight;
//Float_t nbtag;
int nbtag;
int pdgid;
//UInt_t dsid; // new
int dsid;  // old
int mc;
int wtag, ztag;
std::string *re = 0;
Bool_t pass_OneLepton;
Bool_t pass_Trigger;
Bool_t pass_SigJetORFatJet;
Bool_t pass_isVBF;
Bool_t pass_Merg_AtLeastOneFatJet;
Bool_t pass_Merg_MET_100GeV;
Bool_t pass_Merg_METoverWpT;
Bool_t pass_Merg_pTW_200GeV;
Bool_t pass_Merg_GGF_pTW_MlvJ;
Bool_t pass_Merg_GGF_pTJ_MlvJ;
Bool_t pass_Merg_bVeto;
Bool_t passZtagMass;
Bool_t passWtagMass;

bool passPreSelection = false;
TString preCutfor1Lep = "Pass_OneLepton && Pass_Trigger && Pass_SigJetORFatJet && !Pass_isVBF && Pass_Merg_AtLeastOneFatJet && Pass_Merg_MET_100GeV && Pass_Merg_METoverWpT && Pass_Merg_pTW_200GeV && Pass_Merg_GGF_pTJ_MlvJ && Pass_Merg_GGF_pTW_MlvJ && Pass_Merg_bVeto && abs(weight)<10";
TString preCutfor1Lep_loose = "Pass_OneLepton && Pass_Trigger && Pass_SigJetORFatJet && !Pass_isVBF && Pass_Merg_AtLeastOneFatJet && Pass_Merg_MET_100GeV && Pass_Merg_pTW_200GeV && abs(weight)<10";
TString preCutfor1Lep_veryloose = "Pass_OneLepton && Pass_Trigger && Pass_Merg_AtLeastOneFatJet && abs(weight)<10";

//TString preCutfor1Lep = "Pass_OneLepton && Pass_Trigger && Pass_SigJetORFatJet && !Pass_isVBF && Pass_Merg_AtLeastOneFatJet && Pass_Merg_MET_100GeV && Pass_Merg_METoverWpT && Pass_Merg_pTW_200GeV && Pass_Merg_GGF_pTJ_MlvJ && Pass_Merg_GGF_pTW_MlvJ && Pass_Merg_bVeto && !(fatjet_pt>1900 && abs(weight)>0.04) && !(fatjet_pt<1500 && fatjet_pt>1100 && abs(weight)>0.15) && !(fatjet_pt>2000 && weight<-0.01) && !(fatjet_pt>1150&&fatjet_pt<1600&&weight>0.25)";


std::vector<Color_t> colorList={
	kRed-3,
	kBlue-3,
	kGreen-3,
	kOrange-3,
	kBlack,
	kBlue-9,
	kOrange-9,
	kCyan-9,
	kGreen-9,
	kBlue,
	kRed,
	kOrange,
};
std::vector<Color_t> colorListSig={
	kBlue-3,
	kRed-3,
	kGreen-3,
	kOrange-3,
	kAzure-3,
	kYellow-3,
	kOrange-3,
};

TCut returnTCut(TString step="forPlot", TString CH=ch)
{
	TCut precut = "";
	if( step=="beforeTargetCut" ){
		if( CH!="1Lep" ){
			//precut = "fatJet_nBTags"+sigModel+">=-9 && "+ // new
			precut = "fatJet_D2"+sigModel+">0 &&" +
				 "fatJet_m"+sigModel+">=50 &&" +
				 "fatJet_m"+sigModel+"<=150 &&" + 
				 "fatJet_D2"+sigModel+"<=5";
		}else{
			//precut = "fatjet_D2>0 & fatjet_nBTags>=-99 & fatjet_D2<=3 & fatjet_m>=50 & fatjet_m<=150 & "+
			precut = "fatjet_D2>0 & fatjet_D2<=5 & fatjet_m>=50 & fatjet_m<=150 & "+
				 preCutfor1Lep;
		}
	}
	else if( step=="beforeBosonTag" ){
		if( CH!="1Lep" ){
			//precut = "fatJet_nBTags"+sigModel+">=-99 && "+ // New
		  	precut = "fatJet_m"+sigModel+">=50 &&" +
				 "fatJet_D2"+sigModel+">0";
		}else{
			//precut = "fatjet_D2>0 && fatjet_nBTags>=-99 && fatjet_m >=50 &&" +
			precut = "fatjet_D2>0 && fatjet_m >=50 &&" + // Old
				 preCutfor1Lep;
		}
	}
	else if( step=="forPlot" ){
		if( CH!="1Lep" ){
			precut = "fatJet_nBTags"+sigModel+">=-99 & "+
				 "fatJet_D2"+sigModel+">0 &" +
				 "fatJet_m"+sigModel+">=50 &" +
				 "fatJet_m"+sigModel+"<=150 &" + 
				 "fatJet_D2"+sigModel+"<=5";
		}else{
			precut = "fatjet_D2>0 & fatjet_nBTags>=-99 & fatjet_D2<=5 & fatjet_m>=50 & fatjet_m<=150 & "+
				 preCutfor1Lep;
		}
	}
	return precut;

}

std::vector<TString> retrieveWPs(std::string file){
	std::vector<TString> wps={"","","",""};
	ifstream ifs(file);
	std::string line;
	while( getline(ifs, line) ){
		TString str = TString(line);
		if( str.Index("MassCutLow")!=-1 ){
			str = str(str.Index(":")+1, str.Length());
			wps.at(0) = str;
		}
		else if( str.Index("MassCutHigh")!=-1 ){
			str = str(str.Index(":")+1, str.Length());
			wps.at(1) = str;
		}
		else if( str.Index("D2Cut")!=-1 ){
			str = str(str.Index(":")+1, str.Length());
			wps.at(2) = str;
		}
		else if( str.Index("KtDRCut")!=-1 ){
			str = str(str.Index(":")+1, str.Length());
			wps.at(3) = str;
		}
		
	}
	return wps;
}
	

void SetBranch(TTree* tr, TString ch="0Lep"){
	if( ch=="0Lep"){
	  tr->SetBranchAddress("vvJ_mt"+sigModel,&vvm);
	  tr->SetBranchAddress("fatJet_m"+sigModel,&m);
	  tr->SetBranchAddress("fatJet_pt"+sigModel,&pt);
	  tr->SetBranchAddress("fatJet_D2"+sigModel,&d2);
	  tr->SetBranchAddress("mu"+sigModel,&mu);
	  tr->SetBranchAddress("fatJet_nBTags"+sigModel,&nbtag); // New
	  //tr->SetBranchAddress("fatJet_nBtag"+sigModel,&nbtag); // Old
	  tr->SetBranchAddress("fatJet_WTag50"+sigModel,&wtag);
	  tr->SetBranchAddress("fatJet_ZTag50"+sigModel,&ztag);
	  tr->SetBranchAddress("Region"+sigModel,&re);
	  //tr->SetBranchAddress("DSID",&dsid); //New
	  tr->SetBranchAddress("DSID"+sigModel,&dsid);
	  //tr->SetBranchAddress("weight",&weight);//New
	  tr->SetBranchAddress("weight"+sigModel,&weight);
	}else if( ch=="2Lep" ){
	  tr->SetBranchAddress("llJ_m"+sigModel,&vvm);
	  tr->SetBranchAddress("fatJet_m"+sigModel,&m);
	  tr->SetBranchAddress("fatJet_pt"+sigModel,&pt);
	  tr->SetBranchAddress("fatJet_D2"+sigModel,&d2);
	  tr->SetBranchAddress("fatJet_KtDR"+sigModel,&ktdr);
	  tr->SetBranchAddress("weight",&weight);
	  tr->SetBranchAddress("mu"+sigModel,&mu);
	  //tr->SetBranchAddress("fatJet_nBTags"+sigModel,&nbtag);
	  tr->SetBranchAddress("fatJet_nBtag"+sigModel,&nbtag); // Old
	  tr->SetBranchAddress("fatJet_truthID"+sigModel,&pdgid);
	  tr->SetBranchAddress("fatJet_WTag50"+sigModel,&wtag);
	  tr->SetBranchAddress("fatJet_ZTag50"+sigModel,&ztag);
	  tr->SetBranchAddress("Region"+sigModel,&re);
	  //tr->SetBranchAddress("MCChannelNumber",&mc);
	  tr->SetBranchAddress("DSID"+sigModel,&dsid);
	}else if( ch=="1Lep" ){
	  tr->SetBranchAddress("X_boosted_m",&vvm);
	  tr->SetBranchAddress("fatjet_m",&m);
	  tr->SetBranchAddress("fatjet_D2",&d2);
	  tr->SetBranchAddress("fatjet_pt",&pt);
	  tr->SetBranchAddress("fatjet_nBTags",&nbtag);
	  //tr->SetBranchAddress("fatjet_nBtag",&nbtag); // OLd
	  //tr->SetBranchAddress("MCChannelNumber",&mc);
	  tr->SetBranchAddress("Pass_OneLepton",&pass_OneLepton);
	  tr->SetBranchAddress("Pass_Trigger",&pass_Trigger);
	  tr->SetBranchAddress("Pass_SigJetORFatJet",&pass_SigJetORFatJet);
	  tr->SetBranchAddress("Pass_isVBF",&pass_isVBF);
	  tr->SetBranchAddress("Pass_Merg_AtLeastOneFatJet",&pass_Merg_AtLeastOneFatJet);
	  tr->SetBranchAddress("Pass_Merg_MET_100GeV",&pass_Merg_MET_100GeV);
	  tr->SetBranchAddress("Pass_Merg_METoverWpT",&pass_Merg_METoverWpT);
	  tr->SetBranchAddress("Pass_Merg_pTW_200GeV",&pass_Merg_pTW_200GeV);
	  tr->SetBranchAddress("Pass_Merg_GGF_pTW_MlvJ",&pass_Merg_GGF_pTW_MlvJ);
	  tr->SetBranchAddress("Pass_Merg_GGF_pTJ_MlvJ",&pass_Merg_GGF_pTJ_MlvJ);
	  tr->SetBranchAddress("Pass_Merg_bVeto",&pass_Merg_bVeto);
	  tr->SetBranchAddress("Pass_WTaggerMassCut",&passWtagMass);
	  tr->SetBranchAddress("Pass_ZTaggerMassCut",&passZtagMass);
	  tr->SetBranchAddress("DSID",&dsid);
	  tr->SetBranchAddress("weight",&weight);
	  tr->SetBranchAddress("mu",&mu); //mu or AvgMu
	}
}

float recom_sig(float s, float b, float be){

        float n = s + b;
        float z = 0;
        if( be==0. ){
                z += n*log(n/b);
                z -= (n-b);
                z *= 2.;
                z = sqrt(z);
                if( n>=b ) z *= 1.;
                else z *= -1.;
        }else{
                z += n*log( (n*(b+be)) / (b*b+n*be) );
                z -= b*b/be*log( 1 + (be*(n-b))/(b*(b+be)) );
                z *= 2.;
                z = sqrt(z);
                if( n>=b ) z *= 1.;
                else z *= -1.;
        }
	if( b==0 ) z = 0;
        return z;

}


// Asysmptotic significance New
float significance(TH1F *h00, TH1F *h01, TH1F *h1, TString type="eff", TString value="value") // type="eff" or "mass"
{
	
	if( type=="mass" ){
		float significance=0.0;
		float error=0.0;
		for(int i=1; i<=h01->GetNbinsX(); i++){

        	        float v0 = h01->GetBinContent(i);
        	        float v1 = h1->GetBinContent(i);
        	        float v0e = h01->GetBinError(i);
        	        float v1e = h1->GetBinError(i);

			float com_par = 0.5 * pow(2*((v0+v1)*log(1+v0/v1)-v0), -1/2.);
                        float s_part = 2 * log(1+v0/v1);
                        float b_part = 2 * (log(1+v0/v1)-v0/v1);
                        float sig_bin_er = sqrt( pow(com_par*s_part*v0e, 2)+pow(com_par*b_part*v1e,2) );

			//float sig_bin = float( 2*( (v0+v1)*log(1+v0/v1)-v0 ) );
			float sig_bin = pow(recom_sig(v0, v1, v1e*v1e),2);
			if( isnan(sig_bin) || isinf(sig_bin) || sig_bin < 0 ) significance += 0.0;
			else significance += sig_bin;
			if( isnan(sig_bin_er) || isinf(sig_bin_er) || sig_bin_er < 0 ) error += 0.0;
			else error += pow(sig_bin_er,2);

			if( v1==0 ) std::cout << rebin21[i-1] << ":Empty bin:" << v0 << "," << v1 << "," << sig_bin << std::endl;
			//if( isnan(sig_bin) || isinf(sig_bin) ) std::cout << rebin21[i-1] << ":Wrong filling:" << v0 << "," << v1 << "," << sig_bin << std::endl;
			//if( value!="error" ) std::cout << rebin21[i-1] << ":signal:" << v0 << ", backgrond:" << v1 << ", significance:" << sig_bin << std::endl;
        	}
		if( value=="error") return sqrt(error);
		//if( value=="error") return 0.0;
		else return sqrt(significance);
	}
	else if( type=="eff" ){
        	float effS = h01->Integral()/h00->Integral();
        	float sumB = h1->Integral();
		float effSig = effS / (3./2+sqrt(sumB));
		return effSig;
	}
	return 0.;
	
}


// Asysmptotic significance new
float significance(TH2F *h00, TH2F *h01, TH2F *h1, int ptBin, TString type="eff", TString value="value")
{

	if( type=="mass" ){
		float significance=0.0;
		float error=0.0;
		for(int i=1; i<=h01->GetNbinsX(); i++){
	
	                float v0 = h01->GetBinContent(i, ptBin);
	                float v1 = h1->GetBinContent(i, ptBin);
	                float v0e = h01->GetBinError(i, ptBin);
	                float v1e = h1->GetBinError(i, ptBin);
		
			v0 = v0 < 0 ? 0.0 : v0;
			v1 = v1 < 0 ? 0.0 : v1;

			float com_par = 0.5 * pow(2*((v0+v1)*log(1+v0/v1)-v0), -1/2.);
                        float s_part = 2 * log(1+v0/v1);
                        float b_part = 2 * (log(1+v0/v1)-v0/v1);
                        float sig_bin_er = sqrt( pow(com_par*s_part*v0e, 2)+pow(com_par*b_part*v1e,2) );
		
			//float sig_bin = float( 2*( (v0+v1)*log(1+v0/v1)-v0 ) );
			float sig_bin = pow(recom_sig(v0,v1,v1e*v1e),2);
			if( isnan(sig_bin) || isinf(sig_bin) || sig_bin < 0 ) significance += 0.0;
			else significance += sig_bin;
			if( isnan(sig_bin_er) || isinf(sig_bin_er) || sig_bin_er < 0 ) error += 0.0;
			else error += pow(sig_bin_er,2);
	        }
		if( value=="error" ) return sqrt(error);
		else return sqrt(significance);
	}
	else if( type=="eff" ){
		float sumSb=0.0;
		float sumSa=0.0;
		float sumBa=0.0;
		for(int i=1; i<=h00->GetNbinsX(); i++){
	                sumSb += h00->GetBinContent(i, ptBin);
	                sumSa += h01->GetBinContent(i, ptBin);
	                sumBa += h1->GetBinContent(i, ptBin);
		}
        	float effS = sumSa/sumSb;
		float effSig = effS / (2/3.+sqrt(sumBa));// VV->JJ : 1.5
		return effSig;
	}
	else if( type=="logli" ){
		float v0=0.0;
		float v1=0.0;
		float v1er=0.0;
		for(int i=1; i<=h01->GetNbinsX(); i++){
	                v0 += h01->GetBinContent(i, ptBin);
	               	v1 += h1->GetBinContent(i, ptBin);
	               	v1er += pow(h1->GetBinError(i, ptBin),2);
		}
		//float significance = float( 2*( (v0+v1)*log(1+v0/v1)-v0 ) );
		float significance = recom_sig(v0, v1, v1er);
		//std::cout << significance << "," << v0 << "," << v1 << "," << v1er << std::endl;
		return significance;
	}
	return 0.;
		
}


// Pt flat
float ptFlat(TH1F *h, float pt)
{
	float ave = h->Integral()/h->GetNbinsX();
	float scale = ave/h->GetBinContent(h->FindBin(pt));
	return scale;
}

// Efficiency
void setEfficiency(TH1F *h0, TH1F *h1, TH1F*he, TString type="efficiency") // type="eff" or "mass"
{
	
	for(int i=1; i<=h0->GetNbinsX(); i++){

	        float v0 = h0->GetBinContent(i);
	        float v1 = h1->GetBinContent(i);
	        float v0e = h0->GetBinError(i);
	        float v1e = h1->GetBinError(i);
	
		v0 = v0 < 0 ? 0.0 : v0;
		v1 = v1 < 0 ? 0.0 : v1;
		
		float eff = (v0!=0) ? v1/v0 : -99;
		if( type == "rejection" ) eff = (v0!=0) ? 1/(v1/v0) : -99;
		float effEr = 0;
		if( type=="efficiency" ){
			effEr = (v0!=0) ? sqrt( pow(v1/v0/v0*v0e,2) + pow(1./v0*v1e,2) ) : -99;
		}else{
			effEr = (v1!=0) ? sqrt( pow(v0/v1/v1*v1e,2) + pow(1./v1*v0e,2) ) : -99;
		}
			
		he->SetBinContent(i, eff);
		he->SetBinError(i, effEr);
		
	}
	
}

// Ratio
void setRatio(TH1F *hn, TH1F *hd, TH1F*hr) // type="eff" or "mass"
{
	
	for(int i=1; i<=hd->GetNbinsX(); i++){

	        float vn = hn->GetBinContent(i);
	        float vd = hd->GetBinContent(i);
	        float vne = hn->GetBinError(i);
	        float vde = hd->GetBinError(i);
	
		vn = vn < 0 ? 0.0 : vn;
		vd = vd < 0 ? 0.0 : vd;
		
		float ratio = (vd!=0) ? vn/vd : -99;
		float ratioEr = 0;
		ratioEr = (vd!=0) ? sqrt( pow(vn/vd/vd*vde,2) + pow(1./vd*vne,2) ) : -99;
			
		hr->SetBinContent(i, ratio);
		hr->SetBinError(i, ratioEr);
		
	}
	
}


void setHisto(TH1F* h, Color_t fillColor, Color_t lineColor, float xMin, float xMax, float yMin, float yMax, TString value, int lineStyle=1, bool isFill=false){

	h->SetStats(0);
	h->SetTitle("");
	h->GetYaxis()->SetTitleSize(0.04);
	h->GetYaxis()->SetTitleOffset(1.24);
	h->GetXaxis()->SetTitle(value);
	h->GetXaxis()->SetRangeUser(xMin,xMax);
	h->GetYaxis()->SetTitle("A.U.");
	h->GetYaxis()->SetRangeUser(yMin,yMax);
	h->SetLineColor(lineColor);
	h->SetLineWidth(2);
	h->SetMarkerColor(lineColor);
	h->SetMarkerStyle(8);
	h->SetMarkerSize(0.8);
	h->GetYaxis()->SetMaxDigits(2);
	h->SetLineStyle(lineStyle);
	if( isFill ){
		h->SetFillColor(fillColor);
		h->SetFillStyle(3002);
	}

}

/*
bool dsidMatch(int DSID, string region="0Lep_Z"){

	bool isMatch=false;
	for(int i=0; i<int( sigDSID[region].size() ); i++){
		if( DSID==sigDSID[region].at(i) ) isMatch = true;
	}
	return isMatch;

}
*/
