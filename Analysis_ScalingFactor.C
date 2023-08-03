//  Analysis_ScalingFactor.c
//
//  Created by Andreas Erhart on 27.07.23.
//  for the analysis of the CMV light output data - light output quantification via global, model-independent scaling parameter

// available data sets (2022/05/05)
// 300K (total acquisition time: 22132sec)
// 800mK (total acquisition time: 44621sec)



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// process the FADC root files (event-by-event) - 800mK (return TH1D)
// for the 800mK measurement: apply NPE calibration (event-by-event) & scaling factor (event-by-event)

TH1D* data_in_hist_800mK(TString fileName, float channelNr, TString histTitle, TString Header, float runtime, double scalingParam, TString outFileName){
    
    TString baseName = "fadc0channel";
    TString chainName = baseName+channelNr;
    
    TH1D *hist_800mK = new TH1D("hist",Header,1000,-50,310);
    
    TChain *ch = new TChain(chainName);
    ch->Add(fileName);
    float nentries = ch->GetEntries();
    
    vector<unsigned int>* acc = new vector<unsigned int>;
    TBranch* b_acc;
    ch->SetBranchAddress("accumulator",&acc,&b_acc);
    for (int i = 0; i<nentries; i++) {
        ch->GetEntry(i);
        // default threshold: cut events below 60NPE
        if ((acc->at(4)-3.57024e5)/72.84 > 60.){
            hist_800mK->Fill(((acc->at(4)-3.57024e5)/72.84)*scalingParam);   // apply pedestal subtraction, NPE calibration (1/72.84) and scaling factor (1*scalingParam)
        }
    }
    
    hist_800mK->SetTitle(Header);
    hist_800mK->Scale(1/runtime);
    
    hist_800mK->SetLineColor(kBlue);
    hist_800mK->SetLineWidth(2);
    
    hist_800mK->GetYaxis()->SetTitle("rate [Hz]");
    hist_800mK->GetXaxis()->SetTitle("NPE");
    
    return hist_800mK;
    
}


//-------It is possible to do it without looping over the events: it will reduce largely the processing time ;)
TH1D* data_in_hist_800mK_fast(TString fileName, float channelNr, TString histTitle, TString Header, float runtime, double scalingParam, TString outFileName){
    
    TFile* filein = new TFile(outFileName, "RECREATE");
    
    TString baseName = "fadc0channel";
    
    TString chainName = baseName+channelNr;
    TChain *ch = new TChain(chainName);
    ch->Add(fileName);
    
    ch->Draw(Form("(accumulator[4]-3.57024e5)/72.84*%f>>hist_800mK(1000,-50,310)", scalingParam), "(accumulator[4]-3.57024e5)/72.84 > 60.");
    TH1D* hist_800mK = (TH1D*) filein->Get("hist_800mK");
    
    hist_800mK->SetTitle(Header);
    hist_800mK->Scale(1./runtime);
    
    hist_800mK->SetLineColor(kBlue);
    hist_800mK->SetLineWidth(2);
    
    hist_800mK->GetYaxis()->SetTitle("rate [Hz]");
    hist_800mK->GetXaxis()->SetTitle("NPE");
    
    return hist_800mK;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// process the FADC root files (event-by-event) - 300K (return TH1D)
// for the 300K measurement: apply NPE calibration (event-by-event), reference spectrum - no scaling applied (set scalingParam = 1)

TH1D* data_in_hist_300K(TString fileName, float channelNr, TString histTitle, TString Header, float runtime, double scalingParam, TString outFileName){
    
    TString baseName = "fadc0channel";
    TString chainName = baseName+channelNr;
    
    TH1D *hist_300K = new TH1D("hist",Header,1000,-50,310);
    
    TChain *ch = new TChain(chainName);
    ch->Add(fileName);
    float nentries = ch->GetEntries();
    
    vector<unsigned int>* acc = new vector<unsigned int>;
    TBranch* b_acc;
    ch->SetBranchAddress("accumulator",&acc,&b_acc);
    
    for (int i = 0; i<nentries; i++) {
        
        ch->GetEntry(i);
        
        // default threshold: cut events below 60NPE
        
        if ((acc->at(4)-3.57019e5)/83.46 > 60.){
            
            hist_300K->Fill(((acc->at(4)-3.57019e5)/83.46)*scalingParam);   // apply pedestal subtraction, NPE calibration (1/83.46) and scaling factor (1*scalingParam)
    
        }
        
    }
    
    
    hist_300K->SetTitle(Header);
    hist_300K->Scale(1/runtime);
    
    hist_300K->SetLineColor(kBlue);
    hist_300K->SetLineWidth(2);
    
    hist_300K->GetYaxis()->SetTitle("rate [Hz]");
    hist_300K->GetXaxis()->SetTitle("NPE");
    
    return hist_300K;
    
}


//-------It is possible to do it without looping over the events: it will reduce largely the processing time ;)
TH1D* data_in_hist_300K_fast(TString fileName, float channelNr, TString histTitle, TString Header, float runtime, double scalingParam, TString outFileName){
    
    TFile* filein = new TFile(outFileName, "RECREATE");
    
    TString baseName = "fadc0channel";
    
    TString chainName = baseName+channelNr;
    TChain *ch = new TChain(chainName);
    ch->Add(fileName);
    
    ch->Draw(Form("(accumulator[4]-3.57019e5)/83.46*%f>>hist_300K(1000,-50,310)", scalingParam), "(accumulator[4]-3.57019e5)/83.46 > 60.");
    TH1D* hist_300K = (TH1D*) filein->Get("hist_300K");
    
    hist_300K->SetTitle(Header);
    hist_300K->Scale(1./runtime);
    
    hist_300K->SetLineColor(kBlue);
    hist_300K->SetLineWidth(2);
    
    hist_300K->GetYaxis()->SetTitle("rate [Hz]");
    hist_300K->GetXaxis()->SetTitle("NPE");
    
    return hist_300K;
}



/*
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//alternatively: gain calibration applied on the spectrum

double ScaleX_300K(Double_t x)
{
    Double_t v;
    //calibration function ugl: y[keV] = 0.92*x[MCA a.u.] + 2.36
    v = ((x/83.46));
    return v;
}

double ScaleX_800mK(Double_t x)
{
    Double_t v;
    //calibration function ugl: y[keV] = 0.92*x[MCA a.u.] + 2.36
    v = ((x/72.84));
    return v;
}

void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t))
{
    if (!a) return; // just a precaution
    if (a->GetXbins()->GetSize())
    {
        // an axis with variable bins
        // note: bins must remain in increasing order, hence the "Scale"
        // function must be strictly (monotonically) increasing
        TArrayD X(*(a->GetXbins()));
        for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
        a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
    }
    else
    {
        // an axis with fix bins
        // note: we modify Xmin and Xmax only, hence the "Scale" function
        // must be linear (and Xmax must remain greater than Xmin)
        a->Set( a->GetNbins(),
               Scale(a->GetXmin()), // new Xmin
               Scale(a->GetXmax()) ); // new Xmax
    }
    return;
}

void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t))
{
    if (!h) return; // just a precaution
    ScaleAxis(h->GetXaxis(), Scale);
    return;
}
*/



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// step 1 - plot the full spectrum

void plotSPECTRA(){
    
    TCanvas *c1 = new TCanvas("cHistos","cHistos", 50, 50, 800, 600);
    gPad->SetGridx(0); gPad->SetGridy(0); gPad->SetLogy();
    
    TH1D* hist_300K = data_in_hist_300K_fast("/Volumes/Data_Drive/Studium/Promotionsstudium/2022-05-05_Prototype_polystyrene_Cryostat/300K/Trigger_SiPM_noCoinc_800mV_thr90_cryostat-prototypePS_300K_*", 1, "hist", "hist", 22132., 1., "test_git");

    
    TH1D* hist_800mK = data_in_hist_800mK_fast("/Volumes/Data_Drive/Studium/Promotionsstudium/2022-05-05_Prototype_polystyrene_Cryostat/800mK/Trigger_SiPM_noCoinc_800mV_thr90_cryostat-prototypePS_800mK_*", 1, "hist", "hist", 44621., 0.61, "test_git");
    
    
    c1->cd(1);
    //ScaleXaxis(hist_300K, ScaleX_300K); // for alternative gain calibration
    hist_300K->SetLineColor(kGreen-3);
    hist_300K->Rebin(4);
    hist_300K->SetLineWidth(2);
    hist_300K->SetStats(kFALSE);
    hist_300K->SetTitle("Full Spectrum");
    hist_300K->Draw("HIST");
    
    hist_800mK->SetLineColor(kBlue+2);
    hist_800mK->Rebin(4);
    hist_800mK->SetLineWidth(2);
    hist_800mK->SetStats(kFALSE);
    hist_800mK->SetTitle("Full Spectrum");
    hist_800mK->Draw("HISTsame");
    
    hist_300K->GetXaxis()->SetRangeUser(0,310);
    hist_300K->GetYaxis()->SetRangeUser(5e-5,1);
    
    TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->SetBorderSize(0);
    leg->AddEntry(hist_300K,"Full Spectrum, 300K","l");
    leg->AddEntry(hist_800mK,"Full Spectrum, 0.8K, scal.param. =0.62","l");
    leg->Draw();
    
    
    //c1->Print(outFileName+".png");
    //hist_300K->SaveAs("Plot_Trigger_disc75mm_noCoinc_fullSpectrum.root");
}


Double_t chi2_from_TH1Ds(TH1D* histo_1, TH1D* histo_2, Double_t PE_min, Double_t PE_max){
    int index_PE_min = histo_1->FindBin(PE_min);
    int index_PE_max = histo_1->FindBin(PE_max);
    
    Double_t chi2, Content_1, Content_2, Error_1, Error_2;
    for (int i=index_PE_min; i<index_PE_max+1; i++){
        Content_1 = histo_1->GetBinContent(i+1);
        Error_1 = histo_1->GetBinError(i+1);
        Content_2 = histo_2->GetBinContent(i+1);
        Error_2 = histo_2->GetBinError(i+1);
        if (Content_1>0 || Content_2>0){
            chi2 += ((Content_1-Content_2)*(Content_1-Content_2))/(Error_1*Error_1+Error_2*Error_2);
        }
    }
    return chi2;
}

int Non_zero_NBins(TH1D* histo_1, TH1D* histo_2, Double_t PE_min, Double_t PE_max){
    int index_PE_min = histo_1->FindBin(PE_min);
    int index_PE_max = histo_1->FindBin(PE_max);
    
    int Non_zero;
    Double_t Content_1, Content_2, Error_1, Error_2;
    for (int i=index_PE_min; i<index_PE_max+1; i++){
        Content_1 = histo_1->GetBinContent(i+1);
        Error_1 = histo_1->GetBinError(i+1);
        Content_2 = histo_2->GetBinContent(i+1);
        Error_2 = histo_2->GetBinError(i+1);
        if (Content_1>0 || Content_2>0){
            Non_zero+=1;
        }
    }
    return Non_zero;
}

///// Fit function
///
ROOT::Fit::FitResult Fit_Data_Simu(TString File_800mK= "../Files_800mK/Trigger_SiPM_noCoinc_800mV_thr90_cryostat-prototypePS_800mK_*", TString File_300K="../Files_300K/Trigger_SiPM_noCoinc_800mV_thr90_cryostat-prototypePS_300K_*") {

    Double_t PE_min=45.;
    Double_t PE_max=300.;
    
    TH1D* h_300K = data_in_hist_300K_fast(File_300K, 1, "hist_300K", "hist_300K", 1976.6, 1., "file300k.root");
    h_300K->Rebin(10);
    TH1D* h_800mK = new TH1D();
    
    //cf. https://root.cern.ch/doc/master/fitCircle_8C.html
    auto chi2Function_Fit = [&](const Double_t *par) {
        //minimisation function computing the sum of squares of chi2
        Double_t Scaling = par[0];
        
        //Scale 800mK
        h_800mK = data_in_hist_800mK_fast(File_800mK, 1, "hist_800mK", "hist_800mK", 7475.1, Scaling, "file800mK.root");
        h_800mK->Rebin(10);
        /// ------------ Total Chi2
        double chi2 = chi2_from_TH1Ds(h_300K, h_800mK, PE_min, PE_max);
        return chi2;
    };
    
    // wrap chi2 function in a function object for the fit
    // 4 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn_Fit(chi2Function_Fit, 1);
    ROOT::Fit::Fitter fitter_Fit;
    
    // Initial fit parameters
    Double_t pStart_Fit[1] = {0.65};

    fitter_Fit.SetFCN(fcn_Fit, pStart_Fit);
    fitter_Fit.Config().ParSettings(0).SetLimits(0.5, 0.7);
    fitter_Fit.Config().ParSettings(0).SetName("Sacling Factor");
    
    // do the fit
    cout<<"Fit starts"<<endl;
    bool ok_Fit = fitter_Fit.FitFCN();
    cout<<"Fit End"<<endl;
    
    // test fit
    if(!ok_Fit) {
        Error("fit", "%s", "Fit failed");
    }
    
    const ROOT::Fit::FitResult &result_Fit = fitter_Fit.Result();
    //    result_Fit.Print(std::cout);
    
    Double_t Fitted_Scaling = result_Fit.Parameter(0);
    Double_t Fitted_Scaling_Error = result_Fit.ParError(0);
    
    double Chi2_Fit = result_Fit.MinFcnValue();
    int N_dof = Non_zero_NBins(h_300K, h_800mK, PE_min, PE_max) - 1;
    
    cout << "\n**********************************************************************"<<endl;
    cout<<" Fit between "<<PE_min<<" and "<<PE_max<<" PE"<<endl;
    cout << "Chi2/dof = " << Chi2_Fit << "/" << N_dof << " = " << Chi2_Fit/N_dof<<endl;
    cout <<" Scaling Factor = " << Fitted_Scaling << " +- " << Fitted_Scaling_Error <<endl;
    cout << "**********************************************************************"<<endl;
    
    h_300K->Draw();
    h_800mK->SetLineColor(kRed);
    h_800mK->Draw("SAME");
    return result_Fit;
}

// Residual histogram
TH1D* Residuals_histos(TH1D* histo_data, TH1D* Model, Double_t Min, Double_t Max, bool pc = true){
    TH1D* hRes = (TH1D*) histo_data->Clone();
    double Value_i, E_i, Diff_i, Error_i;
    double chi2=0;
    for (int i = 1; i<hRes->GetNbinsX(); i++){
        Value_i = hRes->GetBinContent(i);
        E_i = hRes->GetBinCenter(i);
        if (Value_i>0&&E_i>Min&&E_i<Max) {
            Diff_i = histo_data->GetBinContent(i)-Model->GetBinContent(i);
            Error_i = sqrt(histo_data->GetBinError(i)*histo_data->GetBinError(i)+Model->GetBinError(i)*Model->GetBinError(i));
            if (pc) {
                hRes->SetBinContent(i, Diff_i/Value_i*100);
                hRes->SetBinError(i, Error_i/Value_i*100);
            }
            else {
                hRes->SetBinContent(i, Diff_i/Error_i);
                hRes->SetBinError(i, 1);
            }
            chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
        }
        else {
            hRes->SetBinContent(i, 0);
            hRes->SetBinError(i, 0);}
    }

    cout<<"\nChi2/ndf = "<<chi2<<"/"<<(Max-Min)/Model->GetBinWidth(1)<<" = "<<chi2/((Max-Min)/Model->GetBinWidth(1))<<endl;
    hRes->SetLineColor(Model->GetLineColor());
    hRes->SetTitle("");
    if (pc){
        hRes->GetYaxis()->SetRangeUser(-100, 100);
        hRes->GetYaxis()->SetTitle("Residuals [%]");
    }
    else {
        hRes->GetYaxis()->SetRangeUser(-20, 20);
        hRes->GetYaxis()->SetTitle("Residuals [#sigma]");
    }
    
    
    return hRes;
}
