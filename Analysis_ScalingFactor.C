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
    
    TH1D* hist_300K = data_in_hist_300K("/Volumes/Data_Drive/Studium/Promotionsstudium/2022-05-05_Prototype_polystyrene_Cryostat/300K/Trigger_SiPM_noCoinc_800mV_thr90_cryostat-prototypePS_300K_*", 1, "hist", "hist", 22132., 1., "test_git");

    
    TH1D* hist_800mK = data_in_hist_800mK("/Volumes/Data_Drive/Studium/Promotionsstudium/2022-05-05_Prototype_polystyrene_Cryostat/800mK/Trigger_SiPM_noCoinc_800mV_thr90_cryostat-prototypePS_800mK_*", 1, "hist", "hist", 44621., 0.62, "test_git");
    
    
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






