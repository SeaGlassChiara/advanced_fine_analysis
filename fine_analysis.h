#ifndef ALCOR_UTILITYFINEANALYSIS_H
#define ALCOR_UTILITYFINEANALYSIS_H
//
//  --- --- --- ---
//
//  Fine tune for ALCOR readout utilities
//
//  Authors: Chiara Fraticelli, Nicola Rubini
//
//  --- --- --- ---
//
#include "analysis_utils.h"
using namespace analysis_utils;
//
//  Global Fine tune data Structure
std::map<global_identifier_t,std::tuple<double, double>>    kFineTuneParameters;        // Is detector unit oriented
std::map<TString,TH2F*>                                     kFineTuneRawHistograms;     // Is Run oriented
//
//  Global Fine tune variables
const TString   output_preprocess_fine_analysis_directory       = output_preprocess_directory                       + TString("/FineAnalysis/");
const TString   output_preprocess_fine_analysis_file_raw        = output_preprocess_fine_analysis_directory         + TString("/%s/FineTuneRaw.root");
const TString   output_preprocess_fine_analysis_file_rslt       = output_preprocess_fine_analysis_directory         + TString("/%s/FineTuneResults.root");
const TString   output_preprocess_fine_analysis_file_check      = output_preprocess_fine_analysis_directory         + TString("/%s/FineTuneResultsCheck.root");
const TString   output_preprocess_fine_analysis_graphics_dir    = output_preprocess_fine_analysis_directory         + TString("/%s/FitCheck/");
TF1*            fine_analysis_fit_function                      = new TF1("fine_analysis_fit_function", "[0]*(1./(exp((x-[1])/[2])+1))*(1./(exp(([3]-x)/[4])+1))", 20., 150.);
//
//  Global Fine tune functions
//  --- Calibration basic function
double  calculate_calibrated_phase      ( int kFineParameter, double kCalibrationMin, double kCalibrationMax );
double  calculate_calibrated_phase      ( int kFineParameter, global_identifier_t kFilters );
double  calculate_calibrated_phase      ( int kFineParameter, TString kRunTag, int kGlobalIndex ) { return calculate_calibrated_phase(kFineParameter,{kRunTag,kGlobalIndex}); };
template< typename TH1_Type = TH1F >
int     get_fit_bump                    ( TH1_Type* histo, double critical_value, int bin_step );
//  --- Getters
std::tuple<double,double>
        get_fine_calibration            ( global_identifier_t kFilters );
double  get_fine_calibration_min        ( global_identifier_t kFilters );
double  get_fine_calibration_max        ( global_identifier_t kFilters );
template< typename TH2_Type = TH2F >
TH2_Type*
        get_fine_tune_parameters        ( TString kRunTag );
//  --- Setters
void    set_fine_calibration            ( global_identifier_t kFilters, double kCalibrationMin, double kCalibrationMax );
void    set_fine_calibration            ( global_identifier_t kFilters );
void    set_fine_analysis_fit_function  ( TH1D* histo, int riseThreshold, int dropThreshold, int max_step, int min_step );
//  --- Fine Tune Raw utilities
template< typename TH2_Type = TH2F >
TH2_Type*
        build_fine_tune_raw_histogram   ( std::vector<TString> kInputFileNames, TString kRunTag, TString kOutputFileName, bool kRecalculate );
template< typename TH2_Type = TH2F >
TH2_Type*
        build_fine_tune_raw_histogram   ( TString kRunTag, TString kOutputFileName, bool kRecalculate );
template< typename TH2_Type = TH2F >
void    set_fine_tune_raw_histogram     ( global_identifier_t kFilters, TH2_Type* kFineTuneRawHistogram );
void    set_fine_tune_raw_histogram     ( global_identifier_t kFilters );
TH1D*   get_fine_tune_raw_histogram     ( global_identifier_t kFilters );
template< typename TH2_Type = TH2F >
void    set_fine_tune_raw_histogram     ( TString kRunTag, TH2_Type* kFineTuneRawHistogram );
void    set_fine_tune_raw_histogram     ( TString kRunTag );
template< typename TH2_Type = TH2F >
TH2_Type*
        get_fine_tune_raw_histogram     ( TString kRunTag );
//  --- Fine Tune Analysis utilities
template< typename TH2_Type = TH2F >
TH2_Type*
        run_fine_tune_analysis          ( TString kRunTag, TString kOutputFileName, TString kOutputGraphics, bool kRecalculate );
//
//  Calibration basic function
double
calculate_calibrated_phase
( int kFineParameter, double kCalibrationMin, double kCalibrationMax ) {
    double kHalfCut         = 0.5 * ( kCalibrationMin + kCalibrationMax);
    double kNormalisation   = ( kCalibrationMax - kCalibrationMin );
    double kCurrentPhase    = ( kFineParameter - kCalibrationMin ) / kNormalisation;
    return kFineParameter < kHalfCut ? kCurrentPhase : kCurrentPhase -= 1;
}
double
calculate_calibrated_phase
( int kFineParameter, global_identifier_t kFilters ) { return calculate_calibrated_phase(kFineParameter, get_fine_calibration_min(kFilters), get_fine_calibration_max(kFilters)); }
template< typename TH1_Type = TH1F >
int
get_fit_bump
( TH1_Type* histo, double critical_value, int bin_step ) {
    if ( bin_step > 0 ) {
        for ( int j=1; j<=histo->GetNbinsX()-bin_step; j+=bin_step ) {
            double y1 = histo->GetBinContent(j);
            double y2 = histo->GetBinContent(j+bin_step);
            if ( fabs(y1-y2) > critical_value ){
                return (int)(j+bin_step*0.5);
            }
        }
    } else if ( bin_step < 0 ) {
        for ( int j=histo->GetNbinsX(); j>(-bin_step); j+=bin_step ) {
            double y1 = histo->GetBinContent(j);
            double y2 = histo->GetBinContent(j+bin_step);
            if ( fabs(y1-y2) > critical_value ){
                return (int)(j+bin_step*0.5);
            }
        }
    }
    return -1;
}
//
//  Getters
double  get_fine_calibration_min        ( global_identifier_t kFilters ) { return get<0>(get_fine_calibration(kFilters));   };
double  get_fine_calibration_max        ( global_identifier_t kFilters ) { return get<1>(get_fine_calibration(kFilters));   };
std::tuple<double,double>
get_fine_calibration
 ( global_identifier_t kFilters ) {
    if ( kFineTuneParameters.find( kFilters ) == kFineTuneParameters.end() ) set_fine_calibration( kFilters );
    return kFineTuneParameters[kFilters];
}
//
//  Setters
void    set_fine_calibration            ( global_identifier_t kFilters, double kCalibrationMin, double kCalibrationMax ) { kFineTuneParameters[kFilters] = { kCalibrationMin, kCalibrationMax }; };
void
set_fine_calibration
( global_identifier_t kFilters ) {
    auto kCurrentFineAnalysis = get_fine_tune_parameters( get<0>(kFilters) );
    auto kCurrentFineMinimum = kCurrentFineAnalysis->GetBinContent(get<1>(kFilters),4);
    auto kCurrentFineMaximum = kCurrentFineAnalysis->GetBinContent(get<1>(kFilters),2);
    set_fine_calibration( kFilters, kCurrentFineMinimum, kCurrentFineMaximum );
}
void
set_fine_analysis_fit_function
 ( TH1D* histo, int riseThreshold =1, int dropThreshold =1, int max_step = -3, int min_step = 3) {
    double guess = histo->GetBinContent(70);
    dropThreshold = 0.25*guess;
    riseThreshold = 0.25*guess;
    double max = get_fit_bump(histo, dropThreshold, max_step);
    double min = get_fit_bump(histo, riseThreshold, min_step);
    double height   = histo->GetBinContent((int)(0.5*(max+min)));
    double minguess = histo->GetBinCenter(min);
    double maxguess = histo->GetBinCenter(max);
    // Set Parameter 0
    fine_analysis_fit_function->SetParLimits(0, height*0.7, height*1.3);
    fine_analysis_fit_function->SetParameter(0, height );
    // Set Parameter 1
    fine_analysis_fit_function->SetParLimits(1, maxguess*0.7, maxguess*1.3);
    fine_analysis_fit_function->SetParameter(1, maxguess);
    // Set Parameter 2
    fine_analysis_fit_function->SetParLimits(2, 0, 1.);
    fine_analysis_fit_function->SetParameter(2, 0.4);
    // Set Parameter 3
    fine_analysis_fit_function->SetParLimits(3, minguess*0.7, minguess*1.3);
    fine_analysis_fit_function->SetParameter(3, minguess);
    // Set Parameter 4
    fine_analysis_fit_function->SetParLimits(4, 0, 1.);
    fine_analysis_fit_function->SetParameter(4, 0.4);
}
//
//  Fine Tune Raw utilities
template< typename TH2_Type = TH2F >
TH2_Type*
build_fine_tune_raw_histogram
 ( std::vector<TString> kInputFileNames, TString kRunTag, TString kOutputFileName, bool kRecalculate ) {
    //!  Check the output has already been produced
    TFile*  kFileOut = new TFile( Form (output_preprocess_fine_analysis_file_raw,kRunTag.Data()));
    //!  (Re)Calculate if requested
    cout << Form("%s%s",output_preprocess_fine_analysis_directory.Data(), Form (output_preprocess_fine_analysis_file_raw,kRunTag.Data())) << endl;
    cout << "[INFO] Looking if fine tune histogram for tuning is cached" << endl;
    if ( kFileOut->IsOpen() && !kRecalculate )  { cout << "[INFO] Found! Using cached" << endl; return (TH2F*)(kFileOut->Get("hFine_All")); }
    else { cout << "[INFO] Not found! Re-generating" << endl; delete  kFileOut; }
    //! histogram with fine distribution
    TH2_Type* hFine_All = new TH2_Type("hFine_All", ";index;fine", kGlobalIndexRange, 0, kGlobalIndexRange, kFineRange, 0, kFineRange);
    //! Loop on filenames
    for ( auto kCurrentFileName : kInputFileNames ) {
        std::cout << "[INFO] Opening file: " << kCurrentFileName.Data() << std::endl;
        //! Load File
        auto kCurrentFile = TFile::Open(kCurrentFileName);
        if ( !kCurrentFile || !kCurrentFile->IsOpen() ) {
            std::cout << "[WARNING] Opening file: " << kCurrentFileName.Data() << " failed!" << std::endl;
            continue;
        }
        //! Load Tree
        data_t      kCurrentData;
        TTree*      kCurrentTree = (TTree *)kCurrentFile->Get("alcor");
        if ( !kCurrentTree ) {
            std::cout << "[WARNING] No \"alcor\" tree in file: " << kCurrentFileName.Data() << ". Loading failed!" << std::endl;
            continue;
        }
        //! Get Entries and loop to fill histogram
        auto nEvents = kCurrentTree->GetEntries();
        load_tree( kCurrentTree, kCurrentData );
        for (int iEv = 0; iEv < nEvents; ++iEv) {
            kCurrentTree->GetEntry(iEv);
            int iCurrentIndex = get_global_index( kCurrentData.fifo, kCurrentData.pixel, kCurrentData.column, kCurrentData.tdc );
            hFine_All->Fill( iCurrentIndex, kCurrentData.fine);
        }
        kCurrentFile->Close();
    }
    //! Save custom output
    if ( kOutputFileName.Length() != 0 ) {
        // system(Form("mkdir -p %s", kOutputFileName.Data())); //! TODO: Create folder (?)
        kFileOut = new TFile( kOutputFileName, "RECREATE" );
        hFine_All->Write();
        kFileOut->Close();
    } else {
        system(Form("mkdir -p %s/%s",output_preprocess_fine_analysis_directory.Data(),kRunTag.Data()));
        kFileOut = new TFile( Form (output_preprocess_fine_analysis_file_raw,kRunTag.Data()), "RECREATE" );
        hFine_All->Write();
        kFileOut->Close();
    }
    return hFine_All;
}
template< typename TH2_Type = TH2F >
TH2_Type*
build_fine_tune_raw_histogram
 ( TString kRunTag, TString kOutputFileName, bool kRecalculate ) {
    std::vector<TString> kInputFileNames;
    for ( Int_t iFile = 0; iFile < 24; iFile++ ) { kInputFileNames.push_back(Form(intput_rawdata_decoded_file,kRunTag.Data(),iFile)); }
    return build_fine_tune_raw_histogram<TH2_Type>( kInputFileNames, kRunTag, kOutputFileName, kRecalculate );
}
template< typename TH2_Type = TH2F >
void    set_fine_tune_raw_histogram     ( global_identifier_t kFilters, TH2_Type* kFineTuneRawHistogram ) { set_fine_tune_raw_histogram( get<0>(kFilters), kFineTuneRawHistogram ); }
void    set_fine_tune_raw_histogram     ( global_identifier_t kFilters ) { set_fine_tune_raw_histogram( get<0>(kFilters) ); }
TH1D*   get_fine_tune_raw_histogram     ( global_identifier_t kFilters ) { return get_fine_tune_raw_histogram( get<0>(kFilters) )->ProjectionY( Form("hFine_%i", get<1>(kFilters) ), get<1>(kFilters) + 1, get<1>(kFilters) + 1); }
template< typename TH2_Type = TH2F >
void    set_fine_tune_raw_histogram     ( TString kRunTag, TH2_Type* kFineTuneRawHistogram ) { kFineTuneRawHistograms[kRunTag] = kFineTuneRawHistogram; }
void    set_fine_tune_raw_histogram     ( TString kRunTag ) { kFineTuneRawHistograms[kRunTag] = build_fine_tune_raw_histogram( kRunTag, "", false ); };
template< typename TH2_Type = TH2F >
TH2_Type*
get_fine_tune_raw_histogram
 ( TString kRunTag ) {
    if (kFineTuneRawHistograms.find(kRunTag) == kFineTuneRawHistograms.end()) set_fine_tune_raw_histogram(kRunTag);
    return kFineTuneRawHistograms[kRunTag];
}
//
//  Fine Tune Analysis
template< typename TH2_Type = TH2F >
TH2_Type*
run_fine_tune_analysis
 ( TString kRunTag, TString kOutputFileName, TString kOutputGraphics, bool kRecalculate ) {
    //!  Check the output has already been produced
    TFile*  kFileOut = new TFile( Form (output_preprocess_fine_analysis_file_rslt,kRunTag.Data()) );
    //!  (Re)Calculate if requested
    cout << "[INFO] Looking if fine tune histogram for tuning is cached" << endl;
    if ( kFileOut->IsOpen() && !kRecalculate )  { cout << "[INFO] Found! Using cached" << endl; return (TH2F*)(kFileOut->Get("kFine_All_Tune_Params")); }
    else { cout << "[INFO] Not found! Re-generating" << endl; delete  kFileOut; }
    //! Run in Batch mode
    gROOT->SetBatch(true);
    //! Load Files
    auto hFine_All = get_fine_tune_raw_histogram( kRunTag );
    //! Create output
    TH2_Type*   kFine_All_Tune_Params   = new TH2_Type( "kFine_All_Tune_Params", "kFine_All_Tune_Params", kGlobalIndexRange, -0.5, kGlobalIndexRange-0.5, 5, 0., 5.);
    TH1F*       kMaximumDistribution    = new TH1F( "kMaximumDistribution", "kMaximumDistribution", 200, 0, 200 );
    TH1F*       kMinimumDistribution    = new TH1F( "kMinimumDistribution", "kMinimumDistribution", 200, 0, 200 );
    kFine_All_Tune_Params->GetYaxis()->SetBinLabel(1,"Normalisation");
    kFine_All_Tune_Params->GetYaxis()->SetBinLabel(2,"Maximum");
    kFine_All_Tune_Params->GetYaxis()->SetBinLabel(3,"Sigma Max");
    kFine_All_Tune_Params->GetYaxis()->SetBinLabel(4,"Minimum");
    kFine_All_Tune_Params->GetYaxis()->SetBinLabel(5,"Sigma Min");
    //! Loop on found infos
    for ( Int_t iIndex = 0; iIndex < kGlobalIndexRange; iIndex++ ) {
        global_identifier_t kCurrentFilter;
        get<0>(kCurrentFilter)  = kRunTag;
        get<1>(kCurrentFilter)  = iIndex;
        auto kCurrentInfos      = uGetInfos(iIndex);
        auto kCurrentFineHisto  = get_fine_tune_raw_histogram( kCurrentFilter );
        if ( kCurrentFineHisto->GetEntries() <= 100 ) {
            cout << "[WARNING] Skipping empty histogram " << iIndex << " chip:" << get<0>(kCurrentInfos) << " pixel:" <<  get<1>(kCurrentInfos) << " column:" <<  get<2>(kCurrentInfos) << " TDC:" << get<3>(kCurrentInfos)  << endl;
            continue;
        }
        set_fine_analysis_fit_function( kCurrentFineHisto );
        kCurrentFineHisto->Fit( fine_analysis_fit_function, "SQ", "" );
        kFine_All_Tune_Params->SetBinContent( iIndex+1, 1, fine_analysis_fit_function->GetParameter(0) );
        kFine_All_Tune_Params->SetBinContent( iIndex+1, 2, fine_analysis_fit_function->GetParameter(1) );
        kFine_All_Tune_Params->SetBinContent( iIndex+1, 3, fine_analysis_fit_function->GetParameter(2) );
        kFine_All_Tune_Params->SetBinContent( iIndex+1, 4, fine_analysis_fit_function->GetParameter(3) );
        kFine_All_Tune_Params->SetBinContent( iIndex+1, 5, fine_analysis_fit_function->GetParameter(4) );
        kMaximumDistribution->Fill(fine_analysis_fit_function->GetParameter(1));
        kMinimumDistribution->Fill(fine_analysis_fit_function->GetParameter(3));
        if ( !(kOutputGraphics.Length() == 0) ) {
            system(Form("mkdir -p %s",kOutputGraphics.Data()));
            TLatex* kLatex = new TLatex();
            auto kMinimum = fine_analysis_fit_function->GetParameter(3);
            auto kMaximum = fine_analysis_fit_function->GetParameter(1);
            TCanvas* c1 = new TCanvas("", "", 1000, 300);
            c1->Divide(3,1);
            c1->cd(1);
            kCurrentFineHisto->GetXaxis()->SetRangeUser(kMinimum -10, kMinimum +20);
            kCurrentFineHisto->DrawCopy();
            kLatex->DrawLatexNDC(0.6, 0.5, Form("Min: %.2f", kMinimum));
            c1->cd(2);
            kCurrentFineHisto->GetXaxis()->SetRangeUser(kMinimum -10, kMaximum +10);
            kCurrentFineHisto->DrawCopy();
            kLatex->DrawLatexNDC(0.4, 0.5, Form("Max: %.2f", kMaximum));
            kLatex->DrawLatexNDC(0.4, 0.45, Form("Min: %.2f", kMinimum));
            c1->cd(3);
            kCurrentFineHisto->GetXaxis()->SetRangeUser(kMaximum -20, kMaximum +10);
            kCurrentFineHisto->DrawCopy();
            kLatex->DrawLatexNDC(0.3, 0.5, Form("Max: %.2f", kMaximum));
            c1->SaveAs(Form("%s/%s.pdf",kOutputGraphics.Data(),kCurrentFineHisto->GetName()));
        } else {
            system(Form("mkdir -p %s", Form( output_preprocess_fine_analysis_graphics_dir.Data(),kRunTag.Data())));
            TLatex* kLatex = new TLatex();
            auto kMinimum = fine_analysis_fit_function->GetParameter(3);
            auto kMaximum = fine_analysis_fit_function->GetParameter(1);
            TCanvas* c1 = new TCanvas("", "", 1000, 300);
            c1->Divide(3,1);
            c1->cd(1);
            kCurrentFineHisto->GetXaxis()->SetRangeUser(kMinimum -10, kMinimum +20);
            kCurrentFineHisto->DrawCopy();
            kLatex->DrawLatexNDC(0.6, 0.5, Form("Min: %.2f", kMinimum));
            c1->cd(2);
            kCurrentFineHisto->GetXaxis()->SetRangeUser(kMinimum -10, kMaximum +10);
            kCurrentFineHisto->DrawCopy();
            kLatex->DrawLatexNDC(0.4, 0.5, Form("Max: %.2f", kMaximum));
            kLatex->DrawLatexNDC(0.4, 0.45, Form("Min: %.2f", kMinimum));
            c1->cd(3);
            kCurrentFineHisto->GetXaxis()->SetRangeUser(kMaximum -20, kMaximum +10);
            kCurrentFineHisto->DrawCopy();
            kLatex->DrawLatexNDC(0.3, 0.5, Form("Max: %.2f", kMaximum));
            c1->SaveAs(Form("%s/%s.pdf", Form( output_preprocess_fine_analysis_graphics_dir.Data(),kRunTag.Data()),kCurrentFineHisto->GetName()));
        }
    }
    //!  Save custom output
    if ( kOutputFileName.Length() != 0 ) {
        // system(Form("mkdir -p %s", kOutputFileName.Data())); //! TODO: Create folder (?)
        kFileOut = new TFile( kOutputFileName, "RECREATE" );
        kMaximumDistribution->Write();
        kMinimumDistribution->Write();
        hFine_All->Write();
        kFine_All_Tune_Params->Write();
        kFileOut->Close();
    } else {
        system(Form("mkdir -p %s/%s",output_preprocess_fine_analysis_directory.Data(),kRunTag.Data()));
        kFileOut = new TFile( Form (output_preprocess_fine_analysis_file_rslt,kRunTag.Data()), "RECREATE" );
        kMaximumDistribution->Write();
        kMinimumDistribution->Write();
        hFine_All->Write();
        kFine_All_Tune_Params->Write();
        kFileOut->Close();
    }
    //! Run in Batch mode
    gROOT->SetBatch(false);
    return kFine_All_Tune_Params;
}
template< typename TH2_Type = TH2F >
TH2_Type*
get_fine_tune_parameters ( TString kRunTag ) {
    //!  Check the output has already been produced
    TFile*  kFileOut = new TFile( Form (output_preprocess_fine_analysis_file_rslt,kRunTag.Data()) );
    //!  (Re)Calculate if requested
    if ( kFileOut->IsOpen() )  { return (TH2F*)(kFileOut->Get("kFine_All_Tune_Params")); }
    return run_fine_tune_analysis( kRunTag, "", "", false );
}
template< typename TH2_Type = TH2F >
TH2_Type*
check_normalisation_fine_tune
 ( std::vector<TString> kInputFileNames, TString kRunTag, TString kOutputFileName="", bool kRecalculate=false ) {
    //!  Check the output has already been produced
    TFile*  kFileOut = new TFile( Form (output_preprocess_fine_analysis_file_check,kRunTag.Data()));
    //!  (Re)Calculate if requested
    cout << Form("%s%s",output_preprocess_fine_analysis_directory.Data(), Form (output_preprocess_fine_analysis_file_raw,kRunTag.Data())) << endl;
    cout << "[INFO] Looking if fine tune histogram for tuning is cached" << endl;
    if ( kFileOut->IsOpen() && !kRecalculate )  { cout << "[INFO] Found! Using cached" << endl; return (TH2F*)(kFileOut->Get("hFine_All")); }
    else { cout << "[INFO] Not found! Re-generating" << endl; delete  kFileOut; }
    //! histogram with fine distribution
    TH2_Type* hFine_All_Tuned = new TH2_Type("hFine_All_Tuned", ";index;fine", kGlobalIndexRange, 0, kGlobalIndexRange, 189, -1, 2 );
    //! Loop on filenames
    for ( auto kCurrentFileName : kInputFileNames ) {
        std::cout << "[INFO] Opening file: " << kCurrentFileName.Data() << std::endl;
        //! Load File
        auto kCurrentFile = TFile::Open(kCurrentFileName);
        if ( !kCurrentFile || !kCurrentFile->IsOpen() ) {
            std::cout << "[WARNING] Opening file: " << kCurrentFileName.Data() << " failed!" << std::endl;
            continue;
        }
        //! Load Tree
        data_t      kCurrentData;
        TTree*      kCurrentTree = (TTree *)kCurrentFile->Get("alcor");
        if ( !kCurrentTree ) {
            std::cout << "[WARNING] No \"alcor\" tree in file: " << kCurrentFileName.Data() << ". Loading failed!" << std::endl;
            continue;
        }
        //! Get Entries and loop to fill histogram
        auto nEvents = kCurrentTree->GetEntries();
        load_tree( kCurrentTree, kCurrentData );
        for (int iEv = 0; iEv < nEvents; ++iEv) {
            kCurrentTree->GetEntry(iEv);
            int iCurrentIndex = get_global_index( kCurrentData.fifo, kCurrentData.pixel, kCurrentData.column, kCurrentData.tdc );
            hFine_All_Tuned->Fill( iCurrentIndex, calculate_calibrated_phase(kCurrentData.fine,kRunTag,iCurrentIndex) );
        }
        kCurrentFile->Close();
    }
    //! Save custom output
    if ( kOutputFileName.Length() != 0 ) {
        // system(Form("mkdir -p %s", kOutputFileName.Data())); //! TODO: Create folder (?)
        kFileOut = new TFile( kOutputFileName, "RECREATE" );
        hFine_All_Tuned->Write();
        kFileOut->Close();
    } else {
        system(Form("mkdir -p %s/%s",output_preprocess_fine_analysis_directory.Data(),kRunTag.Data()));
        kFileOut = new TFile( Form (output_preprocess_fine_analysis_file_check,kRunTag.Data()), "RECREATE" );
        hFine_All_Tuned->Write();
        kFileOut->Close();
    }
    return hFine_All_Tuned;
    
    
    /*
    //! Create output
    TH2_Type* hFine_All_Tuned = new TH2_Type("hFine_All_Tuned", "hFine_All_Tuned", kGlobalIndexRange, -0.5, kGlobalIndexRange-0.5, 189, -1, 2 );
    //! Loop on filenames
    for ( auto kCurrentFileName : kInputFileNames ) {
        std::cout << "[INFO] Opening file: " << kCurrentFileName.Data() << std::endl;
        //! Load FileF
        auto kCurrentFile = TFile::Open(kCurrentFileName);
        if ( !kCurrentFile || !kCurrentFile->IsOpen() ) {
            std::cout << "[WARNING] Opening file: " << kCurrentFileName.Data() << " failed!" << std::endl;
            continue;
        }
        //! Load Tree
        data_t    kCurrentData;
        TTree*      kCurrentTree = (TTree *)kCurrentFile->Get("alcor");
        if ( !kCurrentTree ) {
            std::cout << "[WARNING] No \"alcor\" tree in file: " << kCurrentFileName.Data() << ". Loading failed!" << std::endl;
            continue;
        }
        //! Get Entries and loop to fill histogram
        auto nEvents = kCurrentTree->GetEntries();
        LoadTree( kCurrentTree, kCurrentData );
        //!
        //! kFine\_All\_Tune\_Params con i dati calibrati
        for (int iEv = 0; iEv < nEvents; ++iEv) {
            kCurrentTree->GetEntry(iEv);
            int     iCurrentIndex   = get_global_index( kCurrentData.fifo, kCurrentData.pixel, kCurrentData.column, kCurrentData.tdc );
            if ( iCurrentIndex < 0 ) {
                std::cout << "[WARNING] Invalid global index " << iCurrentIndex << " found! Skipping..." << std::endl;
                continue;
            }
            hFine_All_Tuned->Fill( iCurrentIndex, calculate_calibrated_phase(kCurrentData.fine, {kRunTag,iCurrentIndex}) );
        }
        kCurrentFile->Close();
    }
    //! Save output
    if ( !(kOutputFileName.Length() == 0) ) {
        TFile*  kFileOut = new TFile( kOutputFileName, "RECREATE" );
        kFine_All_Tune_Params->Write();
        hFine_All_Tuned->Write();
        kFileOut->Close();
    }
    //! Run in Batch mode
    gROOT->SetBatch(false);
    return hFine_All_Tuned;
     */
}
template< typename TH2_Type = TH2F >
TH2_Type*
check_normalisation_fine_tune
 ( TString kRunTag, TString kOutputFileName="", bool kRecalculate=false ) {
    std::vector<TString> kInputFileNames;
    for ( Int_t iFile = 0; iFile < 24; iFile++ ) { kInputFileNames.push_back(Form(intput_rawdata_decoded_file,kRunTag.Data(),iFile)); }
    return check_normalisation_fine_tune<TH2_Type>( kInputFileNames, kRunTag, kOutputFileName, kRecalculate );
}
#endif
