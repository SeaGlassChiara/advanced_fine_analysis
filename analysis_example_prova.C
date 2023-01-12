/**
   root [0] .L analysis_example.C
   root [1] analysis_example("../data/20221015-234954/decoded")
**/

#include "analysis_utils.h"
#include "fine_analysis.h"

int frame_size = 1024; // [clock cycles]
const double reference_delay = 10.; // [ns]

namespace au = analysis_utils;

void
analysis_example_prova(std::string dirname, std::vector<std::string> filenames = au::all_filenames)
{

  /** 
   ** HISTOGRAMS
   **/
      
  auto hN = new TH2F("hN", ";N_{TIME0} fired pixels;N_{TIME1} fired pixels;", 33, 0., 33., 33, 0., 33.);
  TH1 *hDelta[6] = {nullptr};
  TH2 *hMap[6] = {nullptr};
  
  TH2F *t_vs_fine = new TH2F("t_vs_fine", "t_vs_fine", (int)(0.25*kFineRange), 0, (int)(0.25*kFineRange), 2 * 60 * 25, -0 * au::coarse_to_ns, 7 * au::coarse_to_ns);
  TH2F *t_vs_phase = new TH2F("t_vs_phase", "t_vs_phase", 100, -0.6, 0.6, 2 * 60 * 25, -0 * au::coarse_to_ns, 7 * au::coarse_to_ns);
  TH2F *traw_vs_fine = new TH2F("traw_vs_fine", "traw_vs_fine", (int)(0.25*kFineRange), 0, (int)(0.25*kFineRange), 2 * 60 * 25, -0 * au::coarse_to_ns, 7 * au::coarse_to_ns);
  TH2F *traw_vs_phase = new TH2F("traw_vs_phase", "traw_vs_phase", 100, -0.6, 0.6, 2 * 60 * 25, -0 * au::coarse_to_ns, 7 * au::coarse_to_ns);
  
  TH1F *hDelta_t_chip4_chip5 = new TH1F("hDelta_t_chip4_chip5", "hDelta_t_chip4_chip5", 100, -15, 15); 
  
  
  for (int ichip = 0; ichip < 6; ++ichip) {
    hMap[ichip] = new TH2F(Form("hMap_%d", ichip), ";pixel row;pixel column", 8, 0., 8., 4, 0., 4.);
    hDelta[ichip] = new TH1F(Form("hDelta_%d", ichip), "hit - reference time (ns)", 2 * frame_size, -frame_size * au::coarse_to_ns, frame_size * au::coarse_to_ns);
  }
  
  /** 
   ** POPULATE FRAMED DATA
   **/
  
  au::framed_data_t framed_data;
  while (au::populate_framed_data(framed_data, dirname, filenames, frame_size)) {
  
  /** 
   ** POST-PROCESSING ANALYSIS 
   **/
      
  std::cout << " --- post processing " << std::endl;

  /** 
   ** LOOP OVER DATA
   **/
      
  /** loop over spills **/
  for (auto &spill_data : framed_data) {
    auto spill = spill_data.first;
    auto &frames = spill_data.second;

    /** loop over frames **/
    for (auto &frame_data : frames) {
      auto frame = frame_data.first;
      auto &chips = frame_data.second;

      /** 
       ** COMPUTE REFERENCE TIME 
       **/
      
      double average_time[6] = {0.};
      double average_time_raw[6] = {0.};

      /** loop over chips **/
      for (auto &chip_data : chips) {
        auto chip = chip_data.first;
        auto &channels = chip_data.second;
        
        /** loop over channels **/
        for (auto &channel_data : channels) {
          auto channel = channel_data.first;
          auto &hits = channel_data.second;

          /** compute time of first hit **/
          auto &hit = hits[0];
          double time_raw = hit.coarse * au::coarse_to_ns + hit.rollover * au::rollover_to_ns; // [ns]
          average_time_raw[chip] += time_raw;
          
          double time = hit.coarse * au::coarse_to_ns + hit.rollover * au::rollover_to_ns + calculate_calibrated_phase( hit.fine, "20221015-234954", au::get_global_index( hit.fifo, hit.pixel, hit.column, hit.tdc ) ) * au::coarse_to_ns; // [ns] 
          average_time[chip] += time;
          
        } /** end of loop over channels **/

        if (channels.size() > 0)
          average_time[chip] /= channels.size();
          average_time_raw[chip] /= channels.size();
        
      } /** end of loop over chips **/
      
      /** request 32 fired pixels on both timing scintillators **/
      hN->Fill(chips[4].size(), chips[5].size());
      if (chips[4].size() != 32 || chips[5].size() != 32) continue;

      /** request compatible timing between timing scintillators **/      
      auto delta = average_time[4] - average_time[5];
      //riempi histo con questo delta
      hDelta_t_chip4_chip5->Fill(delta);
      
      if (fabs(delta) > 1.) continue;

      /** compute reference time as average of scintillators 
          keeping in mind they are downstream and are delayed 
          with respect to cherenkov photons **/
      auto reference_raw = 0.5 * (average_time_raw[4] + average_time_raw[5]) - reference_delay; // [ns]
      auto reference = 0.5 * (average_time[4] + average_time[5]) - reference_delay; // [ns]

      /** 
       ** CHERENKOV PHOTON TIMING 
       **/
      
      /** loop over chips **/
      for (auto &chip_data : chips) {
        auto chip = chip_data.first;
        auto &channels = chip_data.second;
        
      //  if (chip != 4) continue;

        /** loop over channels **/
        for (auto &channel_data : channels) {
          auto channel = channel_data.first;
          auto &hits = channel_data.second;
          
         if (channel != 0) continue;

          /** time of the first hit **/
          auto &hit = hits[0];
          
         if ( hit.tdc != 0 ) continue;
          
          auto index = au::get_index(hit.pixel, hit.column);
          double time_raw = hit.coarse * au::coarse_to_ns + hit.rollover * au::rollover_to_ns; // [ns]
          double time = hit.coarse * au::coarse_to_ns + hit.rollover * au::rollover_to_ns + calculate_calibrated_phase( hit.fine, "20221015-234954", au::get_global_index( hit.fifo, hit.pixel, hit.column, hit.tdc ) ) * au::coarse_to_ns; // [ns]
          double delta_raw = time_raw - reference;
          double delta = time - reference; 
          
          hDelta[chip]->Fill(delta);
          
           t_vs_fine->Fill(hit.fine, delta);
  	  t_vs_phase->Fill(calculate_calibrated_phase( hit.fine, "20221015-234954", au::get_global_index( hit.fifo, hit.pixel, hit.column, hit.tdc )), delta );
	  traw_vs_fine->Fill(hit.fine, delta_raw);
  	  traw_vs_phase->Fill(calculate_calibrated_phase( hit.fine, "20221015-234954", au::get_global_index( hit.fifo, hit.pixel, hit.column, hit.tdc ) ), delta_raw);

          /** request hit time to be compatible with reference timing from scintillators **/
          if (std::fabs(delta) < 10.)
            hMap[chip]->Fill(index.first, index.second);
          
        } /** end of loop over channels **/

      } /** end of loop over chips **/
            
    } /** end of loop over frames **/
    
  } /** end of loop over spills **/

	// break;
  
  }
  /** 
   ** WRITE OUTPUT TO FILE
   **/
      
  auto fout = TFile::Open("analysis_example_prova.root", "RECREATE");
  hN->Write();
  
    t_vs_fine->Write();
    t_vs_phase->Write();
    traw_vs_fine->Write();
    traw_vs_phase->Write();
    hDelta_t_chip4_chip5->Write();
    
  for (int ichip = 0; ichip < 6; ++ichip) {
    hMap[ichip]->Write();
    hDelta[ichip]->Write();
    
    
  }
  
  fout->Close();
  
  return;

}
