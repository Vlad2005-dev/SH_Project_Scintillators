// extract_histograms.C
// Saves every histogram stored by pmt_events in each *_pulse_events.root file of a folder,
// one folder per histogram, so the same plot at every HV sits side by side:
//   <folder>/histograms/<histogram>/log_y/<run>.png
//   <folder>/histograms/<histogram>/linear_y/<run>.png
// e.g. histograms/hCharge/log_y/2000V_2000V_20min_6_5_int1_wave_0.png
//
// Run (from the code folder):
//     root -l extract_histograms.C                                   // uses DATA_DIR below
//     root -l 'extract_histograms.C("D:/results/calibration_data")'  // another folder

#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TList.h"

// ---- settings ----
const char *DATA_DIR = "../results/calibration_data";   // folder with the *_pulse_events.root files

void extract_histograms(const char *dir = DATA_DIR)
{
   gStyle->SetOptStat(1110);   // entries, mean, std dev
   std::string d = dir;
   std::string outDir = d + "/histograms/";
   const std::string ending = "_pulse_events.root";

   // ---- find the files ----
   std::vector<std::string> names;
   TSystemDirectory sd("d", dir);
   TList *files = sd.GetListOfFiles();
   if (!files) { printf("Cannot read folder %s\n", dir); return; }
   for (TObject *o : *files) {
      std::string n = o->GetName();
      if (n.size() > ending.size() && n.substr(n.size() - ending.size()) == ending) names.push_back(n);
   }
   std::sort(names.begin(), names.end());
   if (names.empty()) { printf("No *%s files in %s\n", ending.c_str(), dir); return; }

   TCanvas *c = new TCanvas("cHist", "histogram", 900, 600);
   int nSaved = 0;
   for (const std::string &n : names) {
      TFile *f = TFile::Open((d + "/" + n).c_str());
      if (!f || f->IsZombie()) { printf("Cannot open %s\n", n.c_str()); continue; }
      std::string stem = n.substr(0, n.size() - ending.size());
      printf("%s:", stem.c_str());

      // every histogram in the file (the pulses tree is skipped)
      for (TObject *k : *f->GetListOfKeys()) {
         TH1 *h = dynamic_cast<TH1 *>(((TKey *)k)->ReadObj());
         if (!h) continue;
         std::string hDir = outDir + h->GetName() + "/";
         gSystem->mkdir((hDir + "log_y").c_str(), true);
         gSystem->mkdir((hDir + "linear_y").c_str(), true);
         h->SetTitle(Form("%s  (%s)", h->GetTitle(), stem.c_str()));
         c->cd();
         h->Draw("hist");
         c->SetLogy(1); c->SaveAs((hDir + "log_y/" + stem + ".png").c_str());
         c->SetLogy(0); c->SaveAs((hDir + "linear_y/" + stem + ".png").c_str());
         printf(" %s", h->GetName());
         nSaved++;
      }
      printf("\n");
      f->Close();
   }
   printf("%d histograms from %zu files, each with log and linear y -> %s\n", nSaved, names.size(), outDir.c_str());
}
