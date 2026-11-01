// extract_gain.C
// Takes the gain histogram (hGain) out of every *_pulse_events.root file in a folder and saves,
// in <folder>/gain_plots/:
//   log_y/, linear_y/      every plot twice, with a log and a linear y axis:
//      gain_<file>.png        one plot per file
//      gain_overlay.png       all files on one plot (same bins, normalised to area 1)
//      gain_vs_hv.png         mean and median gain against HV, with a fit G = a * HV^b
//   gain_vs_hv.csv         the numbers behind that plot
// HV is read from the start of the file name (e.g. 2000V_..._pulse_events.root -> 2000 V).
//
// Run (from the code folder):
//     root -l extract_gain.C                                   // uses DATA_DIR below
//     root -l 'extract_gain.C("D:/results/calibration_data")'  // another folder

#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TList.h"

// ---- settings ----
const char  *DATA_DIR   = "../results/calibration_data";  // folder with the *_pulse_events.root files
const double GAIN_MAX   = 4e7;     // overlay x range: 0 to this many electrons
const int    N_BINS     = 200;     // overlay bins
const double E_CHARGE_PC = 1.602176634e-7;   // electron charge in pC

// Save a canvas twice: <outDir>log_y/<name> and <outDir>linear_y/<name>
static void saveLogAndLinear(TCanvas *c, const std::string &outDir, const std::string &name)
{
   c->SetLogy(1); c->SaveAs((outDir + "log_y/" + name).c_str());
   c->SetLogy(0); c->SaveAs((outDir + "linear_y/" + name).c_str());
}

void extract_gain(const char *dir = DATA_DIR)
{
   gStyle->SetOptStat(1110);   // entries, mean, std dev
   std::string d = dir;
   std::string outDir = d + "/gain_plots/";
   gSystem->mkdir((outDir + "log_y").c_str(), true);
   gSystem->mkdir((outDir + "linear_y").c_str(), true);

   // ---- find the files, sorted by HV ----
   struct Run { double hv; std::string name; };
   std::vector<Run> runs;
   TSystemDirectory sd("d", dir);
   TList *files = sd.GetListOfFiles();
   if (!files) { printf("Cannot read folder %s\n", dir); return; }
   for (TObject *o : *files) {
      std::string n = o->GetName();
      if (n.size() < 18 || n.substr(n.size() - 18) != "_pulse_events.root") continue;
      runs.push_back({atof(n.c_str()), n});   // atof stops at the 'V': "2000V_..." -> 2000
   }
   std::sort(runs.begin(), runs.end(), [](const Run &a, const Run &b) { return a.hv < b.hv; });
   if (runs.empty()) { printf("No *_pulse_events.root files in %s\n", dir); return; }

   TCanvas *c1 = new TCanvas("cGain", "gain", 900, 600);
   TCanvas *c2 = new TCanvas("cOverlay", "gain overlay", 900, 600);
   TLegend *leg = new TLegend(0.70, 0.45, 0.88, 0.88);
   std::vector<double> hv, mean, median;
   FILE *csv = fopen((outDir + "gain_vs_hv.csv").c_str(), "w");
   fprintf(csv, "hv_V,pulses,mean_gain,median_gain\n");
   printf("\n%-8s %8s %12s %12s   file\n", "HV (V)", "pulses", "mean gain", "median gain");

   for (size_t i = 0; i < runs.size(); i++) {
      std::string path = d + "/" + runs[i].name;
      TFile *f = TFile::Open(path.c_str());
      if (!f || f->IsZombie()) { printf("Cannot open %s\n", path.c_str()); continue; }
      TH1D *h = (TH1D *)f->Get("hGain");
      TTree *t = (TTree *)f->Get("pulses");
      if (!h || !t) { printf("%s: no hGain or pulses tree\n", runs[i].name.c_str()); f->Close(); continue; }
      h->SetDirectory(nullptr);   // keep the histogram after the file is closed

      // per-file plot, as saved by pmt_events
      std::string stem = runs[i].name.substr(0, runs[i].name.size() - 18);
      c1->cd();
      h->SetTitle(Form("Gain, %g V;gain (electrons);pulses", runs[i].hv));
      h->Draw("hist");
      saveLogAndLinear(c1, outDir, "gain_" + stem + ".png");

      // gain of every pulse from the tree: median, and the overlay with common bins
      std::vector<double> g;
      double q; t->SetBranchAddress("charge_pC", &q);
      TH1D *ho = new TH1D(Form("hO%zu", i), "", N_BINS, 0, GAIN_MAX);
      ho->SetDirectory(nullptr);
      for (Long64_t k = 0; k < t->GetEntries(); k++) {
         t->GetEntry(k);
         g.push_back(q / E_CHARGE_PC);
         ho->Fill(q / E_CHARGE_PC);
      }
      std::sort(g.begin(), g.end());
      double med = g.empty() ? 0 : g[g.size() / 2];
      double avg = h->GetMean();

      ho->Scale(ho->Integral() > 0 ? 1.0 / ho->Integral() : 1.0);
      ho->SetLineColor(gStyle->GetColorPalette(i * (gStyle->GetNumberOfColors() - 1) / std::max<size_t>(1, runs.size() - 1)));
      ho->SetLineWidth(2);
      ho->SetStats(false);
      ho->SetTitle("Gain at each HV (normalised);gain (electrons);fraction of pulses");
      c2->cd();
      ho->Draw(i == 0 ? "hist" : "hist same");
      leg->AddEntry(ho, Form("%g V", runs[i].hv), "l");

      hv.push_back(runs[i].hv); mean.push_back(avg); median.push_back(med);
      fprintf(csv, "%g,%lld,%.4g,%.4g\n", runs[i].hv, (long long)t->GetEntries(), avg, med);
      printf("%-8g %8lld %12.3e %12.3e   %s\n", runs[i].hv, (long long)t->GetEntries(), avg, med, runs[i].name.c_str());
      f->Close();
   }
   fclose(csv);
   c2->cd(); leg->Draw();
   saveLogAndLinear(c2, outDir, "gain_overlay.png");

   // ---- gain against HV, fit G = a * HV^b ----
   TCanvas *c3 = new TCanvas("cHV", "gain vs HV", 900, 600);
   TGraph *gMean = new TGraph(hv.size(), hv.data(), mean.data());
   TGraph *gMed  = new TGraph(hv.size(), hv.data(), median.data());
   gMean->SetTitle("Gain against HV;HV (V);gain (electrons)");
   gMean->SetMarkerStyle(20); gMean->SetMarkerColor(kBlue + 1); gMean->SetLineColor(kBlue + 1);
   gMed->SetMarkerStyle(24);  gMed->SetMarkerColor(kRed + 1);   gMed->SetLineColor(kRed + 1);
   gMean->GetXaxis()->SetLimits(hv.front() - 100, hv.back() + 100);
   gMean->Draw("AP");
   gMed->Draw("P same");
   TF1 *fit = new TF1("fitHV", "[0]*pow(x,[1])", hv.front() * 0.95, hv.back() * 1.05);
   fit->SetParameters(mean.front() / pow(hv.front(), 7), 7);
   fit->SetLineColor(kBlue + 1);
   gMean->Fit(fit, "Q", "", hv.front() * 0.95, hv.back() * 1.05);
   TLegend *l3 = new TLegend(0.14, 0.72, 0.55, 0.88);
   l3->AddEntry(gMean, "mean gain", "p");
   l3->AddEntry(gMed, "median gain", "p");
   l3->AddEntry(fit, Form("fit to mean: G #propto HV^{%.2f}", fit->GetParameter(1)), "l");
   l3->Draw();
   c3->SetLogy(1); c3->SaveAs((outDir + "log_y/gain_vs_hv.png").c_str());
   gMean->SetMinimum(0);   // linear axis from 0 (a log axis can't start at 0, so set only here)
   c3->SetLogy(0); c3->SaveAs((outDir + "linear_y/gain_vs_hv.png").c_str());

   printf("\nfit to mean gain: G = %.3g * HV^%.2f  (+/- %.2f)\n", fit->GetParameter(0),
          fit->GetParameter(1), fit->GetParError(1));
   printf("plots and csv -> %s\n", outDir.c_str());
}
