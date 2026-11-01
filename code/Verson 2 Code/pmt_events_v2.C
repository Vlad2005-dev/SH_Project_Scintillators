// pmt_events_v2.C
// ROOT macro version of pmt_events_v2.py
//
// Usage:
//     root -l pmt_events_v2.C                                       // uses DATA_FILE below
//     root -l 'pmt_events_v2.C("wave_0__16_10_2026.dat")'           // any file name
//     root -l 'pmt_events_v2.C("/data/run5/wave_0.dat", 50000)'     // full path, first 50000 events
//     root -l 'pmt_events_v2.C("D:/data/23_10_2026/wave_0.dat", -1, "D:/tests")'   // outputs elsewhere
//     root -l pmt_events_v2.C+                             // compiled (faster) version
//
// Outputs go in a folder named after the folder holding the data, created next to this macro,
// with one sub-folder per data file:
//     data   D:/data/23_10_2026/wave_0.dat
//     macro  C:/pmt/pmt_events_v2.C
//     ->     C:/pmt/23_10_2026_results/wave_0/23_10_2026_wave_0_pulse_events.root  (and the plots)
// The third argument puts the <data folder>_results/ folder somewhere else instead.
// Works from any folder ROOT is started in, with / or \ paths.
//
// File layout: [24-byte header][1024 float32 samples] per event, little-endian.
// Samples are 12-bit ADC counts, converted to mV with mV = counts / 4096 * 1000.
// 1 time bin = 200 ps, so an event is 1024 bins = 204.8 ns.
//
// Per event:
//   baseline    = mean of bins BASE_LO..BASE_HI-1
//   noise rms   = standard deviation of those bins
//   time gate    = GATE_CENTER_NS +/- GATE_WIDTH_NS/2 (bins GATE_LO..GATE_HI-1)
//   pulse height = (baseline - lowest point) in the time gate
//   pulse if height > THRESH_SIGMA * rms; only these go into the tree, csv, waveforms and
//                height/charge/gain/peak-time plots.
//   NOT limited to the time gate: the noise rms and baseline histograms/stats (all events), and the
//                empty events and the "all pulses" peak-time plot, which look for the lowest point
//                in the whole record (0-204.8 ns).
//   integral    = sum of (baseline - V) * dt over the samples on the pulse side of the baseline
//                 around the peak: from where the signal crosses the baseline before the peak
//                 (walking back from the peak) to where it crosses it again after the peak
//   charge      = integral / R_OHM.  With V in mV and t in ns, (mV*ns)/ohm = pC.
//
// Outputs:
//   pulse_events.root   TTree "pulses" + the histograms
//   pulse_events.csv    one row per pulse event (only if WRITE_CSV = true; same content as the TTree)
//   pulse_waveforms     first N_SHOW pulse waveforms
//   pulse_stats         rms / baseline / height / charge / gain histograms
//   pulse_times         peak time of all pulses (whole record, time gate marked)
//                       and of the pulses in the time gate
//   each plot is saved as .png (if SAVE_PNG = true) and/or .pdf (if SAVE_PDF = true)

#include <cstdio>
#include <fstream>
#include <cstdint>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

// ---------------- settings you can change ----------------
const char  *DATA_FILE   = "/Users/vladimirsheydorov/SH_Project/raw_data/wave_0__23_10_2026.dat";  // default input (full or relative path) wave_0__23_10_2026.dat wave_0__16_10_2026.dat wave_0.dat
const char  *OUT_BASE    = "";      // where <data folder>_results/ is made; "" = the folder this macro is in
const bool   PREFIX_OUT  = true;    // true: output files start with <data folder>_<input name>_
const bool   WRITE_CSV   = false;   // true: also write the pulse table as a .csv (same content as the TTree)
const bool   SAVE_PNG    = false;    // true: save every plot as a .png
const bool   SAVE_PDF    = false;    // true: save every plot as a .pdf (both on = both saved)
const double PDF_LINE_SCALE = 1.5;  // line width multiplier in .pdf files (ROOT default 3; 1 = same as .png)
const double PNG_SCALE   = 2.0;     // .png resolution multiplier (1 = canvas size in pixels; 2 = sharp; 3 = print)
const int    N_SHOW      = 12;      // how many pulse waveforms to draw
const int    N_SAMP      = 1024;    // bins per event
const double BIN_NS      = 0.2;     // 200 ps per bin
const double POLARITY    = -1.0;    // -1 = negative pulses
const int    BASE_LO     = 10;      // baseline region: bins 10..209
const int    BASE_HI     = 210;
const double THRESH_SIGMA = 4.0;
const double GATE_CENTER_NS = 100.0;  // look for the pulse in a time gate around this time ...
const double GATE_WIDTH_NS  = 30.0;   // ... of this total width (center +/- width/2)
const double ADC_TO_MV   = 1000.0 / 4096.0;
const double R_OHM       = 60.0;    // input impedance used for Q = integral(V dt) / R
const double E_CHARGE    = 1.602176634e-19;   // electron charge in C, for gain = Q / e
// ----------------------------------------------------------

// time gate in bins (from GATE_CENTER_NS and GATE_WIDTH_NS): GATE_LO..GATE_HI-1, both edges included
const int GATE_LO = (int)std::lround((GATE_CENTER_NS - 0.5 * GATE_WIDTH_NS) / BIN_NS);
const int GATE_HI = (int)std::lround((GATE_CENTER_NS + 0.5 * GATE_WIDTH_NS) / BIN_NS) + 1;

// One record in the file: 6 x uint32 header, then 1024 floats = 4120 bytes
struct Record {
   uint32_t hdr[6];
   float    x[1024];
};

// Largest excursion from the baseline (pulse flipped up) in bins lo..hi-1, and where it is
static void findPeak(const double *v, double base, int lo, int hi, double &height, int &peak)
{
   height = -1e30;
   peak = lo;
   for (int i = lo; i < hi; i++) {
      double s = (v[i] - base) * POLARITY;   // flipped so the pulse points up
      if (s > height) { height = s; peak = i; }
   }
}

// Baseline, rms, pulse height, peak bin, integration window and charge for one waveform
// (voltages in mV), looking for the pulse in bins lo..hi-1 (default: the time gate).
// The charge is returned as a positive number in pC (pulse flipped up).
static void analyse(const double *v, double &base, double &rms, double &height, int &peak,
                    int &iStart, int &iEnd, double &charge, int lo = GATE_LO, int hi = GATE_HI)
{
   const int nb = BASE_HI - BASE_LO;
   double sum = 0;
   for (int i = BASE_LO; i < BASE_HI; i++) sum += v[i];
   base = sum / nb;

   double ss = 0;
   for (int i = BASE_LO; i < BASE_HI; i++) ss += (v[i] - base) * (v[i] - base);
   rms = std::sqrt(ss / nb);

   findPeak(v, base, lo, hi, height, peak);

   // integration window: start = first bin after the signal crosses the baseline before the peak,
   //                     end   = last bin before the signal crosses the baseline after the peak
   iStart = peak;
   while (iStart > 0 && (v[iStart - 1] - base) * POLARITY > 0) iStart--;
   iEnd = peak;
   while (iEnd < N_SAMP - 1 && (v[iEnd + 1] - base) * POLARITY > 0) iEnd++;

   double area = 0;                          // mV * ns
   for (int i = iStart; i <= iEnd; i++) area += (v[i] - base) * POLARITY * BIN_NS;
   charge = area / R_OHM;                    // mV*ns/ohm = pC
}

// Histogram range covering the central part of the data (0.1%..99.9%), with padding
static void zoomRange(std::vector<double> v, double &lo, double &hi)
{
   std::sort(v.begin(), v.end());
   auto pct = [&](double p) {
      double pos = p / 100.0 * (v.size() - 1);
      size_t i = (size_t)pos;
      double fr = pos - i;
      return (i + 1 < v.size()) ? v[i] * (1 - fr) + v[i + 1] * fr : v[i];
   };
   lo = pct(0.1);
   hi = pct(99.9);
   if (hi <= lo) {                       // (almost) constant data
      lo = v.front(); hi = v.back();
      if (hi <= lo) { lo -= 0.5; hi += 0.5; }
   }
   double w = hi - lo;
   lo -= 0.1 * w;
   hi += 0.1 * w;
}

static TH1D *makeHist(const char *name, const char *title, const std::vector<double> &v, int nbins)
{
   double lo, hi;
   zoomRange(v, lo, hi);
   TH1D *h = new TH1D(name, title, nbins, lo, hi);
   h->SetDirectory(nullptr);             // keep it alive independently of the output file
   for (double x : v) h->Fill(x);
   return h;
}

// Save a canvas as <name>.png and/or <name>.pdf, depending on SAVE_PNG / SAVE_PDF
static void savePlot(TCanvas *c, const std::string &name)
{
   if (SAVE_PNG) c->SaveAs((name + ".png").c_str());
   if (SAVE_PDF) c->SaveAs((name + ".pdf").c_str());
}

static void printStats(const char *name, std::vector<double> v)
{
   double mean = 0;
   for (double x : v) mean += x;
   mean /= v.size();
   double ss = 0;
   for (double x : v) ss += (x - mean) * (x - mean);
   double sd = std::sqrt(ss / v.size());
   std::sort(v.begin(), v.end());
   double med = v[v.size() / 2];
   printf("%-15s%12.3f%10.3f%12.3f%12.3f%12.3f\n", name, mean, sd, med, v.front(), v.back());
}

// Folders of a path, with "." and ".." resolved; relative paths start from the working directory
static std::vector<std::string> pathParts(const char *file)
{
   std::string p = gSystem->IsAbsoluteFileName(file) ? file
                 : std::string(gSystem->WorkingDirectory()) + "/" + file;
   std::vector<std::string> parts;
   size_t i = 0;
   while (i <= p.size()) {
      size_t j = p.find_first_of("/\\", i);
      if (j == std::string::npos) j = p.size();
      std::string s = p.substr(i, j - i);
      if (s == "..") { if (!parts.empty()) parts.pop_back(); }
      else if (!s.empty() && s != ".") parts.push_back(s);
      i = j + 1;
   }
   return parts;
}

void pmt_events_v2(const char *dataFile = DATA_FILE, Long64_t maxEvents = -1, const char *outBase = OUT_BASE)
{
   // output names: <macro folder>/<data folder>_results/<input name>/<data folder>_<input name>_pulse_events.root etc.
   std::vector<std::string> parts = pathParts(dataFile);
   std::string stem    = parts.back();
   if (stem.find_last_of('.') != std::string::npos) stem = stem.substr(0, stem.find_last_of('.'));
   std::string dataDir = parts.size() > 1 ? parts[parts.size() - 2] : "data";   // e.g. 23_10_2026
   std::string macro   = __FILE__;                                           // full path of this macro
   std::string base    = outBase[0] ? std::string(outBase) + "/"
                       : macro.substr(0, macro.find_last_of("/\\") + 1);
   std::string outDir  = base + dataDir + "_results/" + stem + "/";         // one folder per data file
   std::string pre = outDir + (PREFIX_OUT ? dataDir + "_" + stem + "_" : "");
   std::string fRoot = pre + "pulse_events.root", fCsv = pre + "pulse_events.csv",
               fWave = pre + "pulse_waveforms", fStat = pre + "pulse_stats",   // plots: extension added on save
               fTime = pre + "pulse_times";

   FILE *f = fopen(dataFile, "rb");
   if (!f) {
      printf("Cannot open %s\n", dataFile);
      if (!gSystem->IsAbsoluteFileName(dataFile))
         printf("  (relative paths start from the folder ROOT was started in: %s)\n"
                "  Give the full path, e.g. root -l 'pmt_events_v2.C(\"D:/data/23_10_2026/wave_0.dat\")'\n",
                gSystem->WorkingDirectory());
      return;
   }
   gSystem->mkdir(outDir.c_str(), true);     // true = also make missing parent folders
   printf("outputs -> %s\n", outDir.c_str());

   // file size via ifstream: ftell is 32-bit on Windows and breaks for files > 2 GB
   long long size = (long long)std::ifstream(dataFile, std::ios::binary | std::ios::ate).tellg();
   Long64_t nEvents = size / (long long)sizeof(Record);
   printf("%s: %lld bytes -> %lld events of %d bins\n", dataFile, size, nEvents, N_SAMP);
   printf("threshold: %g\u03c3  (pulse = more than %g\u03c3 below the baseline, \u03c3 = noise rms of each event)\n",
          THRESH_SIGMA, THRESH_SIGMA);
   if (size % sizeof(Record)) printf("WARNING: %lld leftover bytes (layout may differ).\n", size % (long long)sizeof(Record));
   if (maxEvents > 0 && maxEvents < nEvents) nEvents = maxEvents;
   if (nEvents == 0) { printf("Nothing to read.\n"); fclose(f); return; }

   // ---- output file, tree and csv ----
   TFile *out = new TFile(fRoot.c_str(), "RECREATE");
   TTree *tree = new TTree("pulses", "PMT events containing a pulse");
   Long64_t evt;
   UInt_t   counter, timetag;
   double   bBase, bRms, bHeight, bPeakNs, bStartNs, bEndNs, bWidthNs, bCharge;
   tree->Branch("event",     &evt,     "event/L");
   tree->Branch("counter",   &counter, "counter/i");
   tree->Branch("timetag",   &timetag, "timetag/i");
   tree->Branch("baseline",  &bBase,   "baseline/D");
   tree->Branch("rms",       &bRms,    "rms/D");
   tree->Branch("height",    &bHeight, "height/D");
   tree->Branch("peak_ns",   &bPeakNs, "peak_ns/D");
   tree->Branch("start_ns",  &bStartNs,"start_ns/D");
   tree->Branch("end_ns",    &bEndNs,  "end_ns/D");
   tree->Branch("width_ns",  &bWidthNs,"width_ns/D");
   tree->Branch("charge_pC", &bCharge, "charge_pC/D");

   FILE *csv = WRITE_CSV ? fopen(fCsv.c_str(), "w") : nullptr;
   if (csv) fprintf(csv, "event,counter,timetag,baseline_mV,noise_rms_mV,pulse_height_mV,peak_time_ns,start_ns,end_ns,width_ns,charge_pC\n");

   // ---- scan all events ----
   std::vector<double> allBase, allRms, pHeight, pPeakNs, pCharge;
   std::vector<double> anyPeakNs;            // peak time of every pulse in the whole record
   std::vector<std::vector<double>> shown;   // waveforms of the first N_SHOW pulse events
   std::vector<Long64_t> shownIdx;
   allBase.reserve(nEvents);
   allRms.reserve(nEvents);

   Record rec;
   double v[N_SAMP];
   Long64_t nRead = 0, nPulseAny = 0;
   for (Long64_t i = 0; i < nEvents; i++) {
      if (fread(&rec, sizeof(Record), 1, f) != 1) break;
      nRead++;
      if (i % 10000 == 0) { printf("  scanned %lld / %lld events, %zu with a pulse\r", i, nEvents, pHeight.size()); fflush(stdout); }
      for (int k = 0; k < N_SAMP; k++) v[k] = rec.x[k] * ADC_TO_MV;

      double base, rms, height, charge;
      int peak, iStart, iEnd;
      analyse(v, base, rms, height, peak, iStart, iEnd, charge);
      allBase.push_back(base);                            // noise: all events, no filter
      allRms.push_back(rms);

      double hAny;                                        // whole record: empty events and the all-pulses plot
      int pAny;
      findPeak(v, base, 0, N_SAMP, hAny, pAny);
      if (hAny > THRESH_SIGMA * rms) {
         nPulseAny++;
         anyPeakNs.push_back(pAny * BIN_NS);
      }

      if (height > THRESH_SIGMA * rms) {                  // pulse in the time gate
         evt = i; counter = rec.hdr[4]; timetag = rec.hdr[5];
         bBase = base; bRms = rms; bHeight = height; bPeakNs = peak * BIN_NS;
         bStartNs = iStart * BIN_NS; bEndNs = iEnd * BIN_NS;
         bWidthNs = (iEnd - iStart + 1) * BIN_NS; bCharge = charge;
         tree->Fill();
         if (csv) fprintf(csv, "%lld,%u,%u,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n",
                 evt, counter, timetag, bBase, bRms, bHeight, bPeakNs,
                 bStartNs, bEndNs, bWidthNs, bCharge);
         pHeight.push_back(height);
         pPeakNs.push_back(bPeakNs);
         pCharge.push_back(charge);
         if ((int)shown.size() < N_SHOW) {
            shown.push_back(std::vector<double>(v, v + N_SAMP));
            shownIdx.push_back(i);
         }
      }
   }
   fclose(f);
   if (csv) fclose(csv);
   printf("\n");

   Long64_t nEmpty = nRead - nPulseAny;      // whole record, not limited to the time gate
   printf("processed file      : %s\n", dataFile);
   printf("\n");
   printf("threshold           : %g\u03c3\n", THRESH_SIGMA);
   printf("empty events        : %lld of %lld (%.1f%%)  (no pulse > %g\u03c3 anywhere in 0-%g ns)\n",
          nEmpty, nRead, 100.0 * nEmpty / nRead, THRESH_SIGMA, N_SAMP * BIN_NS);
   printf("pulses in %g-%g ns time gate : %zu  (of %lld events with a pulse anywhere, > %g\u03c3)\n\n",
          GATE_LO * BIN_NS, (GATE_HI - 1) * BIN_NS, pHeight.size(), nPulseAny, THRESH_SIGMA);

   if (pHeight.empty()) {
      printf("No pulses in the time gate (height > %g\u03c3 with the minimum in %g-%g ns).\n",
             THRESH_SIGMA, GATE_LO * BIN_NS, (GATE_HI - 1) * BIN_NS);
      out->Close();
      return;
   }

   // ---- summary table ----
   printf("%-15s%12s%10s%12s%12s%12s   (mV)\n", "", "mean", "std", "median", "min", "max");
   printStats("pulse height", pHeight);
   printStats("charge (pC)", pCharge);
   printStats("noise rms", allRms);
   printStats("noise average", allBase);
   printf("\n");

   gStyle->SetLineScalePS(PDF_LINE_SCALE);   // only affects .pdf output
   gStyle->SetImageScaling(PNG_SCALE);       // only affects .png output

   // ---- waveform plot ----
   int n = shown.size();
   TCanvas *cw = new TCanvas("cw", "Pulse waveforms", 1400, 900);
   cw->Divide(3, (n + 2) / 3);
   std::vector<double> t(N_SAMP);
   for (int k = 0; k < N_SAMP; k++) t[k] = k * BIN_NS;
   for (int k = 0; k < n; k++) {
      cw->cd(k + 1);
      double base, rms, height, charge;
      int peak, iStart, iEnd;
      analyse(shown[k].data(), base, rms, height, peak, iStart, iEnd, charge);
      TGraph *g = new TGraph(N_SAMP, t.data(), shown[k].data());
      g->SetTitle(Form("event %lld: height %.1f mV, rms %.2f mV, Q %.3f pC;time (ns);voltage (mV)",
                       shownIdx[k], height, rms, charge));
      g->SetLineWidth(1);
      g->Draw("AL");
      TLine *bl = new TLine(0, base, N_SAMP * BIN_NS, base);
      bl->SetLineStyle(2);
      bl->Draw();
      TMarker *pk = new TMarker(peak * BIN_NS, shown[k][peak], 23);
      pk->SetMarkerColor(kRed);
      pk->Draw();
      // integration window (green dashed lines)
      gPad->Update();
      double ylo = gPad->GetUymin(), yhi = gPad->GetUymax();
      TLine *l1 = new TLine(iStart * BIN_NS, ylo, iStart * BIN_NS, yhi);
      TLine *l2 = new TLine(iEnd * BIN_NS,   ylo, iEnd * BIN_NS,   yhi);
      for (TLine *l : {l1, l2}) { l->SetLineColor(kGreen + 2); l->SetLineStyle(3); l->Draw(); }
   }
   savePlot(cw, fWave);

   // ---- distributions ----
   // pulse histograms use only the filtered pulses, noise histograms use all events: say so in the titles
   std::string sel = Form("(pulses > %g#sigma, peak %g-%g ns)", THRESH_SIGMA,
                          GATE_CENTER_NS - 0.5 * GATE_WIDTH_NS, GATE_CENTER_NS + 0.5 * GATE_WIDTH_NS);
   std::vector<double> pGain;                // gain = number of electrons = Q / e
   for (double q : pCharge) pGain.push_back(q * 1e-12 / E_CHARGE);

   TH1D *hR = makeHist("hRms",    "Noise RMS (all events, no filter);noise rms (mV);events",     allRms,  100);
   TH1D *hB = makeHist("hBase",   "Baseline (all events, no filter);baseline (mV);events",       allBase, 100);
   TH1D *hH = makeHist("hHeight", ("Pulse height " + sel + ";pulse height (mV);events").c_str(), pHeight, 300);
   TH1D *hQ = makeHist("hCharge", ("Pulse charge " + sel + ";charge (pC);events").c_str(),       pCharge, 300);
   TH1D *hG = makeHist("hGain",   ("Gain = Q/e "   + sel + ";gain (electrons);events").c_str(),  pGain,   300);
   TH1D *hP = makeHist("hPeak",   ("Peak time "    + sel + ";peak time (ns);events").c_str(),    pPeakNs, 100);

   // peak time of all pulses over the whole record: one bin per sample, so no binning spikes
   TH1D *hA = new TH1D("hPeakAll",
                       Form("Peak time, all pulses > %g#sigma (whole record);peak time (ns);events", THRESH_SIGMA),
                       N_SAMP, -0.5 * BIN_NS, (N_SAMP - 0.5) * BIN_NS);
   hA->SetDirectory(nullptr);
   for (double x : anyPeakNs) hA->Fill(x);

   gStyle->SetOptStat("emr");            // entries, mean, rms
   TCanvas *cs = new TCanvas("cs", "Pulse stats", 1500, 800);
   cs->Divide(3, 2);
   cs->cd(1); hR->Draw();
   cs->cd(2); hB->Draw();
   cs->cd(3); hH->Draw();
   cs->cd(4); hQ->Draw();
   cs->cd(5); hG->Draw();
   savePlot(cs, fStat);

   // ---- peak times, in their own window ----
   TCanvas *ct = new TCanvas("ct", "Pulse times", 1400, 550);
   ct->Divide(2, 1);
   ct->cd(1); hA->Draw();
   gPad->Update();                       // time gate (red dashed lines)
   for (double x : {GATE_CENTER_NS - 0.5 * GATE_WIDTH_NS, GATE_CENTER_NS + 0.5 * GATE_WIDTH_NS}) {
      TLine *l = new TLine(x, gPad->GetUymin(), x, gPad->GetUymax());
      l->SetLineColor(kRed); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
   }
   ct->cd(2); hP->Draw();
   savePlot(ct, fTime);

   // ---- save to the .root file ----
   out->cd();
   tree->Write();
   hR->Write(); hB->Write(); hH->Write(); hQ->Write(); hG->Write(); hA->Write(); hP->Write();
   out->Close();
   const char *fmt = SAVE_PNG && SAVE_PDF ? ".png + .pdf" : SAVE_PNG ? ".png" : SAVE_PDF ? ".pdf" : "not saved";
   printf("saved %s%s%s\nplots (%s): %s, %s, %s\n", fRoot.c_str(), WRITE_CSV ? ", " : "", WRITE_CSV ? fCsv.c_str() : "",
          fmt, fWave.c_str(), fStat.c_str(), fTime.c_str());
   if (!SAVE_PNG && !SAVE_PDF) printf("WARNING: SAVE_PNG and SAVE_PDF are both false, no plot files were written.\n");
}