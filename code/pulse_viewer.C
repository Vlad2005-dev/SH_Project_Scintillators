// pulse_viewer.C
// Step through the pulses one at a time to see what they look like.
// Uses the settings (DATA_FILE, THRESH_SIGMA, time gate, ...) and the analysis of pmt_events.C,
// which must be in the same folder, so it finds exactly the same pulses as that macro.
// aq
// Usage:
//     root -l pulse_viewer.C                                        // DATA_FILE from pmt_events.C
//     root -l 'pulse_viewer.C("D:/data/23_10_2026/wave_0.dat")'     // another file
//     root -l 'pulse_viewer.C(DATA_FILE, 3.0)'                      // only pulses taller than 3 mV
//     root -l 'pulse_viewer.C(DATA_FILE, 0, 488)'                   // only pulses peaking in bin 488
//     root -l 'pulse_viewer.C(DATA_FILE, 0, -1, true)'              // tallest point in the whole record, no time gate
//
// Controls (buttons in the small "Pulses" window, or keys with the mouse over the plot):
//     next: n or right arrow      previous: p or left arrow      +/-100: N / P      save: s
// At the ROOT prompt:
//     showPulse(123)    pulse number 123 of the selection (counting from 0)
//     showEvent(4567)   event number 4567 of the file (the "event" branch of the tree)
//
// Top plot:    whole waveform; time gate yellow (where the pulse is looked for); empty-event count of the
//              whole file (events with no pulse anywhere in 0-204.8 ns, time gate not applied)
// Both plots: fixed voltage range around the baseline if FIXED_Y = true
// Bottom plot: zoom around the peak, every sample marked; baseline blue,
//              THRESH_SIGMA threshold red, integration window green
// Don't .L pmt_events.C in the same ROOT session first (it would be defined twice).

#include "pmt_events.C"
#include <fstream>
#include "TROOT.h"
#include "TBox.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TControlBar.h"
#include "Buttons.h"
#include "KeySymbols.h"

// ---------------- settings you can change ----------------
const bool   VIEW_ALL       = false;  // true: tallest point anywhere in the record (0-204.8 ns); false: pulse in the time gate (as in the tree)
const double MIN_HEIGHT_MV  = 0;      // only show pulses at least this tall (0 = all)
const int    ONLY_PEAK_BIN  = -1;     // only show pulses peaking in this bin, e.g. 488 (-1 = any bin)
const bool   SORT_BY_HEIGHT = true;  // true: tallest pulse first; false: in file order
const double ZOOM_BEFORE_NS = 5.0;    // bottom plot: from this long before the peak ...
const double ZOOM_AFTER_NS  = 15.0;   // ... to this long after it
const bool   FIXED_Y        = false;   // true: same voltage range for every pulse (relative to its baseline); false: fit each waveform
const double Y_BELOW_MV     = 8.0;    // fixed range: from this far below the baseline ...
const double Y_ABOVE_MV     = 4.0;    // ... to this far above it (anything outside is cut off at the plot edge)
// ----------------------------------------------------------

static std::string           gPvFile;     // data file being viewed
static std::vector<Long64_t> gPvEvents;   // event numbers of the selected pulses
static int                   gPvCur = 0;  // pulse currently shown
static std::string           gPvSel;      // selection text for the titles
static std::string           gPvEmpty;    // empty-event count for the plot
static int                   gPvLo = GATE_LO, gPvHi = GATE_HI;   // where the pulse is looked for
static TControlBar          *gPvBar = nullptr;   // the button window

static bool readEvent(Long64_t i, double *v)
{
   std::ifstream in(gPvFile, std::ios::binary);   // streamoff is 64-bit, so files > 2 GB work on Windows too
   in.seekg((std::streamoff)i * sizeof(Record));
   Record rec;
   if (!in.read((char *)&rec, sizeof(Record))) return false;
   for (int k = 0; k < N_SAMP; k++) v[k] = rec.x[k] * ADC_TO_MV;
   return true;
}

void pvKey(Int_t event, Int_t px, Int_t py, TObject *);

static TCanvas *viewerCanvas()
{
   TCanvas *c = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("pv");
   if (c) return c;
   c = new TCanvas("pv", "Pulse viewer", 1200, 850);
   c->Divide(1, 2);
   c->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", nullptr, nullptr, "pvKey(Int_t,Int_t,Int_t,TObject*)");
   return c;
}

// Graph of bins lo..hi-1 of v, deleted automatically when the pad is cleared
static TGraph *waveGraph(const double *v, int lo, int hi, double base)
{
   TGraph *g = new TGraph(hi - lo);
   for (int k = lo; k < hi; k++) g->SetPoint(k - lo, k * BIN_NS, v[k]);
   if (FIXED_Y) {
      g->SetMinimum(base - Y_BELOW_MV);
      g->SetMaximum(base + Y_ABOVE_MV);
   } else {
      double ymin = *std::min_element(v + lo, v + hi), ymax = *std::max_element(v + lo, v + hi);
      double r = std::max(ymax - ymin, 1.0);
      g->SetMinimum(ymin - 0.05 * r);
      g->SetMaximum(ymax + 0.40 * r);    // room for the legend above the trace
   }
   g->SetEditable(kFALSE);               // clicking the plot does not move the points
   g->SetBit(kCanDelete);
   return g;
}

// Box, line or marker drawn on top of the plot: deleted with the pad, and not movable with the mouse
static void decoration(TObject *o)
{
   o->SetBit(kCanDelete);
   o->SetBit(kCannotPick);
   o->Draw();
}

static void legendStyle(TLegend *leg)
{
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.04);
   decoration(leg);
}

void showPulse(int k)
{
   if (gPvEvents.empty()) { printf("No pulses selected.\n"); return; }
   k = std::max(0, std::min(k, (int)gPvEvents.size() - 1));
   gPvCur = k;
   Long64_t ev = gPvEvents[k];
   double v[N_SAMP];
   if (!readEvent(ev, v)) { printf("Cannot read event %lld from %s\n", ev, gPvFile.c_str()); return; }

   double base, rms, height, charge;
   int peak, iStart, iEnd;
   analyse(v, base, rms, height, peak, iStart, iEnd, charge, gPvLo, gPvHi);
   const double tPeak = peak * BIN_NS;
   const double thrV  = base + POLARITY * THRESH_SIGMA * rms;   // threshold level in mV

   TCanvas *c = viewerCanvas();

   // ---- top: whole waveform ----
   c->cd(1);
   gPad->Clear();
   TGraph *g = waveGraph(v, 0, N_SAMP, base);
   g->SetTitle(Form("pulse %d of %zu   (event %lld)   %s;time (ns);voltage (mV)",
                    k, gPvEvents.size(), ev, gPvSel.c_str()));
   g->SetLineColor(kBlack);
   g->Draw("AL");
   gPad->Update();
   double y0 = gPad->GetUymin(), y1 = gPad->GetUymax();
   TBox *bw = new TBox((GATE_LO - 0.5) * BIN_NS, y0, (GATE_HI - 0.5) * BIN_NS, y1);
   bw->SetFillColor(TColor::GetColor("#fff2a8"));   // solid colours: transparent ones upset the line colours in .png
   decoration(bw);
   TLine *lb = new TLine(0, base, N_SAMP * BIN_NS, base);
   lb->SetLineColor(kBlue);
   decoration(lb);
   TLine *lt = new TLine(0, thrV, N_SAMP * BIN_NS, thrV);
   lt->SetLineColor(kRed);
   decoration(lt);
   g->Draw("L");                                       // waveform on top of the boxes
   gPad->RedrawAxis();
   TMarker *m = new TMarker(tPeak, v[peak], 20);
   m->SetMarkerColor(kRed);
   decoration(m);
   TLegend *leg = new TLegend(0.10, 0.75, 0.70, 0.89);
   leg->SetNColumns(3);
   leg->AddEntry(g,  "waveform", "l");
   leg->AddEntry(lb, "baseline (2-42 ns average)", "l");
   leg->AddEntry(lt, Form("threshold (%g#sigma)", THRESH_SIGMA), "l");
   leg->AddEntry(m,  "peak", "p");
   leg->AddEntry(bw, Form("time gate (%g-%g ns)", GATE_LO * BIN_NS, (GATE_HI - 1) * BIN_NS), "f");
   legendStyle(leg);
   TLatex *te = new TLatex(0.72, 0.85, gPvEmpty.c_str());   // whole file, no time gate
   te->SetNDC();
   te->SetTextFont(42);                                   // same font as the legend
   te->SetTextSize(0.04);
   decoration(te);

   // ---- bottom: zoom around the peak ----
   c->cd(2);
   gPad->Clear();
   int lo = std::max(0, peak - (int)std::lround(ZOOM_BEFORE_NS / BIN_NS));
   int hi = std::min(N_SAMP, peak + (int)std::lround(ZOOM_AFTER_NS / BIN_NS) + 1);
   TGraph *z = waveGraph(v, lo, hi, base);
   z->SetTitle(Form("height %.2f mV = %.1f #sigma,   rms %.2f mV,   peak %.1f ns,   Q %.4f pC,   integration %.1f ns"
                    ";time (ns);voltage (mV)",
                    height, height / rms, rms, tPeak, charge, (iEnd - iStart + 1) * BIN_NS));
   z->SetMarkerStyle(20);
   z->SetMarkerSize(0.6);
   z->Draw("APL");
   gPad->Update();
   y0 = gPad->GetUymin(); y1 = gPad->GetUymax();
   TBox *bi = new TBox((iStart - 0.5) * BIN_NS, y0, (iEnd + 0.5) * BIN_NS, y1);
   bi->SetFillColor(TColor::GetColor("#c6efc6"));
   decoration(bi);
   lb = new TLine(lo * BIN_NS, base, (hi - 1) * BIN_NS, base);
   lb->SetLineColor(kBlue);
   decoration(lb);
   lt = new TLine(lo * BIN_NS, thrV, (hi - 1) * BIN_NS, thrV);
   lt->SetLineColor(kRed);
   decoration(lt);
   z->Draw("PL");
   gPad->RedrawAxis();
   leg = new TLegend(0.10, 0.80, 0.80, 0.89);
   leg->SetNColumns(4);
   leg->AddEntry(z,  "samples (every 0.2 ns)", "lp");
   leg->AddEntry(lb, "baseline", "l");
   leg->AddEntry(lt, Form("threshold (%g#sigma)", THRESH_SIGMA), "l");
   leg->AddEntry(bi, "integration window (charge)", "f");
   legendStyle(leg);

   c->GetPad(1)->Modified();                 // make sure both plots repaint when called from a button or key
   c->GetPad(2)->Modified();
   c->Modified();
   c->Update();
   gSystem->ProcessEvents();
   printf("pulse %d/%zu  event %lld  height %.2f mV (%.1f\u03c3)  peak %.1f ns (bin %d)  Q %.4f pC\n",
          k, gPvEvents.size(), ev, height, height / rms, tPeak, peak, charge);
}

void nextPulse()         { showPulse(gPvCur + 1); }
void prevPulse()         { showPulse(gPvCur - 1); }
void jumpPulses(int d)   { showPulse(gPvCur + d); }

void showEvent(Long64_t ev)
{
   auto it = std::find(gPvEvents.begin(), gPvEvents.end(), ev);
   if (it == gPvEvents.end()) { printf("Event %lld is not in the selection.\n", ev); return; }
   showPulse((int)(it - gPvEvents.begin()));
}

// Save the plot shown as <macro folder>/<data folder>_results/<data file>/pulses/event_<n>.png
void savePulse()
{
   if (gPvEvents.empty()) return;
   std::vector<std::string> parts = pathParts(gPvFile.c_str());
   std::string stem = parts.back();
   if (stem.find_last_of('.') != std::string::npos) stem = stem.substr(0, stem.find_last_of('.'));
   std::string dataDir = parts.size() > 1 ? parts[parts.size() - 2] : "data";
   std::string macro = __FILE__;
   std::string dir = macro.substr(0, macro.find_last_of("/\\") + 1) + dataDir + "_results/" + stem + "/pulses/";
   gSystem->mkdir(dir.c_str(), true);
   viewerCanvas()->SaveAs(Form("%sevent_%lld.png", dir.c_str(), gPvEvents[gPvCur]));
}

// Keys pressed with the mouse over the plot
void pvKey(Int_t event, Int_t px, Int_t py, TObject *)
{
   if (event != kKeyPress) return;
   if (px == 'n' || py == kKey_Right)     nextPulse();
   else if (px == 'p' || py == kKey_Left) prevPulse();
   else if (px == 'N')                    jumpPulses(100);
   else if (px == 'P')                    jumpPulses(-100);
   else if (px == 's')                    savePulse();
}

void pulse_viewer(const char *dataFile = DATA_FILE, double minHeight = MIN_HEIGHT_MV,
                  int peakBin = ONLY_PEAK_BIN, bool viewAll = VIEW_ALL)
{
   gPvFile = dataFile;
   FILE *f = fopen(dataFile, "rb");
   if (!f) { printf("Cannot open %s\n", dataFile); return; }

   // ---- find the pulses (same cuts as pmt_events.C, plus the extra ones above) ----
   gPvLo = viewAll ? 0 : GATE_LO;            // where the pulse is looked for
   gPvHi = viewAll ? N_SAMP : GATE_HI;
   gPvEvents.clear();
   std::vector<double> heights;
   Record rec;
   double v[N_SAMP];
   Long64_t i = 0, nPulseAny = 0;
   while (fread(&rec, sizeof(Record), 1, f) == 1) {
      if (i % 50000 == 0) { printf("  scanned %lld events, %zu pulses selected\r", i, gPvEvents.size()); fflush(stdout); }
      for (int k = 0; k < N_SAMP; k++) v[k] = rec.x[k] * ADC_TO_MV;
      double base, rms, height, charge;
      int peak, iStart, iEnd;
      analyse(v, base, rms, height, peak, iStart, iEnd, charge, gPvLo, gPvHi);
      double hAny;                          // for the empty events: whole record, no time gate or viewer cuts
      int pAny;
      findPeak(v, base, 0, N_SAMP, hAny, pAny);
      if (hAny > THRESH_SIGMA * rms) nPulseAny++;
      bool keep = height > THRESH_SIGMA * rms
               && height >= minHeight
               && (peakBin < 0 || peak == peakBin);
      if (keep) { gPvEvents.push_back(i); heights.push_back(height); }
      i++;
   }
   fclose(f);

   if (SORT_BY_HEIGHT) {
      std::vector<size_t> idx(gPvEvents.size());
      for (size_t j = 0; j < idx.size(); j++) idx[j] = j;
      std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return heights[a] > heights[b]; });
      std::vector<Long64_t> sorted;
      for (size_t j : idx) sorted.push_back(gPvEvents[j]);
      gPvEvents = sorted;
   }

   gPvSel = Form("> %g#sigma", THRESH_SIGMA);
   gPvSel += viewAll ? ", whole record" : Form(", peak %g-%g ns", GATE_LO * BIN_NS, (GATE_HI - 1) * BIN_NS);
   if (minHeight > 0) gPvSel += Form(", > %g mV", minHeight);
   if (peakBin >= 0)  gPvSel += Form(", bin %d only", peakBin);
   if (SORT_BY_HEIGHT) gPvSel += ", tallest first";
   TString selText = gPvSel;                // terminal version: ROOT's #sigma -> real sigma sign
   selText.ReplaceAll("#sigma", "\u03c3");
   printf("\nprocessed file: %s\n", dataFile);
   printf("threshold: %g\u03c3  (pulse = more than %g\u03c3 below the baseline, \u03c3 = noise rms of each event)\n",
          THRESH_SIGMA, THRESH_SIGMA);
   printf("%lld events scanned, %zu pulses selected (%s)\n", i, gPvEvents.size(), selText.Data());
   Long64_t nEmpty = i - nPulseAny;
   printf("empty events: %lld of %lld (%.1f%%)  (no pulse > %g\u03c3 anywhere in 0-%g ns; time gate not applied)\n",
          nEmpty, i, i ? 100.0 * nEmpty / i : 0.0, THRESH_SIGMA, N_SAMP * BIN_NS);
   gPvEmpty = Form("empty events: %lld of %lld (%.1f%%)", nEmpty, i, i ? 100.0 * nEmpty / i : 0.0);
   if (gPvEvents.empty()) return;

   if (!gROOT->IsBatch()) {
      delete gPvBar;                       // old button window from an earlier run
      gPvBar = new TControlBar("vertical", "Pulses", 20, 20);
      TControlBar *bar = gPvBar;
      bar->AddButton("next  >",       "nextPulse()",     "next pulse (key n or right arrow)");
      bar->AddButton("<  previous",   "prevPulse()",     "previous pulse (key p or left arrow)");
      bar->AddButton("+100  >>",      "jumpPulses(100)", "100 pulses forward (key N)");
      bar->AddButton("<<  -100",      "jumpPulses(-100)","100 pulses back (key P)");
      bar->AddButton("save png",      "savePulse()",     "save this plot (key s)");
      bar->Show();
      printf("keys over the plot: n/p or arrows = next/previous, N/P = +/-100, s = save.\n"
             "prompt: showPulse(k), showEvent(event number)\n");
   }
   showPulse(0);
}
