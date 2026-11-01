// dark_count.C
// Dark count vs discriminator threshold.
//
//     root -l dark_count.C        (run from the code/ folder)
//
// Error bars are Poisson: sqrt(N).
// Set COUNT_TIME_S to the counting time per point to plot a rate in Hz instead of raw counts.
//
// Output: ../results/dark_count/dark_count_rate.png

#include <cmath>
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TSystem.h"

// ---------------- settings you can change ----------------
const double COUNT_TIME_S = 60;      // counting time per point in s; 0 = plot raw counts
const char  *OUT_DIR      = "../results/dark_count";   // relative to the folder ROOT is started in
// ----------------------------------------------------------

void dark_count()
{
   const int n = 10;
   double thr[n]    = { -2, -4, -6, -8, -10, -12, -14, -16, -18, -20 };   // mV
   double counts[n] = { 73, 26, 14, 18,  11,  10,   5,   7,  12,   6 };   // -4 mV: second reading was 37

   double y[n], ey[n], ex[n];
   const double scale = COUNT_TIME_S > 0 ? 1.0 / COUNT_TIME_S : 1.0;
   for (int i = 0; i < n; i++) {
      y[i]  = counts[i] * scale;
      ey[i] = std::sqrt(counts[i]) * scale;
      ex[i] = 0;
   }

   TCanvas *c = new TCanvas("c", "Dark count rate", 900, 650);
   c->SetGrid();
   TGraphErrors *g = new TGraphErrors(n, thr, y, ex, ey);
   g->SetTitle(Form("Dark counts vs threshold;threshold (mV);%s",
                    COUNT_TIME_S > 0 ? "dark count rate (Hz)" : "counts"));
   g->SetMarkerStyle(20);
   g->SetMarkerSize(1.2);
   g->SetLineWidth(2);
   g->Draw("APL");
   g->GetYaxis()->SetRangeUser(0, (73 + 3 * std::sqrt(73.0)) * scale);
   gSystem->mkdir(OUT_DIR, true);
   c->SaveAs(Form("%s/dark_count_rate.png", OUT_DIR));
}
