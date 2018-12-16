/**
 *  @file   thesis-plots/chapter_modelling_dQdx/deltaX_plots.C
 *
 *  @brief  Some thesis plots.
 *
 *  $Log: $
 */

#include "TAxis.h"
#include "TLatex.h"
#include "TF2.h"

#include "bethe-faster/BetheFaster.h"
R__LOAD_LIBRARY(libbethe-faster.so)

TCanvas * PlotThreeDDeltaX();
TCanvas * PlotTwoDDeltaX();
TCanvas * PlotRecombinationR();
TCanvas * PlotRecombinationInvR(const double lowX, const double highX, const double lowY, const double highY);
TCanvas * GetNewCanvas(const std::size_t width = 800UL, const std::size_t height = 600UL);

std::size_t g_canvasCounter = 0UL;

//-----------------------------------------------------------------------------------------------------------------------------------------

int deltaX_plots()
{
    bf::PlotHelper::SetGlobalPlotStyle();

    // Set up the propagator
    const auto detector = bf::DetectorHelper::GetMicroBooNEDetector();
    const auto propagator = bf::Propagator{detector};

    PlotTwoDDeltaX()->SaveAs("deltaX_twoDdeltaX.eps");

    PlotRecombinationR()->SaveAs("recomb_R.eps");
    PlotRecombinationInvR(0., 2400., -40., 120.)->SaveAs("recomb_invR.eps");

    gStyle->SetPalette(kBrownCyan);
    PlotThreeDDeltaX()->SaveAs("deltaX_threeDdeltaX.eps");

    return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotThreeDDeltaX()
{
    TCanvas *pCanvas = GetNewCanvas();

    const double pitch = 0.03;
    const double width = 0.06;

    TF2 *pThreeDFunction = new TF2{"threeDFunction", "std::min([0] / std::abs(std::sin(x) * sin(y)), [1] / std::abs(std::cos(x) * sin(y)))", 0., M_PI / 2.,  0.2, M_PI / 2.};
    pThreeDFunction->SetParameter(0, width);
    pThreeDFunction->SetParameter(1, pitch);

    pThreeDFunction->SetLineWidth(1);
    pThreeDFunction->SetLineColor(kBlack);
    pThreeDFunction->SetMaximum(0.3);

    TF1 *pCopy = pThreeDFunction->DrawCopy("surf1");    

    pCopy->GetZaxis()->SetTitle("\\Delta x_i \\text{ (cm)}");
    pCopy->GetYaxis()->SetTitle("\\theta_i \\text{ (radians)}");
    pCopy->GetXaxis()->SetTitle("\\phi_i \\text{ (radians)}");

    pCopy->GetZaxis()->SetTitleOffset(1.5);
    pCopy->GetXaxis()->SetTitleOffset(1.9);
    pCopy->GetYaxis()->SetTitleOffset(1.9);

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotTwoDDeltaX()
{
    TCanvas *pCanvas = GetNewCanvas();

    const double pitch = 0.03;
    const double width = 0.06;

    TF1 cosTermFunction{"cosTermFunction", "[0] / std::abs(std::cos(x))", 0., M_PI / 2.};
    cosTermFunction.SetParameter(0, pitch);
    cosTermFunction.SetLineColor(bf::PlotHelper::GetSchemeColour(5UL));
    cosTermFunction.SetLineStyle(2);
    TF1 *pCopyFunction = cosTermFunction.DrawCopy();

    pCopyFunction->SetMaximum(0.2);
    pCopyFunction->GetYaxis()->SetTitle("L_{xz,i} \\text{ (cm)}");
    pCopyFunction->GetXaxis()->SetTitle("\\phi_i \\text{ (radians)}");

    TF1 sinTermFunction{"sinTermFunction", "[0] / std::abs(std::sin(x))", 0., M_PI / 2.};
    sinTermFunction.SetParameter(0, width);
    sinTermFunction.SetLineColor(bf::PlotHelper::GetSchemeColour(6UL));
    sinTermFunction.SetLineStyle(2);
    sinTermFunction.DrawClone("same");

    TF1 fullFunction{"fullFunction", "std::min([0] / std::abs(std::sin(x)), [1] / std::abs(std::cos(x)))", 0., M_PI / 2.};
    fullFunction.SetParameter(0, width);
    fullFunction.SetParameter(1, pitch);
    fullFunction.SetLineColor(bf::PlotHelper::GetSchemeColour(0UL));
    fullFunction.SetLineStyle(1);
    fullFunction.DrawClone("same");

    TLegend legend{0.66, 0.66, 0.80, 0.88};
    legend.SetBorderSize(1);
    legend.AddEntry(&fullFunction, "L_{xz,i}\\", "l");
    legend.AddEntry(&cosTermFunction, "\\frac{p}{|\\cos \\phi_i|}", "l");
    legend.AddEntry(&sinTermFunction, "\\frac{w_i}{|\\sin \\phi_i|}", "l");
    legend.DrawClone("same");

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotRecombinationR()
{
    TCanvas *pCanvas = GetNewCanvas();

    const double birksA = 0.8;
    const double recombK = 0.0468;
    const double recombRho = 1.4;
    const double recombEpsilon = 0.273;

    TF1 birksFunction{"birksFunction", "[0] / (1. + [1] / ([2] * [3]) * x)", 0., 18.};
    birksFunction.SetParameter(0, birksA);
    birksFunction.SetParameter(1, recombK);
    birksFunction.SetParameter(2, recombRho);
    birksFunction.SetParameter(3, recombEpsilon);
    birksFunction.SetLineColor(bf::PlotHelper::GetSchemeColour(0UL));
    TF1 *pCopyFunc = birksFunction.DrawCopy();

    pCopyFunc->GetYaxis()->SetTitle("R\\");
    pCopyFunc->GetXaxis()->SetTitle("\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}");
    pCopyFunc->SetMinimum(0.);
    pCopyFunc->SetMaximum(1.);
    pCopyFunc->GetXaxis()->SetRangeUser(0., 18.);

    const double modboxA = 0.93;
    const double modboxB = 0.212;

    TF1 modBoxFunction{"modBoxFunction", "std::log([0] + [1] / ([2] * [3]) * x) / ([1] / ([2] * [3]) * x)", 0., 100.};
    modBoxFunction.SetParameter(0, modboxA);
    modBoxFunction.SetParameter(1, modboxB);
    modBoxFunction.SetParameter(2, recombRho);
    modBoxFunction.SetParameter(3, recombEpsilon);
    modBoxFunction.SetLineColor(bf::PlotHelper::GetSchemeColour(1UL));
    modBoxFunction.SetLineStyle(1);
    modBoxFunction.DrawClone("same");

    TLegend legend{0.56, 0.75, 0.88, 0.86};
    legend.SetBorderSize(1);
    legend.AddEntry(&birksFunction, "Birks' model", "l");
    legend.AddEntry(&modBoxFunction, "Modified box model", "l");
    legend.DrawClone("same");

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotRecombinationInvR(const double lowX, const double highX, const double lowY, const double highY)
{
    TCanvas *pCanvas = GetNewCanvas();

    const double birksA = 0.8;
    const double recombK = 0.0468;
    const double recombRho = 1.4;
    const double recombEpsilon = 0.273;
    const double wIon = 23.6e-6;
    const double cParam = 5.076e-3;

    TF1 birksFunction{"birksFunction", "[4] * x / ([0] - [1] / ([2] * [3]) * [4] * x)", lowX, highX};
    birksFunction.SetParameter(0, birksA);
    birksFunction.SetParameter(1, recombK);
    birksFunction.SetParameter(2, recombRho);
    birksFunction.SetParameter(3, recombEpsilon);
    birksFunction.SetParameter(4, wIon / cParam);
    birksFunction.SetLineColor(bf::PlotHelper::GetSchemeColour(0UL));
    TF1 *pCopyFunc = birksFunction.DrawCopy();

    pCopyFunc->GetYaxis()->SetTitle("\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}");
    pCopyFunc->GetXaxis()->SetTitle("\\mathrm{d}Q/\\mathrm{d}x \\text{ (ADC/cm)}");
    pCopyFunc->SetMinimum(lowY);
    pCopyFunc->SetMaximum(highY);

    const double modboxA = 0.93;
    const double modboxB = 0.212;

    TF1 modBoxFunction{"modBoxFunction", "(std::exp([1] / ([2] * [3]) * [4] * x) - [0])/([1] / ([2] * [3]))", lowX, highX};
    modBoxFunction.SetParameter(0, modboxA);
    modBoxFunction.SetParameter(1, modboxB);
    modBoxFunction.SetParameter(2, recombRho);
    modBoxFunction.SetParameter(3, recombEpsilon);
    modBoxFunction.SetParameter(4, wIon / cParam);
    modBoxFunction.SetLineColor(bf::PlotHelper::GetSchemeColour(1UL));
    modBoxFunction.SetLineStyle(1);
    modBoxFunction.DrawClone("same");

    TF1 norecombFunction{"norecombFunction", "[0] * x", lowX, highX};
    norecombFunction.SetParameter(0, wIon / cParam);
    norecombFunction.SetLineColor(bf::PlotHelper::GetSchemeColour(7UL));
    norecombFunction.SetLineStyle(2);
    norecombFunction.DrawClone("same");

    TLegend legend{0.13, 0.72, 0.45, 0.86};
    legend.SetBorderSize(1);
    legend.AddEntry(&birksFunction, "Birks' model", "l");
    legend.AddEntry(&modBoxFunction, "Modified box model", "l");
    legend.AddEntry(&norecombFunction, "No correction", "l");
    legend.DrawClone("same");
    
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * GetNewCanvas(const std::size_t width, const std::size_t height)
{
    const std::string canvasName = "bfCanvas_" + std::to_string(g_canvasCounter++);
    return new TCanvas{canvasName.c_str(), canvasName.c_str(), static_cast<Int_t>(width), static_cast<Int_t>(height)};
}