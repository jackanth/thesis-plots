/**
 *  @file   thesis-plots/chapter_modelling_dQdx/bethe_plots.C
 *
 *  @brief  Some thesis plots.
 *
 *  $Log: $
 */

#include "TAxis.h"
#include "TLatex.h"

#include "bethe-faster/BetheFaster.h"
R__LOAD_LIBRARY(libbethe-faster.so)

TCanvas * PlotdEdxVersusXMean(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle, const unsigned int colour);
TCanvas * PlotEnergyMean(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle, const unsigned int colour);
TCanvas * PlotEnergies(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode);
TCanvas * PlotdEdxVersusX(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode);
TCanvas * PlotdEdxVersusT(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode);
TCanvas * PlotOverlayGraph(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle, const unsigned int colour, const std::string &darkColour, const std::string &label);
TCanvas * GetNewCanvas(const std::size_t width = 800UL, const std::size_t height = 600UL);

std::size_t g_canvasCounter = 0UL;

//-----------------------------------------------------------------------------------------------------------------------------------------

int bethe_plots()
{
    bf::PlotHelper::SetGlobalPlotStyle();

    // Set up the propagator
    const auto detector = bf::DetectorHelper::GetMicroBooNEDetector();
    const auto propagator = bf::Propagator{detector};

    // Bethe equation discussion
    PlotdEdxVersusXMean(propagator, bf::ParticleHelper::GetMuon(), 0UL)->SaveAs("bethe_dEdxVersusX_mean_muon.eps");
    PlotEnergyMean(propagator, bf::ParticleHelper::GetMuon(), 0UL)->SaveAs("bethe_energy_mean_muon.eps");

    // Modal energy loss discussion
    PlotEnergies(propagator, bf::Propagator::PROPAGATION_MODE::MEAN)->SaveAs("modal_energies_mean.eps");
    PlotEnergies(propagator, bf::Propagator::PROPAGATION_MODE::MODAL)->SaveAs("modal_energies_mode.eps");
    PlotdEdxVersusX(propagator, bf::Propagator::PROPAGATION_MODE::MEAN)->SaveAs("modal_dEdxVersusX_mean.eps");
    PlotdEdxVersusX(propagator, bf::Propagator::PROPAGATION_MODE::MODAL)->SaveAs("modal_dEdxVersusX_mode.eps");
    PlotdEdxVersusT(propagator, bf::Propagator::PROPAGATION_MODE::MODAL)->SaveAs("modal_dEdxVersusT_mode.eps");

    TStyle *pStyle = gROOT->GetStyle("BetheFasterStyle");

    pStyle->SetLabelSize(0.065f, "xyz");
    pStyle->SetTitleSize(0.075f, "xyz");

    pStyle->SetPadBottomMargin(0.2f);
    pStyle->SetPadTopMargin(0.1f);
    pStyle->SetPadRightMargin(0.08f);

    PlotOverlayGraph(propagator, bf::ParticleHelper::GetMuon(), 0UL, "#178262", "\\mu")->SaveAs("modal_overlay_muon.eps");
    PlotOverlayGraph(propagator, bf::ParticleHelper::GetChargedPion(), 1UL, "#B24E02", "\\pi^\\pm")->SaveAs("modal_overlay_charged_pion.eps");
    PlotOverlayGraph(propagator, bf::ParticleHelper::GetChargedKaon(), 2UL, "#605C93", "K^\\pm")->SaveAs("modal_overlay_charged_kaon.eps");
    PlotOverlayGraph(propagator, bf::ParticleHelper::GetProton(), 3UL, "#BE2271", "p\\")->SaveAs("modal_overlay_proton.eps");
    
    return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotEnergies(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode)
{
    // Get the particles
    const auto &spMuon = bf::ParticleHelper::GetMuon();
    const auto &spChargedPion = bf::ParticleHelper::GetChargedPion();
    const auto &spChargedKaon = bf::ParticleHelper::GetChargedKaon();
    const auto &spProton = bf::ParticleHelper::GetProton();

    // Propagate the particles
    while ((spMuon->KineticEnergy() < 1000.) && !spMuon->HasFailed())
        propagator.PropagateBackwards(spMuon, 0.03, mode);

    while ((spChargedPion->KineticEnergy() < 1000.) && !spChargedPion->HasFailed())
        propagator.PropagateBackwards(spChargedPion, 0.03, mode);

    while ((spChargedKaon->KineticEnergy() < 1000.) && !spChargedKaon->HasFailed())
        propagator.PropagateBackwards(spChargedKaon, 0.03, mode);

    while ((spProton->KineticEnergy() < 1000.) && !spProton->HasFailed())
        propagator.PropagateBackwards(spProton, 0.03, mode);

    // Get the graphs
    auto muonGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleEnergyGraph(spMuon, false), "\\mu", 0UL, true, 2};
    auto chargedPionGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleEnergyGraph(spChargedPion, false), "\\pi^\\pm", 1UL, true, 2};
    auto chargedKaonGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleEnergyGraph(spChargedKaon, false), "K^\\pm", 2UL, true, 2};
    auto protonGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleEnergyGraph(spProton, false), "p\\", 3UL, true, 2};

    // Create the graphs
    bf::PlotOptions options;
    options.m_xAxisTitle = "x \\text{ (cm)}";
    options.m_yAxisTitle = "T \\text{ (MeV/cm)}";

    TCanvas *pCanvas = GetNewCanvas();

    TMultiGraph *pMultiGraph = nullptr;
    TLegend *pLegend = nullptr;

    auto graphs = std::vector<std::reference_wrapper<bf::MultiGraphEntry>>{muonGraph, chargedPionGraph, chargedKaonGraph, protonGraph};
    bf::PlotHelper::GetMultiGraph(graphs, options, pMultiGraph, pLegend);

    pMultiGraph->GetXaxis()->SetLimits(0., 700.);
    pMultiGraph->SetMinimum(0.);
    pMultiGraph->SetMaximum(1000.);

    pMultiGraph->DrawClone("A");
    pLegend->DrawClone();

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotdEdxVersusX(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode)
{
    // Get the particles
    const auto &spMuon = bf::ParticleHelper::GetMuon();
    const auto &spChargedPion = bf::ParticleHelper::GetChargedPion();
    const auto &spChargedKaon = bf::ParticleHelper::GetChargedKaon();
    const auto &spProton = bf::ParticleHelper::GetProton();

    // Propagate the particles
    while ((spMuon->KineticEnergy() < 10000.) && !spMuon->HasFailed())
        propagator.PropagateBackwards(spMuon, 0.03, mode);

    while ((spChargedPion->KineticEnergy() < 10000.) && !spChargedPion->HasFailed())
        propagator.PropagateBackwards(spChargedPion, 0.03, mode);

    while ((spChargedKaon->KineticEnergy() < 10000.) && !spChargedKaon->HasFailed())
        propagator.PropagateBackwards(spChargedKaon, 0.03, mode);

    while ((spProton->KineticEnergy() < 10000.) && !spProton->HasFailed())
        propagator.PropagateBackwards(spProton, 0.03, mode);

    // Get the graphs
    auto muonGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spMuon, true), "\\mu", 0UL, true, 2};
    auto chargedPionGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spChargedPion, true), "\\pi^\\pm", 1UL, true, 2};
    auto chargedKaonGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spChargedKaon, true), "K^\\pm", 2UL, true, 2};
    auto protonGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spProton, true), "p\\", 3UL, true, 2};

    // Create the graphs
    bf::PlotOptions options;
    options.m_xAxisTitle = "\\text{Residual range (cm)}";
    options.m_yAxisTitle = "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}";

    TCanvas *pCanvas = GetNewCanvas();

    TMultiGraph *pMultiGraph = nullptr;
    TLegend *pLegend = nullptr;

    auto graphs = std::vector<std::reference_wrapper<bf::MultiGraphEntry>>{muonGraph, chargedPionGraph, chargedKaonGraph, protonGraph};
    bf::PlotHelper::GetMultiGraph(graphs, options, pMultiGraph, pLegend);

    pMultiGraph->GetXaxis()->SetLimits(0., 1000.);
    pMultiGraph->SetMinimum(1.);
    pMultiGraph->SetMaximum(5.);

    pMultiGraph->DrawClone("A");
    pLegend->DrawClone();

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotdEdxVersusT(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode)
{
    // Get the particles
    const auto &spMuon = bf::ParticleHelper::GetMuon();
    const auto &spChargedPion = bf::ParticleHelper::GetChargedPion();
    const auto &spChargedKaon = bf::ParticleHelper::GetChargedKaon();
    const auto &spProton = bf::ParticleHelper::GetProton();

    // Propagate the particles
    while ((spMuon->KineticEnergy() < 10000.) && !spMuon->HasFailed())
        propagator.PropagateBackwards(spMuon, 0.03, mode);

    while ((spChargedPion->KineticEnergy() < 10000.) && !spChargedPion->HasFailed())
        propagator.PropagateBackwards(spChargedPion, 0.03, mode);

    while ((spChargedKaon->KineticEnergy() < 10000.) && !spChargedKaon->HasFailed())
        propagator.PropagateBackwards(spChargedKaon, 0.03, mode);

    while ((spProton->KineticEnergy() < 10000.) && !spProton->HasFailed())
        propagator.PropagateBackwards(spProton, 0.03, mode);

    // Get the graphs
    auto muonGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusTGraph(spMuon), "\\mu", 0UL, true, 2};
    auto chargedPionGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusTGraph(spChargedPion), "\\pi^\\pm", 1UL, true, 2};
    auto chargedKaonGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusTGraph(spChargedKaon), "K^\\pm", 2UL, true, 2};
    auto protonGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusTGraph(spProton), "p\\", 3UL, true, 2};

    // Create the graphs
    bf::PlotOptions options;
    options.m_xAxisTitle = "T \\text{ (MeV)}";
    options.m_yAxisTitle = "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}";

    TCanvas *pCanvas = GetNewCanvas();

    TMultiGraph *pMultiGraph = nullptr;
    TLegend *pLegend = nullptr;

    auto graphs = std::vector<std::reference_wrapper<bf::MultiGraphEntry>>{muonGraph, chargedPionGraph, chargedKaonGraph, protonGraph};
    bf::PlotHelper::GetMultiGraph(graphs, options, pMultiGraph, pLegend);

    pMultiGraph->GetXaxis()->SetLimits(0., 2000.);
    pMultiGraph->SetMinimum(1.);
    pMultiGraph->SetMaximum(5.);

    pMultiGraph->DrawClone("A");
    pLegend->DrawClone();

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotdEdxVersusXMean(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle, const unsigned int colour)
{
    // Propagate the particle
    while ((spParticle->KineticEnergy() < 1000.) && !spParticle->HasFailed())
        propagator.PropagateBackwards(spParticle, 0.03, bf::Propagator::PROPAGATION_MODE::MEAN);

    TCanvas *pCanvas = GetNewCanvas();

    TGraph graph = bf::PlotHelper::GetParticledEdxVersusXGraph(spParticle, false);
    graph.GetYaxis()->SetTitle("-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}");
    graph.GetXaxis()->SetTitle("x\\text{ (cm)}");
    graph.SetLineColor(bf::PlotHelper::GetSchemeColour(colour));
    graph.SetLineWidth(2);

    graph.GetXaxis()->SetLimits(0., 500.);
    graph.SetMinimum(2.);
    graph.SetMaximum(5.);

    graph.DrawClone("AL");

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotEnergyMean(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle, const unsigned int colour)
{
    // Propagate the particle
    while ((spParticle->KineticEnergy() < 1000.) && !spParticle->HasFailed())
        propagator.PropagateBackwards(spParticle, 0.03, bf::Propagator::PROPAGATION_MODE::MEAN);

    TCanvas *pCanvas = GetNewCanvas();

    TGraph graph = bf::PlotHelper::GetParticleEnergyGraph(spParticle, false);
    graph.GetYaxis()->SetTitle("T \\text{ (MeV)}");
    graph.GetXaxis()->SetTitle("x\\text{ (cm)}");
    graph.SetLineColor(bf::PlotHelper::GetSchemeColour(colour));
    graph.SetLineWidth(2);

    graph.GetXaxis()->SetLimits(0., 500.);
    graph.SetMinimum(0.);
    graph.SetMaximum(1000.);

    graph.DrawClone("AL");

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotOverlayGraph(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle, const unsigned int colour, const std::string &darkColour, const std::string &label)
{
    // Get the stochastic dE/dx graph
    while ((spParticle->KineticEnergy() < 1000.) && !spParticle->HasFailed())
        propagator.PropagateBackwards(spParticle, 0.03);

    auto stochasticdEdxGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spParticle, true), "", colour, false, 1};

    // Get the modal dE/dx graph
    spParticle->Reset();

    while ((spParticle->KineticEnergy() < 1000.) && !spParticle->HasFailed())
        propagator.PropagateBackwards(spParticle, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    auto modaldEdxGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spParticle, true), "", colour, true, 2};

    // Get mean mean dE/dx graph
    spParticle->Reset();

    while ((spParticle->KineticEnergy() < 1000.) && !spParticle->HasFailed())
        propagator.PropagateBackwards(spParticle, 0.03, bf::Propagator::PROPAGATION_MODE::MEAN);

    auto meandEdxGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spParticle, true), "", colour, true, 2};

    // Create the graphs
    bf::PlotOptions options;
    options.m_xAxisTitle = "\\text{Residual range (cm)}";
    options.m_yAxisTitle = "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}";

    TCanvas *pCanvas = GetNewCanvas(900UL, 300UL);

    TMultiGraph *pMultiGraph = nullptr;
    TLegend *pLegend = nullptr;

    auto graphs = std::vector<std::reference_wrapper<bf::MultiGraphEntry>>{stochasticdEdxGraph, modaldEdxGraph, meandEdxGraph};
    bf::PlotHelper::GetMultiGraph(graphs, options, pMultiGraph, pLegend);

    // Extra formatting
    const auto darkColor = TColor::GetColor(darkColour.c_str());
    meandEdxGraph.Graph().SetLineStyle(2);
    modaldEdxGraph.Graph().SetLineColor(darkColor);
    meandEdxGraph.Graph().SetLineColor(darkColor);

    pMultiGraph->GetXaxis()->SetLimits(0., 700.);
    pMultiGraph->SetMaximum(10.);
    pMultiGraph->GetYaxis()->SetTitleOffset(0.4);

    pMultiGraph->DrawClone("A");
   
    // Draw the label
    TLatex latex;
    latex.SetTextSize(0.1);
    latex.DrawLatex(635., 8.5, label.c_str());

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * GetNewCanvas(const std::size_t width, const std::size_t height)
{
    const std::string canvasName = "bfCanvas_" + std::to_string(g_canvasCounter++);
    return new TCanvas{canvasName.c_str(), canvasName.c_str(), static_cast<Int_t>(width), static_cast<Int_t>(height)};
}