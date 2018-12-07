/**
 *  @file   thesis-plots/chapter_modelling_dQdx/bethe_plots.C
 *
 *  @brief  Some thesis plots.
 *
 *  $Log: $
 */

#include "TAxis.h"
#include "TLatex.h"
#include "TGaxis.h"

#include "bethe-faster/BetheFaster.h"
R__LOAD_LIBRARY(libbethe-faster.so)

TCanvas * PlotdEdxVersusXMean(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle, const unsigned int colour);
TCanvas * PlotEnergyMean(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle, const unsigned int colour);
TCanvas * PlotEnergies(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode);
TCanvas * PlotdEdxVersusX(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode);
TCanvas * PlotdEdxVersusT(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode);
TGraph GetNoDensityEffectErrorGraph(const bf::Detector &detector, const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle);
TCanvas * PlotNoDensityEffectError(const bf::Detector &detector, const bf::Propagator &propagator);
TGraph GetFermiPlateauErrorGraph(const bf::Detector &detector, const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle);
TCanvas * PlotFermiPlateauError(const bf::Detector &detector, const bf::Propagator &propagator);
TGraph GetResidualRangeVersusScaledTGraph(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle);
TCanvas * PlotResidualRangeVersusScaledT(const bf::Propagator &propagator);
TGraph * GetLowEnergyApproxFractionalErrorGraph(const std::shared_ptr<bf::Particle> &spParticle, const std::function<double(double)> &minusdEdxGetter,const double range);
TCanvas * PlotLowEnergyApproxError(const bf::Detector &detector, const bf::Propagator &propagator, const bf::QuickPidAlgorithm &quickPidAlg, const std::shared_ptr<bf::Particle> &spParticle, const unsigned int colour, const std::string &label, const double deltaX, const double range);
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

    // // Modal energy loss discussion
    PlotEnergies(propagator, bf::Propagator::PROPAGATION_MODE::MEAN)->SaveAs("modal_energies_mean.eps");
    PlotEnergies(propagator, bf::Propagator::PROPAGATION_MODE::MODAL)->SaveAs("modal_energies_mode.eps");
    PlotdEdxVersusX(propagator, bf::Propagator::PROPAGATION_MODE::MEAN)->SaveAs("modal_dEdxVersusX_mean.eps");
    PlotdEdxVersusX(propagator, bf::Propagator::PROPAGATION_MODE::MODAL)->SaveAs("modal_dEdxVersusX_mode.eps");
    PlotdEdxVersusT(propagator, bf::Propagator::PROPAGATION_MODE::MODAL)->SaveAs("modal_dEdxVersusT_mode.eps");

    // // Density effect discussion
    PlotNoDensityEffectError(detector, propagator)->SaveAs("densityeffect_error_nodensityeffect.eps");
    PlotFermiPlateauError(detector, propagator)->SaveAs("densityeffect_error_fermiplateau.eps");

    // // Low-energy epproximation discussion
    const auto quickPidAlg = bf::QuickPidAlgorithm{detector};
    PlotResidualRangeVersusScaledT(propagator)->SaveAs("lowenergyapprox_RVersusScaledT.eps");
    PlotLowEnergyApproxError(detector, propagator, quickPidAlg, bf::ParticleHelper::GetMuon(), 0UL, "\\mu", 0.03, 11.)->SaveAs("lowenergyapprox_lowEnergyApproxError_0.03cm_muon.eps");
    PlotLowEnergyApproxError(detector, propagator, quickPidAlg, bf::ParticleHelper::GetChargedPion(), 1UL, "\\pi^\\pm", 0.03, 11.)->SaveAs("lowenergyapprox_lowEnergyApproxError_0.03cm_charged_pion.eps");
    PlotLowEnergyApproxError(detector, propagator, quickPidAlg, bf::ParticleHelper::GetChargedKaon(), 2UL, "K^\\pm", 0.03, 11.)->SaveAs("lowenergyapprox_lowEnergyApproxError_0.03cm_charged_kaon.eps");
    PlotLowEnergyApproxError(detector, propagator, quickPidAlg, bf::ParticleHelper::GetProton(), 3UL, "p\\", 0.03, 11.)->SaveAs("lowenergyapprox_lowEnergyApproxError_0.03cm_proton.eps");
    PlotLowEnergyApproxError(detector, propagator, quickPidAlg, bf::ParticleHelper::GetMuon(), 0UL, "\\mu", 0.3, 11.)->SaveAs("lowenergyapprox_lowEnergyApproxError_0.3cm_muon.eps");
    PlotLowEnergyApproxError(detector, propagator, quickPidAlg, bf::ParticleHelper::GetChargedPion(), 1UL, "\\pi^\\pm", 0.3, 11.)->SaveAs("lowenergyapprox_lowEnergyApproxError_0.3cm_charged_pion.eps");
    PlotLowEnergyApproxError(detector, propagator, quickPidAlg, bf::ParticleHelper::GetChargedKaon(), 2UL, "K^\\pm", 0.3, 11.)->SaveAs("lowenergyapprox_lowEnergyApproxError_0.3cm_charged_kaon.eps");
    PlotLowEnergyApproxError(detector, propagator, quickPidAlg, bf::ParticleHelper::GetProton(), 3UL, "p\\", 0.3, 11.)->SaveAs("lowenergyapprox_lowEnergyApproxError_0.3cm_proton.eps");

    // Change formatting for the overlay graphs
    TStyle *pStyle = gROOT->GetStyle("BetheFasterStyle");

    pStyle->SetLabelSize(0.065f, "xyz");
    pStyle->SetTitleSize(0.075f, "xyz");

    pStyle->SetPadBottomMargin(0.2f);
    pStyle->SetPadTopMargin(0.1f);
    pStyle->SetPadRightMargin(0.08f);

    // Overlay graphs for modal energy loss discussion
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

TGraph GetNoDensityEffectErrorGraph(const bf::Detector &detector, const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    std::vector<double> energies, noDensityEffectErrors;
    const auto &        history = spParticle->GetHistory();

    if (history.empty())
        throw std::runtime_error{"Cannot plot empty particle history"};

    for (const auto &spState : history)
    {
        const double densityCorrection = propagator.DensityCorrection(spParticle->Mass(), spState->GetKineticEnergy());
        const double beta = bf::ParticleHelper::GetParticleBeta(spParticle->Mass(), spState->GetKineticEnergy());
        const double gamma = bf::ParticleHelper::GetParticleGamma(spParticle->Mass(), spState->GetKineticEnergy());
        const double xi = propagator.Xi(spParticle->Mass(), spState->GetKineticEnergy(), spState->Getdx());
        const double denominator = std::log(2. * bf::PhysicalConstants::m_electronMass * beta * beta * gamma * gamma * xi / (detector.m_avgIonizationEnergy * detector.m_avgIonizationEnergy * 1.e-12)) + 0.2 - beta * beta - densityCorrection;

        const double noDensityEffectError = std::abs(densityCorrection / denominator);

        energies.push_back(spState->GetKineticEnergy());
        noDensityEffectErrors.push_back(noDensityEffectError);
    }

    return TGraph{static_cast<Int_t>(history.size()), energies.data(), noDensityEffectErrors.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotNoDensityEffectError(const bf::Detector &detector, const bf::Propagator &propagator)
{
    // Get the particles
    const auto &spMuon = bf::ParticleHelper::GetMuon();
    const auto &spChargedPion = bf::ParticleHelper::GetChargedPion();
    const auto &spChargedKaon = bf::ParticleHelper::GetChargedKaon();
    const auto &spProton = bf::ParticleHelper::GetProton();

    // Propagate the particles
    while ((spMuon->KineticEnergy() < 10000.) && !spMuon->HasFailed())
        propagator.PropagateBackwards(spMuon, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spChargedPion->KineticEnergy() < 10000.) && !spChargedPion->HasFailed())
        propagator.PropagateBackwards(spChargedPion, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spChargedKaon->KineticEnergy() < 10000.) && !spChargedKaon->HasFailed())
        propagator.PropagateBackwards(spChargedKaon, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spProton->KineticEnergy() < 10000.) && !spProton->HasFailed())
        propagator.PropagateBackwards(spProton, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    // Get the graphs
    auto muonGraph = bf::MultiGraphEntry{GetNoDensityEffectErrorGraph(detector, propagator, spMuon), "\\mu", 0UL, true, 2};
    auto chargedPionGraph = bf::MultiGraphEntry{GetNoDensityEffectErrorGraph(detector, propagator, spChargedPion), "\\pi^\\pm", 1UL, true, 2};
    auto chargedKaonGraph = bf::MultiGraphEntry{GetNoDensityEffectErrorGraph(detector, propagator, spChargedKaon), "K^\\pm", 2UL, true, 2};
    auto protonGraph = bf::MultiGraphEntry{GetNoDensityEffectErrorGraph(detector, propagator, spProton), "p\\", 3UL, true, 2};

    // Create the graphs
    bf::PlotOptions options;
    options.m_xAxisTitle = "T \\text{ (MeV)}";
    options.m_yAxisTitle = "R_0(T)\\";
    options.m_legendX1 = 0.14;
    options.m_legendX2 = 0.22;

    TCanvas *pCanvas = GetNewCanvas();

    TMultiGraph *pMultiGraph = nullptr;
    TLegend *pLegend = nullptr;

    auto graphs = std::vector<std::reference_wrapper<bf::MultiGraphEntry>>{muonGraph, chargedPionGraph, chargedKaonGraph, protonGraph};
    bf::PlotHelper::GetMultiGraph(graphs, options, pMultiGraph, pLegend);

    pMultiGraph->GetXaxis()->SetLimits(0., 10000.);

    pMultiGraph->DrawClone("A");
    pLegend->DrawClone();

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TGraph GetFermiPlateauErrorGraph(const bf::Detector &detector, const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    std::vector<double> energies, fermiPlateauErrors;
    const auto &        history = spParticle->GetHistory();

    if (history.empty())
        throw std::runtime_error{"Cannot plot empty particle history"};

    for (const auto &spState : history)
    {
        const double densityCorrection = propagator.DensityCorrection(spParticle->Mass(), spState->GetKineticEnergy());
        const double beta = bf::ParticleHelper::GetParticleBeta(spParticle->Mass(), spState->GetKineticEnergy());
        const double gamma = bf::ParticleHelper::GetParticleGamma(spParticle->Mass(), spState->GetKineticEnergy());
        const double xi = propagator.Xi(spParticle->Mass(), spState->GetKineticEnergy(), spState->Getdx());

        const double numerator = densityCorrection - 2. * std::log(beta * gamma * detector.m_plasmaEnergy / detector.m_avgIonizationEnergy) + beta * beta;
        const double denominator = std::log(2. * bf::PhysicalConstants::m_electronMass * beta * beta * gamma * gamma * xi / (detector.m_avgIonizationEnergy * detector.m_avgIonizationEnergy * 1.e-12)) + 0.2 - beta * beta - densityCorrection;

        const double fermiPlateauError = std::abs(numerator / denominator);

        energies.push_back(spState->GetKineticEnergy());
        fermiPlateauErrors.push_back(fermiPlateauError);
    }

    return TGraph{static_cast<Int_t>(history.size()), energies.data(), fermiPlateauErrors.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotFermiPlateauError(const bf::Detector &detector, const bf::Propagator &propagator)
{
    // Get the particles
    const auto &spMuon = bf::ParticleHelper::GetMuon();
    const auto &spChargedPion = bf::ParticleHelper::GetChargedPion();
    const auto &spChargedKaon = bf::ParticleHelper::GetChargedKaon();
    const auto &spProton = bf::ParticleHelper::GetProton();

    // Propagate the particles
    while ((spMuon->KineticEnergy() < 10000.) && !spMuon->HasFailed())
        propagator.PropagateBackwards(spMuon, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spChargedPion->KineticEnergy() < 10000.) && !spChargedPion->HasFailed())
        propagator.PropagateBackwards(spChargedPion, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spChargedKaon->KineticEnergy() < 10000.) && !spChargedKaon->HasFailed())
        propagator.PropagateBackwards(spChargedKaon, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spProton->KineticEnergy() < 10000.) && !spProton->HasFailed())
        propagator.PropagateBackwards(spProton, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    // Get the graphs
    auto muonGraph = bf::MultiGraphEntry{GetFermiPlateauErrorGraph(detector, propagator, spMuon), "\\mu", 0UL, true, 2};
    auto chargedPionGraph = bf::MultiGraphEntry{GetFermiPlateauErrorGraph(detector, propagator, spChargedPion), "\\pi^\\pm", 1UL, true, 2};
    auto chargedKaonGraph = bf::MultiGraphEntry{GetFermiPlateauErrorGraph(detector, propagator, spChargedKaon), "K^\\pm", 2UL, true, 2};
    auto protonGraph = bf::MultiGraphEntry{GetFermiPlateauErrorGraph(detector, propagator, spProton), "p\\", 3UL, true, 2};

    // Create the graphs
    bf::PlotOptions options;
    options.m_xAxisTitle = "T \\text{ (MeV)}";
    options.m_yAxisTitle = "R_F(T)\\";

    TCanvas *pCanvas = GetNewCanvas();

    TMultiGraph *pMultiGraph = nullptr;
    TLegend *pLegend = nullptr;

    auto graphs = std::vector<std::reference_wrapper<bf::MultiGraphEntry>>{muonGraph, chargedPionGraph, chargedKaonGraph, protonGraph};
    bf::PlotHelper::GetMultiGraph(graphs, options, pMultiGraph, pLegend);

    pMultiGraph->GetXaxis()->SetLimits(0., 10000.);

    pMultiGraph->DrawClone("A");
    pLegend->DrawClone();

    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TGraph GetResidualRangeVersusScaledTGraph(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    std::vector<double> scaledTs, residualRanges;
    const auto &        history = spParticle->GetHistory();

    if (history.empty())
        throw std::runtime_error{"Cannot plot empty particle history"};

    double boundary1 = 0.1;
    double boundary2 = 0.25;
    double boundary3 = 0.5;
    double boundary4 = 0.75;
    double boundary5 = 1.;

    bool reachedBoundary1 = false;
    bool reachedBoundary2 = false;
    bool reachedBoundary3 = false;
    bool reachedBoundary4 = false;
    bool reachedBoundary5 = false;

    for (const auto &spState : history)
    {
        const double scaledT = spState->GetKineticEnergy() / spParticle->Mass();

        if (!reachedBoundary1 && scaledT > boundary1) 
        {
            std::cout << "[m = " << spParticle->Mass() << "] reached T'=" << boundary1 << " at " << spState->GetResidualRange() << "cm" << std::endl;
            reachedBoundary1 = true;
        }

        if (!reachedBoundary2 && scaledT > boundary2) 
        {
            std::cout << "[m = " << spParticle->Mass() << "] reached T'=" << boundary2 << " at " << spState->GetResidualRange() << "cm" << std::endl;
            reachedBoundary2 = true;
        }

        if (!reachedBoundary3 && scaledT > boundary3) 
        {
            std::cout << "[m = " << spParticle->Mass() << "] reached T'=" << boundary3 << " at " << spState->GetResidualRange() << "cm" << std::endl;
            reachedBoundary3 = true;
        }

        if (!reachedBoundary4 && scaledT > boundary4) 
        {
            std::cout << "[m = " << spParticle->Mass() << "] reached T'=" << boundary4 << " at " << spState->GetResidualRange() << "cm" << std::endl;
            reachedBoundary4 = true;
        }

        if (!reachedBoundary5 && scaledT > boundary5) 
        {
            std::cout << "[m = " << spParticle->Mass() << "] reached T'=" << boundary5 << " at " << spState->GetResidualRange() << "cm" << std::endl;
            reachedBoundary5 = true;
        }

        scaledTs.push_back(spState->GetKineticEnergy() / spParticle->Mass());
        residualRanges.push_back(spState->GetResidualRange());
    }

    return TGraph{static_cast<Int_t>(scaledTs.size()), residualRanges.data(), scaledTs.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotResidualRangeVersusScaledT(const bf::Propagator &propagator)
{
    // Get the particles
    const auto &spMuon = bf::ParticleHelper::GetMuon();
    const auto &spChargedPion = bf::ParticleHelper::GetChargedPion();
    const auto &spChargedKaon = bf::ParticleHelper::GetChargedKaon();
    const auto &spProton = bf::ParticleHelper::GetProton();

    // Propagate the particles finely
    while ((spMuon->KineticEnergy() < 1000.) && !spMuon->HasFailed())
        propagator.PropagateBackwards(spMuon, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spChargedPion->KineticEnergy() < 1000.) && !spChargedPion->HasFailed())
        propagator.PropagateBackwards(spChargedPion, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spChargedKaon->KineticEnergy() < 1000.) && !spChargedKaon->HasFailed())
        propagator.PropagateBackwards(spChargedKaon, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spProton->KineticEnergy() < 1000.) && !spProton->HasFailed())
        propagator.PropagateBackwards(spProton, 0.03, bf::Propagator::PROPAGATION_MODE::MODAL);

    // Get the graphs
    auto muonFineGraph = bf::MultiGraphEntry{GetResidualRangeVersusScaledTGraph(propagator, spMuon), "\\mu", 0UL, true, 2};
    auto chargedPionFineGraph = bf::MultiGraphEntry{GetResidualRangeVersusScaledTGraph(propagator, spChargedPion), "\\pi^\\pm", 1UL, true, 2};
    auto chargedKaonFineGraph = bf::MultiGraphEntry{GetResidualRangeVersusScaledTGraph(propagator, spChargedKaon), "K^\\pm", 2UL, true, 2};
    auto protonFineGraph = bf::MultiGraphEntry{GetResidualRangeVersusScaledTGraph(propagator, spProton), "p\\", 3UL, true, 2};

    spMuon->Reset();
    spChargedPion->Reset();
    spChargedKaon->Reset();
    spProton->Reset();

    // Propagate the particles coarsely
    while ((spMuon->KineticEnergy() < 1000.) && !spMuon->HasFailed())
        propagator.PropagateBackwards(spMuon, 0.3, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spChargedPion->KineticEnergy() < 1000.) && !spChargedPion->HasFailed())
        propagator.PropagateBackwards(spChargedPion, 0.3, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spChargedKaon->KineticEnergy() < 1000.) && !spChargedKaon->HasFailed())
        propagator.PropagateBackwards(spChargedKaon, 0.3, bf::Propagator::PROPAGATION_MODE::MODAL);

    while ((spProton->KineticEnergy() < 1000.) && !spProton->HasFailed())
        propagator.PropagateBackwards(spProton, 0.3, bf::Propagator::PROPAGATION_MODE::MODAL);

    // Get the graphs
    auto muonCoarseGraph = bf::MultiGraphEntry{GetResidualRangeVersusScaledTGraph(propagator, spMuon), "", 0UL, true, 2};
    auto chargedPionCoarseGraph = bf::MultiGraphEntry{GetResidualRangeVersusScaledTGraph(propagator, spChargedPion), "", 1UL, true, 2};
    auto chargedKaonCoarseGraph = bf::MultiGraphEntry{GetResidualRangeVersusScaledTGraph(propagator, spChargedKaon), "", 2UL, true, 2};
    auto protonCoarseGraph = bf::MultiGraphEntry{GetResidualRangeVersusScaledTGraph(propagator, spProton), "", 3UL, true, 2};

    muonCoarseGraph.Graph().SetLineStyle(2);
    chargedPionCoarseGraph.Graph().SetLineStyle(2);
    chargedKaonCoarseGraph.Graph().SetLineStyle(2);
    protonCoarseGraph.Graph().SetLineStyle(2);

    // Plot the graphs
    bf::PlotOptions options;
    options.m_xAxisTitle = "\\text{Residual range (cm)}";
    options.m_yAxisTitle = "T\\text{ '}";
    options.m_legendX1 = 0.14;
    options.m_legendX2 = 0.22;

    TCanvas *pCanvas = GetNewCanvas();

    TMultiGraph *pMultiGraph = nullptr;
    TLegend *pLegend = nullptr;

    auto graphs = std::vector<std::reference_wrapper<bf::MultiGraphEntry>>{muonFineGraph, chargedPionFineGraph, chargedKaonFineGraph, protonFineGraph, muonCoarseGraph, chargedPionCoarseGraph, chargedKaonCoarseGraph, protonCoarseGraph};
    bf::PlotHelper::GetMultiGraph(graphs, options, pMultiGraph, pLegend);

    pMultiGraph->SetMaximum(1.4);
    pMultiGraph->GetXaxis()->SetLimits(0., 70.);

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

TGraph * GetLowEnergyApproxFractionalErrorGraph(const std::shared_ptr<bf::Particle> &spParticle, const std::function<double(double)> &minusdEdxGetter, const double range)
{
    std::vector<double> residualRanges, fractionalErrors;
    const auto &        history = spParticle->GetHistory();

    if (history.empty())
        throw std::runtime_error{"Cannot plot empty particle history"};

    for (const auto &spState : history)
    {
        const double residualRange = spState->GetResidualRange();

        if (residualRange < std::numeric_limits<double>::epsilon() || residualRange > range)
            continue;

        residualRanges.push_back(residualRange);

        const double trueMinusdEdx = -spState->GetdEdx();
        fractionalErrors.push_back((minusdEdxGetter(residualRange) - trueMinusdEdx) / trueMinusdEdx);
    }

    return new TGraph{static_cast<Int_t>(residualRanges.size()), residualRanges.data(), fractionalErrors.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas * PlotLowEnergyApproxError(const bf::Detector &detector, const bf::Propagator &propagator, const bf::QuickPidAlgorithm &quickPidAlg, const std::shared_ptr<bf::Particle> &spParticle, const unsigned int colour, const std::string &label, const double deltaX, const double range)
{    
    TCanvas *c = GetNewCanvas();
    
     TPad *pad1 = new TPad("pad1", "pad1", 0, 0.32, 1, 1.0);
    pad1->SetBottomMargin(0.02); // Upper and lower plot are joined
    //pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad

    // Get the modal dE/dx graph
    while ((spParticle->KineticEnergy() < 1000.) && !spParticle->HasFailed())
        propagator.PropagateBackwards(spParticle, deltaX, bf::Propagator::PROPAGATION_MODE::MODAL);

    auto modaldEdxGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spParticle, true), "Modal", colour, true, 2};

    // Create the graphs
    bf::PlotOptions options;
    options.m_xAxisTitle = "\\text{Residual range (cm)}";
    options.m_yAxisTitle = "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}";
    options.m_legendX1 = 0.65;
    options.m_legendX2 = 0.88;
    options.m_legendY1 = 0.67;
    options.m_legendY2 = 0.86;

    const double xi = 0.00291 * deltaX / 0.03; // MeV
    const double chi = std::log(2. * bf::PhysicalConstants::m_electronMass * xi / (detector.m_avgIonizationEnergy * detector.m_avgIonizationEnergy * 1.e-12)) + 0.2;

    // Draw the first order approx function
    TF1 firstOrderApproxFunc("firstOrderApproxFunc", "0.5 * sqrt([0] * [1] * [2] / ([3] * x))", 0., range);
    firstOrderApproxFunc.SetParameter(0, xi);
    firstOrderApproxFunc.SetParameter(1, chi);
    firstOrderApproxFunc.SetParameter(2, spParticle->Mass());
    firstOrderApproxFunc.SetParameter(3, deltaX);

    firstOrderApproxFunc.SetLineColor(bf::PlotHelper::GetSchemeColour(5UL));
    firstOrderApproxFunc.SetLineStyle(2);

    TMultiGraph *pMultiGraph = nullptr;
    TLegend *pLegend = nullptr;

    auto graphs = std::vector<std::reference_wrapper<bf::MultiGraphEntry>>{modaldEdxGraph};
    bf::PlotHelper::GetMultiGraph(graphs, options, pMultiGraph, pLegend);

    // Extra formatting
    pMultiGraph->GetXaxis()->SetLimits(0., range);
    pMultiGraph->SetMaximum(20.);
    pMultiGraph->SetMinimum(0.);
    // pMultiGraph->GetYaxis()->SetTitleOffset(0.4);
    pMultiGraph->GetYaxis()->SetTitle("-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}");
    pMultiGraph->GetYaxis()->SetTitleSize(0.047f);

    pMultiGraph->GetYaxis()->SetLabelSize(0.);
    pMultiGraph->GetXaxis()->SetLabelSize(0.);
    pMultiGraph->DrawClone("A");

    // Draw the second order approx function
    std::function<Double_t(const Double_t *, const Double_t *)> secondOrderApprox = [&](const Double_t *x, const Double_t *p) -> Double_t
    {
        return quickPidAlg.EstimatedEdx(x[0], spParticle->Mass(), deltaX);
    };

    TF1 secondOrderApproxFunc("secondOrderApproxFunc", secondOrderApprox, 0., range);

    secondOrderApproxFunc.SetLineColor(bf::PlotHelper::GetSchemeColour(6UL));
    secondOrderApproxFunc.SetLineStyle(2);

    // Draw the label
    TLatex latex;
    latex.SetTextSize(0.1);
    latex.DrawLatex(8.8, 12.5, label.c_str());

    // Do not draw the Y axis label on the upper plot and redraw a small
   // axis instead, in order to avoid the first label (0) to be clipped.
   
   TGaxis *axis = new TGaxis( 0., 0., 0., 20., 0., 20., 510,"");
   axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   axis->SetLabelSize(18);
   axis->Draw();

   firstOrderApproxFunc.DrawClone("same");
    secondOrderApproxFunc.DrawClone("same");

    pLegend->AddEntry(&firstOrderApproxFunc, "First order approx.", "l");
    pLegend->AddEntry(&secondOrderApproxFunc, "Second order approx.", "l");
    pLegend->DrawClone();

    ///---------------------

    // PAD 2

    c->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1, 0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.3);
    //pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad

    TGraph *pSecondOrderFractionalErrorGraph = GetLowEnergyApproxFractionalErrorGraph(spParticle, 
    [&](double residualRange)
    {
        return quickPidAlg.EstimatedEdx(residualRange, spParticle->Mass(), deltaX);
    }, range);

    TGraph *pFirstOrderFractionalErrorGraph = GetLowEnergyApproxFractionalErrorGraph(spParticle, 
    [&](double residualRange)
    {
        return firstOrderApproxFunc.Eval(residualRange);
    }, range);

    pFirstOrderFractionalErrorGraph->SetLineColor(bf::PlotHelper::GetSchemeColour(5UL));
    pSecondOrderFractionalErrorGraph->SetLineColor(bf::PlotHelper::GetSchemeColour(6UL));

    pFirstOrderFractionalErrorGraph->GetXaxis()->SetLimits(0., range);
    pFirstOrderFractionalErrorGraph->SetMaximum(0.2);
    pFirstOrderFractionalErrorGraph->SetMinimum(-0.4);
    
    pFirstOrderFractionalErrorGraph->SetLineStyle(1);
    pFirstOrderFractionalErrorGraph->SetLineWidth(2);
    pFirstOrderFractionalErrorGraph->GetXaxis()->SetLabelSize(0.11f);
    pFirstOrderFractionalErrorGraph->GetXaxis()->SetTitleSize(0.11f);
    pFirstOrderFractionalErrorGraph->GetYaxis()->SetNdivisions(6, 4, 0);
    pFirstOrderFractionalErrorGraph->GetYaxis()->SetTitle("Frac. error");
    pFirstOrderFractionalErrorGraph->GetYaxis()->SetTitleOffset(1.1);
    pFirstOrderFractionalErrorGraph->GetXaxis()->SetTitle("\\text{Residual range (cm)}");
    pFirstOrderFractionalErrorGraph->GetYaxis()->SetTitleSize(20);
    pFirstOrderFractionalErrorGraph->GetYaxis()->SetTitleFont(43);
    pFirstOrderFractionalErrorGraph->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    pFirstOrderFractionalErrorGraph->GetYaxis()->SetLabelSize(0);

    pFirstOrderFractionalErrorGraph->Draw("AL");

    TGaxis *firstAxis = new TGaxis( 0., -0.4, 0., 0.2, -0.4, 0.2, 4,"");
    firstAxis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    firstAxis->SetLabelSize(18);
    firstAxis->Draw();
    
    pSecondOrderFractionalErrorGraph->SetLineStyle(1);
    pSecondOrderFractionalErrorGraph->SetLineWidth(2);

    pSecondOrderFractionalErrorGraph->Draw("L same");

    return c;
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