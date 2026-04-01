#ifndef ANALYZE_TRIGGER_GROUPINGS_H
#define ANALYZE_TRIGGER_GROUPINGS_H

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cctype>

namespace ReferenceData {
    // Define first reference dataset
    const std::vector<double> referencePTGamma = {3.36, 4.39, 5.41, 6.42, 7.43, 8.44, 9.80, 11.83, 14.48};
    const std::vector<double> referenceRatio = {0.594, 0.664, 0.626, 0.658, 0.900, 0.715, 0.872, 0.907, 0.802};
    const std::vector<double> referenceStatError = {0.014, 0.028, 0.043, 0.061, 0.113, 0.130, 0.120, 0.190, 0.290};

    // Define second reference data (from the second table - Reference Two)
    const std::vector<double> referenceTwoPTGamma = {3.34, 4.38, 5.40, 6.41, 7.42, 8.43, 9.78, 11.81, 14.41};
    const std::vector<double> referenceTwoRatio = {0.477, 0.455, 0.448, 0.430, 0.338, 0.351, 0.400, 0.286, 0.371};
    const std::vector<double> referenceTwoStatError = {0.0020, 0.0060, 0.012, 0.021, 0.032, 0.053, 0.070, 0.130, 0.180};
}

// Define CutValues, FitParameters, and HistogramData within a dedicated namespace
namespace DataStructures {

    struct RunInfo {
        std::vector<int> runsBeforeFirmwareUpdate;
        std::vector<int> runsAfterFirmwareUpdate;
    };

    struct CutValues {
        float clusEnergy = 0;
        float asymmetry = 0;
        float chi = 0;
        std::string triggerName;
        float pTMin = -1;  // Default to -1 indicating no pT bin
        float pTMax = -1;  // Default to -1 indicating no pT bin
    };

    struct HistogramData {
        CutValues cuts;
        std::string histName;  // Name of the histogram

        // Fitted parameters and their errors
        double meanPi0;
        double meanPi0Error;
        double sigmaPi0;
        double sigmaPi0Error;
        double meanEta;
        double meanEtaError;
        double sigmaEta;
        double sigmaEtaError;

        // Mass ratio and its error
        double massRatio;
        double massRatioError;

        // Signal and background yields for pi0
        double signalPi0Yield;
        double signalPi0Error;
        double backgroundPi0Yield;
        double backgroundPi0Error;
        double signalToBackgroundPi0Ratio;
        double signalToBackgroundPi0Error;

        // Signal and background yields for eta
        double signalEtaYield;
        double signalEtaError;
        double backgroundEtaYield;
        double backgroundEtaError;
        double signalToBackgroundEtaRatio;
        double signalToBackgroundEtaError;

        // Resolution parameters
        double pi0FitResolution;
        double pi0FitResolutionError;
        double etaFitResolution;
        double etaFitResolutionError;
    };

    struct FitParameters {
        // Common parameters
        double amplitudeEstimate;
        double amplitudeMin;
        double amplitudeMax;

        // Sigmoid function parameters
        double slopeEstimate;
        double slopeMin;
        double slopeMax;
        double xOffsetEstimate;
        double xOffsetMin;
        double xOffsetMax;

        // Error function parameters
        double sigmaEstimate;
        double sigmaMin;
        double sigmaMax;
    };

    struct TotalPhotonLog {
        std::string triggerGroupName;
        std::string triggerName;
        float clusEnergy = 0.f;
        float chi       = 0.f;
        float asymmetry = 0.f;
        float pTMin     = -1.f;
        float pTMax     = -1.f;
        int totalEntries = 0;
        std::string massWindowLabel;
    };

    struct PtWeightingLog {
        std::string triggerGroupName;
        std::string triggerName;
        float clusEnergy   = 0.f;
        float chi         = 0.f;
        float asymmetry   = 0.f;
        float pTMin       = -1.f;
        float pTMax       = -1.f;
        double weightedAveragePt = 0.0;
        std::string massWindowLabel;
    };

    struct IsolationData {
        std::string triggerGroupName;
        std::string triggerName;
        
        float clusEnergy = 0.f;
        float chi       = 0.f;
        float asymmetry = 0.f;
        float pTMin     = -1.f;
        float pTMax     = -1.f;
        float isoMin    = 0.f;
        float isoMax    = 0.f;
        
        int isolatedCounts = 0;   // # of isolated photons in this bin
        int totalCounts    = 0;
        double ratio  = 0.0;
        double error  = 0.0;

        double weightedPt         = 0.0;
        double binWidth           = 0.0;
        double binCenter          = 0.0;
        double isolatedYield      = 0.0;
        double isolatedYieldError = 0.0;

        std::string massWindowLabel;    // "inMassWindow" or "outsideMassWindow"
        std::string showerCutLabel;     // "withShowerShapeCuts" or "withoutShowerShapeCuts"
    };


    // ------------------------------------
    // The “with Pt” version used when grouping
    // ------------------------------------
    struct IsolationDataWithPt {
        float ptMin          = -1.f;
        float ptMax          = -1.f;
        double weightedPt    = 0.0;
        double ratio         = 0.0;
        double error         = 0.0;
        double isolatedYield = 0.0;
        double isolatedYieldError = 0.0;
        float isoMin         = 0.f;
        float isoMax         = 0.f;
        std::string triggerName;

        // The full "raw" data from CSV (includes showerCutLabel, counts, etc.)
        IsolationData isoData;
    };

    struct CutCombinationData {
        double clusEnergy;
        double chi;
        double asymmetry;

        std::vector<double> pTCentersPi0;
        std::vector<double> meanPi0Values;
        std::vector<double> meanPi0Errors;
        std::vector<double> sigmaPi0Values;
        std::vector<double> sigmaPi0Errors;
        std::vector<double> resolutionPi0Values;
        std::vector<double> resolutionPi0Errors;
        std::vector<double> signalToBackgroundPi0Ratios;
        std::vector<double> signalToBackgroundPi0Errors;
        std::vector<double> pi0YieldValues;
        std::vector<double> pi0YieldErrors;
        std::vector<std::string> triggersUsedPi0;

        std::vector<double> pTCentersEta;
        std::vector<double> meanEtaValues;
        std::vector<double> meanEtaErrors;
        std::vector<double> sigmaEtaValues;
        std::vector<double> sigmaEtaErrors;
        std::vector<double> resolutionEtaValues;
        std::vector<double> resolutionEtaErrors;
        std::vector<double> signalToBackgroundEtaRatios;
        std::vector<double> signalToBackgroundEtaErrors;
        std::vector<double> etaYieldValues;
        std::vector<double> etaYieldErrors;
        std::vector<std::string> triggersUsedEta;

        std::set<std::string> triggersInData;

        // For overlay plots
        std::map<std::string, std::vector<double>> triggerToPtCentersPi0;
        std::map<std::string, std::vector<double>> triggerToMeanPi0Values;
        std::map<std::string, std::vector<double>> triggerToMeanPi0Errors;
        std::map<std::string, std::vector<double>> triggerToSigmaPi0Values;
        std::map<std::string, std::vector<double>> triggerToSigmaPi0Errors;
        std::map<std::string, std::vector<double>> triggerToResolutionPi0Values;
        std::map<std::string, std::vector<double>> triggerToResolutionPi0Errors;
        std::map<std::string, std::vector<double>> triggerToSignalToBackgroundPi0Ratios;
        std::map<std::string, std::vector<double>> triggerToSignalToBackgroundPi0Errors;
        std::map<std::string, std::vector<double>> triggerToPi0YieldValues;
        std::map<std::string, std::vector<double>> triggerToPi0YieldErrors;


        std::map<std::string, std::vector<double>> triggerToPtCentersEta;
        std::map<std::string, std::vector<double>> triggerToMeanEtaValues;
        std::map<std::string, std::vector<double>> triggerToMeanEtaErrors;
        std::map<std::string, std::vector<double>> triggerToSigmaEtaValues;
        std::map<std::string, std::vector<double>> triggerToSigmaEtaErrors;
        std::map<std::string, std::vector<double>> triggerToResolutionEtaValues;
        std::map<std::string, std::vector<double>> triggerToResolutionEtaErrors;
        std::map<std::string, std::vector<double>> triggerToSignalToBackgroundEtaRatios;
        std::map<std::string, std::vector<double>> triggerToSignalToBackgroundEtaErrors;
        std::map<std::string, std::vector<double>> triggerToEtaYieldValues;
        std::map<std::string, std::vector<double>> triggerToEtaYieldErrors;
    };

    std::vector<std::pair<double, double>> pT_bins = {
        {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, {5.0, 6.0},
        {6.0, 7.0}, {7.0, 8.0}, {8.0, 9.0}, {9.0, 10.0},
        {10.0, 12.0}, {12.0, 15.0}, {15.0, 20.0}, {20.0, 30.0}
    };

} // namespace DataStructures


// Namespace for trigger configurations
namespace TriggerConfig {
    // List of all triggers we're interested in
    const std::vector<std::string> allTriggers = {
        "MBD_NandS_geq_1",
        "Photon_2_GeV_plus_MBD_NS_geq_1",
        "Photon_3_GeV_plus_MBD_NS_geq_1",
        "Photon_4_GeV_plus_MBD_NS_geq_1",
        "Photon_5_GeV_plus_MBD_NS_geq_1",
        "Jet_8_GeV_plus_MBD_NS_geq_1",
        "Jet_10_GeV_plus_MBD_NS_geq_1",
        "Jet_12_GeV_plus_MBD_NS_geq_1"

    };

    // List of Photon triggers (excluding MBD_NandS_geq_1)
    const std::vector<std::string> photonTriggers = {
        "Photon_2_GeV_plus_MBD_NS_geq_1",
        "Photon_3_GeV_plus_MBD_NS_geq_1",
        "Photon_4_GeV_plus_MBD_NS_geq_1",
        "Photon_5_GeV_plus_MBD_NS_geq_1"
    };

    // Map of triggers to colors for plotting
    const std::map<std::string, int> triggerColorMap = {
        {"MBD_NandS_geq_1", kBlack},
        {"Photon_2_GeV_plus_MBD_NS_geq_1", kRed},
        {"Photon_3_GeV_plus_MBD_NS_geq_1", kBlue},
        {"Photon_4_GeV_plus_MBD_NS_geq_1", kGreen + 2},
        {"Photon_5_GeV_plus_MBD_NS_geq_1", kMagenta},
        
        {"Jet_8_GeV_plus_MBD_NS_geq_1", kOrange - 3},
        {"Jet_10_GeV_plus_MBD_NS_geq_1", kAzure - 4},
        {"Jet_12_GeV_plus_MBD_NS_geq_1", kPink + 6}
    };

    // Map of triggers to human-readable names
    const std::map<std::string, std::string> triggerNameMap = {
        {"MBD_NandS_geq_1", "MBD NS #geq 1"},
        {"Photon_2_GeV_plus_MBD_NS_geq_1", "Photon 2 GeV + MBD NS #geq 1"},
        {"Photon_3_GeV_plus_MBD_NS_geq_1", "Photon 3 GeV + MBD NS #geq 1"},
        {"Photon_4_GeV_plus_MBD_NS_geq_1", "Photon 4 GeV + MBD NS #geq 1"},
        {"Photon_5_GeV_plus_MBD_NS_geq_1", "Photon 5 GeV + MBD NS #geq 1"},
        // Jet triggers:
        {"Jet_8_GeV_plus_MBD_NS_geq_1",  "Jet 8 GeV + MBD NS #geq 1"},
        {"Jet_10_GeV_plus_MBD_NS_geq_1", "Jet 10 GeV + MBD NS #geq 1"},
        {"Jet_12_GeV_plus_MBD_NS_geq_1", "Jet 12 GeV + MBD NS #geq 1"}
    };

    const std::map<std::pair<float, float>, int> isoEtRangeColorMap = {
        {{-100, 6}, kRed + 1},
        {{-100, 10}, kRed + 1},
        {{-10, 0}, kGreen + 2},
        {{0, 10}, kMagenta + 2}
    };

    const std::map<std::pair<std::string, std::string>, DataStructures::FitParameters> triggerFitParameters = {
        { {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Photon_2_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            // Sigmoid function parameters
            0.68,  // slopeEstimate
            0.67,  // slopeMin
            0.69,  // slopeMax
            
            4.73,  // xOffsetEstimate
            4.58,  // xOffsetMin
            4.78,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Photon_3_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.54,  // slopeEstimate
            0.53,  // slopeMin
            0.55,  // slopeMax
            
            7.03,  // xOffsetEstimate
            6.93,  // xOffsetMin
            7.13,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Photon_4_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.42,  // slopeEstimate
            0.415,  // slopeMin
            0.425,  // slopeMax
            
            8.45,   // xOffsetEstimate
            8.35,   // xOffsetMin
            8.55,    // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1_beforeTriggerFirmwareUpdate", "Photon_3_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.52,  // slopeEstimate
            0.51,  // slopeMin
            0.53,  // slopeMax
            
            7.25,  // xOffsetEstimate
            7.2,  // xOffsetMin
            7.3,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1_beforeTriggerFirmwareUpdate", "Photon_4_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.47,  // slopeEstimate
            0.46,  // slopeMin
            0.48,  // slopeMax
            
            8.45,   // xOffsetEstimate
            8.35,   // xOffsetMin
            8.55,    // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1_beforeTriggerFirmwareUpdate", "Photon_5_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.45,  // slopeEstimate
            0.44,  // slopeMin
            0.46,  // slopeMax
            
            9.65,   // xOffsetEstimate
            9.64,   // xOffsetMin
            9.66,    // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate", "Photon_3_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.64,  // slopeEstimate
            0.63,  // slopeMin
            0.65,  // slopeMax
            
            7.85,  // xOffsetEstimate
            7.8,  // xOffsetMin
            7.9,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate", "Photon_4_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.57,  // slopeEstimate
            0.52,  // slopeMin
            0.62,  // slopeMax
            
            9.0,   // xOffsetEstimate
            9.1,   // xOffsetMin
            8.9,    // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate", "Photon_5_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.5,  // slopeEstimate
            0.49,  // slopeMin
            0.51,  // slopeMax
            
            10.1,   // xOffsetEstimate
            10,   // xOffsetMin
            10.2,    // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_beforeTriggerFirmwareUpdate", "Photon_3_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.52,  // slopeEstimate
            0.51,  // slopeMin
            0.53,  // slopeMax
            
            7.25,  // xOffsetEstimate
            7.2,  // xOffsetMin
            7.3,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_beforeTriggerFirmwareUpdate", "Photon_4_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.47,  // slopeEstimate
            0.46,  // slopeMin
            0.48,  // slopeMax
            
            8.45,   // xOffsetEstimate
            8.35,   // xOffsetMin
            8.55,    // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate", "Photon_3_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.54,  // slopeEstimate
            0.53,  // slopeMin
            0.55,  // slopeMax
            
            8.0,  // xOffsetEstimate
            7.9,  // xOffsetMin
            7.13,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate", "Photon_4_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.42,  // slopeEstimate
            0.415,  // slopeMin
            0.425,  // slopeMax
            
            8.45,   // xOffsetEstimate
            8.35,   // xOffsetMin
            8.55,    // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Photon_2_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.42,  // slopeEstimate
            0.415,  // slopeMin
            0.425,  // slopeMax
            
            8.45,   // xOffsetEstimate
            8.35,   // xOffsetMin
            8.55,    // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Photon_4_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            
            // Sigmoid function parameters
            0.42,  // slopeEstimate
            0.415,  // slopeMin
            0.425,  // slopeMax
            
            8.45,   // xOffsetEstimate
            8.35,   // xOffsetMin
            8.55,    // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"", "Photon_2_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.25,  // amplitudeEstimate
            1.0,   // amplitudeMin
            1.3,   // amplitudeMax
            // Sigmoid function parameters
            0.72,  // slopeEstimate
            0.71,  // slopeMin
            0.73,  // slopeMax
            4.65,  // xOffsetEstimate
            4.55,  // xOffsetMin
            4.75,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"", "Photon_3_GeV_plus_MBD_NS_geq_1"}, {
            1.25,   // amplitudeEstimate
            0.64,   // slopeEstimate
            7.2,   // xOffsetEstimate
            1.0,  // amplitudeMin
            1.3,  // amplitudeMax
            0.62,   // slopeMin
            0.66,  // slopeMax
            7.0,   // xOffsetMin
            7.5    // xOffsetMax
        } },
        { {"", "Photon_4_GeV_plus_MBD_NS_geq_1"}, {
            1.25,   // amplitudeEstimate
            0.55,   // slopeEstimate
            8.5,   // xOffsetEstimate
            1.0,  // amplitudeMin
            1.3,  // amplitudeMax
            0.53,   // slopeMin
            0.57,  // slopeMax
            8.4,   // xOffsetMin
            8.6    // xOffsetMax
        } },
        
        { {"MBD_NandS_geq_1_Jet_8_GeV_plus_MBD_NS_geq_1_Jet_10_GeV_plus_MBD_NS_geq_1_Jet_12_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate", "Jet_8_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.00,  // amplitudeEstimate
            0.999,   // amplitudeMin
            1.001,   // amplitudeMax
            
            // Sigmoid function parameters
            0.52,  // slopeEstimate
            0.51,  // slopeMin
            0.53,  // slopeMax
            
            10.0,  // xOffsetEstimate
            9.9,  // xOffsetMin
            10.1,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Jet_8_GeV_plus_MBD_NS_geq_1_Jet_10_GeV_plus_MBD_NS_geq_1_Jet_12_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate", "Jet_10_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.00,  // amplitudeEstimate
            0.999,   // amplitudeMin
            1.001,   // amplitudeMax
            
            // Sigmoid function parameters
            0.52,  // slopeEstimate
            0.51,  // slopeMin
            0.53,  // slopeMax
            
            12.0,  // xOffsetEstimate
            11.9,  // xOffsetMin
            12.1,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Jet_8_GeV_plus_MBD_NS_geq_1_Jet_10_GeV_plus_MBD_NS_geq_1_Jet_12_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate", "Jet_12_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.00,  // amplitudeEstimate
            0.999,   // amplitudeMin
            1.001,   // amplitudeMax
            
            // Sigmoid function parameters
            0.52,  // slopeEstimate
            0.51,  // slopeMin
            0.53,  // slopeMax
            
            14.3,  // xOffsetEstimate
            14.2,  // xOffsetMin
            14.4,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Jet_8_GeV_plus_MBD_NS_geq_1_Jet_10_GeV_plus_MBD_NS_geq_1", "Jet_8_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.00,  // amplitudeEstimate
            0.999,   // amplitudeMin
            1.001,   // amplitudeMax
            
            // Sigmoid function parameters
            0.52,  // slopeEstimate
            0.51,  // slopeMin
            0.53,  // slopeMax
            
            10.0,  // xOffsetEstimate
            9.9,  // xOffsetMin
            10.1,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        { {"MBD_NandS_geq_1_Jet_8_GeV_plus_MBD_NS_geq_1_Jet_10_GeV_plus_MBD_NS_geq_1", "Jet_10_GeV_plus_MBD_NS_geq_1"}, {
            // Common parameters
            1.00,  // amplitudeEstimate
            0.999,   // amplitudeMin
            1.001,   // amplitudeMax
            
            // Sigmoid function parameters
            0.52,  // slopeEstimate
            0.51,  // slopeMin
            0.53,  // slopeMax
            
            12.0,  // xOffsetEstimate
            11.9,  // xOffsetMin
            12.1,  // xOffsetMax

            // Error function parameters
            0.5,   // sigmaEstimate (placeholder value)
            0.1,   // sigmaMin (placeholder value)
            1.0    // sigmaMax (placeholder value)
        } },
        
        { {"", "Photon_5_GeV_plus_MBD_NS_geq_1"}, {
            1.0,   // amplitudeEstimate
            0.48,  // slopeEstimate
            10.5,   // xOffsetEstimate
            0.98,  // amplitudeMin
            1.05,  // amplitudeMax
            0.46,   // slopeMin
            0.5,   // slopeMax
            10.4,   // xOffsetMin
            10.6    // xOffsetMax
        } }
    };

    // Define a map from trigger names to photon thresholds
    std::map<std::string, double> triggerThresholds = {
        {"MBD_NandS_geq_1", 0.0},
        {"Photon_2_GeV_plus_MBD_NS_geq_1", 2.0},
        {"Photon_3_GeV_plus_MBD_NS_geq_1", 3.0},
        {"Photon_4_GeV_plus_MBD_NS_geq_1", 4.0},
        {"Photon_5_GeV_plus_MBD_NS_geq_1", 5.0},
        // Add other triggers if necessary
    };

}


// Namespace for trigger combination names
namespace TriggerCombinationNames {
    const std::map<std::string, std::string> triggerCombinationNameMap = {
        {"MBD_NandS_geq_1", "MBD NS #geq 1"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 2 GeV"},
        {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 3 GeV"},
        {"MBD_NandS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 4 GeV"},
        {"MBD_NandS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 5 GeV"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 2, 3 GeV"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 2, 4 GeV"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 2, 5 GeV"},
        {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 3, 4 GeV"},
        {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 3, 5 GeV"},
        {"MBD_NandS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 4, 5 GeV"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 2, 3, 4 GeV"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 2, 3, 4 GeV"},
        {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "MBD NS #geq 1 and Photon 3, 4, 5 GeV"},
    };
}

namespace Utils {
    // Helper function to normalize the trigger combination string for case-insensitive and whitespace-insensitive comparison
    std::string normalizeString(const std::string& str) {
        std::string normalized = str;
        std::transform(normalized.begin(), normalized.end(), normalized.begin(), ::tolower);
        normalized.erase(std::remove_if(normalized.begin(), normalized.end(), ::isspace), normalized.end());
        return normalized;
    }
    // Function to check if a string ends with another string
    bool EndsWith(const std::string& fullString, const std::string& ending) {
        if (fullString.length() >= ending.length()) {
            return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
        } else {
            return false;
        }
    }
    // Helper function to strip firmware update tags from the combination name
    std::string stripFirmwareTag(const std::string& combinationName) {
        std::string strippedName = combinationName;
        const std::vector<std::string> firmwareTags = {
            "_beforeTriggerFirmwareUpdate",
            "_afterTriggerFirmwareUpdate"
        };

        for (const auto& tag : firmwareTags) {
            size_t pos = strippedName.find(tag);
            if (pos != std::string::npos) {
                strippedName.erase(pos, tag.length());
                break;
            }
        }
        return strippedName;
    }

    std::string getTriggerCombinationName(const std::string& combinationName, const std::map<std::string, std::string>& nameMap) {
        std::string strippedCombinationName = stripFirmwareTag(combinationName);
        std::string normalizedCombinationName = normalizeString(strippedCombinationName);

        for (const auto& entry : nameMap) {
            if (normalizeString(entry.first) == normalizedCombinationName) {
                // Append firmware tag back to the human-readable name if present
                if (EndsWith(combinationName, "_beforeTriggerFirmwareUpdate")) {
                    return entry.second + " (Before Firmware Update)";
                } else if (EndsWith(combinationName, "_afterTriggerFirmwareUpdate")) {
                    return entry.second + " (After Firmware Update)";
                } else {
                    return entry.second;
                }
            }
        }
        return combinationName; // Default to combination name if not found
    }

    // Function to format a double to three significant figures as a string
    std::string formatToThreeSigFigs(double value) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(3) << value;
        return out.str();
    }

    // Existing sigmoidFit function remains unchanged
    TF1* sigmoidFit(const std::string& name, double xmin, double xmax,
                    double amplitude, double slope, double xOffset,
                    double amplitudeMin, double amplitudeMax,
                    double slopeMin, double slopeMax,
                    double xOffsetMin, double xOffsetMax) {
        // Define a sigmoid function for fitting
        TF1* fitFunc = new TF1(name.c_str(), "[0]/(1+exp(-[1]*(x-[2])))", xmin, xmax);
        fitFunc->SetParNames("Amplitude", "Slope", "XOffset");

        // Set initial parameters
        fitFunc->SetParameter(0, amplitude);  // Amplitude
        fitFunc->SetParameter(1, slope);      // Slope
        fitFunc->SetParameter(2, xOffset);    // XOffset

        // Set parameter limits
        fitFunc->SetParLimits(0, amplitudeMin, amplitudeMax);  // Amplitude limits
        fitFunc->SetParLimits(1, slopeMin, slopeMax);          // Slope limits
        fitFunc->SetParLimits(2, xOffsetMin, xOffsetMax);      // XOffset limits

        return fitFunc;
    }


    TF1* erfFit(const std::string& name, double xmin, double xmax,
                double amplitude, double xOffset, double sigma,
                double amplitudeMin, double amplitudeMax,
                double xOffsetMin, double xOffsetMax,
                double sigmaMin, double sigmaMax) {
        // Define an error function for fitting
        TF1* fitFunc = new TF1(name.c_str(), "[0]*0.5*(1+TMath::Erf((x-[1])/(sqrt(2)*[2])))", xmin, xmax);
        fitFunc->SetParNames("Amplitude", "XOffset", "Sigma");

        // Set initial parameters
        fitFunc->SetParameter(0, amplitude);  // Amplitude
        fitFunc->SetParameter(1, xOffset);    // XOffset
        fitFunc->SetParameter(2, sigma);      // Sigma

        // Set parameter limits
        fitFunc->SetParLimits(0, amplitudeMin, amplitudeMax);  // Amplitude limits
        fitFunc->SetParLimits(1, xOffsetMin, xOffsetMax);      // XOffset limits
        fitFunc->SetParLimits(2, sigmaMin, sigmaMax);          // Sigma limits

        return fitFunc;
    }

}

#endif // ANALYZE_TRIGGER_GROUPINGS_H
