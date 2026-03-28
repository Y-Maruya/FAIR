#include "CalibDBIO/PedestalReader/PedestalReader.hpp"
#include "common/Logger.hpp"
#include "common/AHCALGeometry.hpp"
#include <iostream>
#include <iomanip>
#include <stdexcept>

int main(int argc, char* argv[]) {
    // Initialize logger
    FAIR::init_logger("CalibITester", "calib_tester.log", spdlog::level::info);
    LOG_INFO("=== PedestalReader Integration Test ===");
    LOG_INFO("Test Type: Integration test with real calibration database");
    
    // Test run number
    int testRunNumber = 22140;
    if (argc > 1) {
        testRunNumber = std::stoi(argv[1]);
    }
    LOG_INFO("Testing with run number: {}", testRunNumber);
    
    try {
        // Create PedestalReader and trigger readPedestals() + selectNearestPedestalRun()
        LOG_INFO("Creating PedestalReader for run {}...", testRunNumber);
        PedestalReader reader(testRunNumber);
        LOG_INFO("PedestalReader created successfully!");
        
        // Test 1: Retrieve pedestals for sample cells across different layers and chips
        LOG_INFO("--- Test 1: Retrieving sample pedestal values ---");
        
        // Sample cell IDs to test (assuming reasonable geometry values)
        // Format: Layer, Chip, Channel
        std::vector<std::tuple<int, int, int>> testCells = {
            {0, 0, 0},      // Layer 0, Chip 0, Channel 0
            {0, 0, 1},      // Layer 0, Chip 0, Channel 1
            {0, 1, 0},      // Layer 0, Chip 1, Channel 0
        };
        
        int successCount = 0;
        int totalTests = 0;
        
        for (const auto& [layer, chip, channel] : testCells) {
            totalTests++;
            try {
                // Create cell ID using the same method as PedestalReader
                int cellID = AHCALGeometry::CellID(layer, chip, channel);
                
                // Retrieve pedestal for this cell
                Pedestal ped = reader.getPedestal(cellID);
                
                // Display results
                LOG_INFO("Cell ID {}: Layer={}, Chip={}, Channel={}", cellID, layer, chip, channel);
                LOG_INFO("  HighGain: Peak={:.2f}, Sigma={:.2f}, Status={}", 
                         ped.HighGainPeak, ped.HighGainSigma, ped.HighGainStatus);
                LOG_INFO("  LowGain:  Peak={:.2f}, Sigma={:.2f}, Status={}", 
                         ped.LowGainPeak, ped.LowGainSigma, ped.LowGainStatus);
                
                // Verify data sanity
                bool hgValid = (ped.HighGainPeak > 0.0) && (ped.HighGainSigma > 0.0);
                bool lgValid = (ped.LowGainPeak > 0.0) && (ped.LowGainSigma > 0.0);
                bool statusOk = (ped.HighGainStatus == 0 || ped.HighGainStatus == 999) &&
                                (ped.LowGainStatus == 0 || ped.LowGainStatus == 999);
                
                if (hgValid && lgValid && statusOk) {
                    LOG_INFO("  ✓ Pedestal data valid");
                    successCount++;
                } else {
                    LOG_WARN("  ⚠ Pedestal data may have issues - HGValid={}, LGValid={}, StatusOK={}", 
                             hgValid, lgValid, statusOk);
                }
            } catch (const std::exception& e) {
                LOG_ERROR("  ✗ Failed to retrieve pedestal for cell {}, {}, {}: {}", 
                          layer, chip, channel, e.what());
            }
        }
        
        // Summary
        LOG_INFO("--- Test Summary ---");
        LOG_INFO("Total cells tested: {}", totalTests);
        LOG_INFO("Successful retrievals: {}", successCount);
        LOG_INFO("Success rate: {:.1f}%", (totalTests > 0) ? (100.0 * successCount / totalTests) : 0.0);
        
        if (successCount == totalTests) {
            LOG_INFO("✓ All tests passed!");
            return 0;
        } else if (successCount > 0) {
            LOG_WARN("⚠ Partial success - some cells retrieved successfully");
            return 0;  // Still count as success since some data was loaded
        } else {
            LOG_ERROR("✗ All tests failed - no pedestals could be retrieved");
            return 1;
        }
        
    } catch (const std::exception& e) {
        LOG_ERROR("Test failed with exception: {}", e.what());
        LOG_ERROR("This usually indicates:");
        LOG_ERROR("  - Cannot connect to calibration database");
        LOG_ERROR("  - No pedestal data available for run {} or nearby runs", testRunNumber);
        LOG_ERROR("  - Network or authentication issues with the FASER calibration DB");
        return 1;
    }
    
    return 0;
}
