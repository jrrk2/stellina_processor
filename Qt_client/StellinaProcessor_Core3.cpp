#include "StellinaProcessor.h"
#include "CoordinateUtils.h"
#include <QApplication>
#include <QDir>
#include <QFileInfo>
#include <QJsonObject>
#include <QJsonDocument>
#include <QJsonArray>
#include <QDateTime>
#include <QSettings>
#include <QSplitter>
#include <QMenuBar>
#include <QStatusBar>
#include <QHeaderView>
#include <QElapsedTimer>
#include <QProcess>
#include <QStandardItemModel>
#include <QTableWidgetItem>
#include <QThread>
#include <QPainter>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>  // for std::fmod if needed

bool StellinaProcessor::extractCoordinates(const QJsonObject &json, double &alt, double &az) {
    if (json.contains("motors")) {
        QJsonObject motors = json["motors"].toObject();
        if (motors.contains("ALT") && motors.contains("AZ")) {
            alt = motors["ALT"].toDouble();
            az = motors["AZ"].toDouble();
            return true;
        }
    }
    
    if (json.contains("altitude") && json.contains("azimuth")) {
        alt = json["altitude"].toDouble();
        az = json["azimuth"].toDouble();
        return true;
    }
    
    if (json.contains("alt") && json.contains("az")) {
        alt = json["alt"].toDouble();
        az = json["az"].toDouble();
        return true;
    }
    
    return false;
}

bool StellinaProcessor::checkStellinaQuality(const QJsonObject &json) {
    if (json.contains("stackingData")) {
        QJsonObject stackingData = json["stackingData"].toObject();
        
        if (stackingData.contains("error") && !stackingData["error"].isNull()) {
            return false;
        }
        
        if (stackingData.contains("liveRegistrationResult")) {
            QJsonObject regResult = stackingData["liveRegistrationResult"].toObject();
            
            if (regResult.contains("status")) {
                int status = regResult["status"].toInt();
                if (status != 0) {
                    return false;
                }
            }
            
            if (regResult.contains("statusMessage")) {
                QString statusMsg = regResult["statusMessage"].toString();
                if (statusMsg != "StackingOk") {
                    return false;
                }
            }
            
            if (regResult.contains("starsUsed")) {
                int starsUsed = regResult["starsUsed"].toInt();
                if (starsUsed < 10) {
                    return false;
                }
            }
        }
    }
    
    if (json.contains("quality")) {
        return json["quality"].toBool();
    }
    
    if (json.contains("used_for_stacking")) {
        return json["used_for_stacking"].toBool();
    }
    
    return true;
}

QString StellinaProcessor::formatProcessingTime(qint64 milliseconds) {
    qint64 seconds = milliseconds / 1000;
    qint64 minutes = seconds / 60;
    qint64 hours = minutes / 60;
    
    seconds %= 60;
    minutes %= 60;
    
    if (hours > 0) {
        return QString("%1h %2m %3s").arg(hours).arg(minutes).arg(seconds);
    } else if (minutes > 0) {
        return QString("%1m %2s").arg(minutes).arg(seconds);
    } else {
        return QString("%1s").arg(seconds);
    }
}

void StellinaProcessor::saveProcessingReport() {
    QString outputDir = getOutputDirectoryForCurrentStage();
    QString reportPath = QDir(outputDir).absoluteFilePath(
        QString("processing_report_%1.txt").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss")));
    
    QFile reportFile(reportPath);
    if (!reportFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
        return;
    }
    
    QTextStream out(&reportFile);
    
    out << "Enhanced Stellina Processor - Processing Report\n";
    out << "==============================================\n\n";
    out << "Date: " << QDateTime::currentDateTime().toString() << "\n";
    out << "Processing Mode: " << m_processingModeCombo->currentText() << "\n";
    out << "Source Directory: " << m_sourceDirectory << "\n";
    out << "Dark Directory: " << m_darkDirectory << "\n\n";
    
    out << "Processing Statistics:\n";
    out << "- Total images: " << m_imagesToProcess.length() << "\n";
    out << "- Successfully processed: " << m_processedCount << "\n";
    out << "- Dark calibrated: " << m_darkCalibratedCount << "\n";
    out << "- Errors: " << m_errorCount << "\n";
    out << "- Skipped: " << m_skippedCount << "\n\n";
    if (!m_finalStackedImage.isEmpty()) {
        out << "Final stacked image: " << m_finalStackedImage << "\n\n";
    }
    
    qint64 totalTime = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
    out << "Total processing time: " << formatProcessingTime(totalTime) << "\n";
    
    reportFile.close();
    logMessage(QString("Processing report saved: %1").arg(QFileInfo(reportPath).fileName()), "blue");
}
// Add these test functions to StellinaProcessor_Core.cpp
// Call testLibnovaConversion() from your constructor or a menu action for verification

void StellinaProcessor::testLibnovaConversion() {
    logMessage("=== Testing libnova coordinate conversion ===", "blue");
    
    // Test case 1: First image from your log
    // Blind solve result: RA=10.9127702516°, Dec=41.2122118775°
    // Your current result: RA=132.720487°, Dec=22.962432°
    testSingleConversion(
        "Test 1 - First image",
        42.0410,           // Alt from log
        286.8526,          // Az from log  
        "2024-01-09T22:13:29",  // DATE-OBS from log
        10.9127702516,     // Expected RA from blind solve
        41.2122118775,     // Expected Dec from blind solve
        132.720487,        // Current RA result
        22.962432          // Current Dec result
    );
    
    // Test case 2: Known good coordinate conversion
    // Use a well-known object at a specific time for validation
    // M31 (Andromeda Galaxy) coordinates: RA=10.684°, Dec=41.269°
    testSingleConversion(
        "Test 2 - M31 reference",
        45.0,              // Approximate alt for M31 from London
        290.0,             // Approximate az for M31 from London
        "2024-01-09T22:00:00",  // Round time
        10.684,            // M31 RA (known)
        41.269,            // M31 Dec (known)
        0.0,               // Will be calculated
        0.0                // Will be calculated
    );
    
    // Test case 3: Different time to check time dependency
    testSingleConversion(
        "Test 3 - Time dependency check",
        42.0410,           // Same Alt as test 1
        286.8526,          // Same Az as test 1
        "2024-01-09T23:13:29",  // 1 hour later
        0.0,               // Expected will be different due to time
        0.0,               // Expected will be different due to time
        0.0,               // Will be calculated
        0.0                // Will be calculated
    );
    
    // Test case 4: Different observer location
    testSingleConversion(
        "Test 4 - Different location (Paris)",
        42.0410,           // Same Alt
        286.8526,          // Same Az
        "2024-01-09T22:13:29",  // Same time
        0.0,               // Expected will be different due to location
        0.0,               // Expected will be different due to location
        0.0,               // Will be calculated
        0.0,               // Will be calculated
        48.8566,           // Paris latitude
        2.3522             // Paris longitude
    );
    
    // Test case 5: Zenith pointing (RA should equal LST, Dec should equal latitude)
    // For 2024-01-09T22:13:29 at London, LST ≈ 17.5 hours ≈ 262.5°
    testSingleConversion(
        "Test 5 - Zenith pointing",
        90.0,              // Pointing straight up
        0.0,               // Azimuth doesn't matter at zenith
        "2024-01-09T22:13:29",
        262.5,             // Expected RA ≈ LST (will verify in test)
        51.5074,           // Dec should equal observer latitude
        0.0,               // Will be calculated
        0.0                // Will be calculated
    );
    
    logMessage("=== libnova coordinate conversion tests complete ===", "blue");
}

void StellinaProcessor::testSingleConversion(const QString &testName,
                                           double alt, double az, 
                                           const QString &dateObs,
                                           double expectedRA, double expectedDec,
                                           double currentRA, double currentDec,
                                           double testLat, double testLon) {
    
    logMessage(QString("--- %1 ---").arg(testName), "green");
    logMessage(QString("Input: Alt=%1°, Az=%2°, Time=%3")
                  .arg(alt, 0, 'f', 4)
                  .arg(az, 0, 'f', 4)
                  .arg(dateObs), "gray");
    
    // Save current observer location
    QString savedLocation = m_observerLocation;
    
    // Set test location if provided
    if (testLat != 0.0 || testLon != 0.0) {
        m_observerLocation = QString("%1,%2").arg(testLat, 0, 'f', 4).arg(testLon, 0, 'f', 4);
        logMessage(QString("Using test location: %1°N, %2°E")
                      .arg(testLat, 0, 'f', 4).arg(testLon, 0, 'f', 4), "gray");
    }
    
    // Perform conversion
    double calculatedRA, calculatedDec;
    bool success = convertAltAzToRaDec(alt, az, dateObs, calculatedRA, calculatedDec);
    
    if (success) {
        logMessage(QString("Calculated: RA=%1°, Dec=%2°")
                      .arg(calculatedRA, 0, 'f', 6)
                      .arg(calculatedDec, 0, 'f', 6), "blue");
        
        // Show in HMS/DMS format too
        double ra_hours = calculatedRA / 15.0;
        int h = static_cast<int>(ra_hours);
        int m = static_cast<int>((ra_hours - h) * 60);
        double s = ((ra_hours - h) * 60 - m) * 60;
        
        int d = static_cast<int>(calculatedDec);
        int am = static_cast<int>(qAbs(calculatedDec - d) * 60);
        double as = (qAbs(calculatedDec - d) * 60 - am) * 60;
        
        logMessage(QString("          RA=%1h%2m%3s, Dec=%4°%5'%6\"")
                      .arg(h).arg(m, 2, 10, QChar('0')).arg(s, 0, 'f', 1)
                      .arg(d).arg(am, 2, 10, QChar('0')).arg(as, 0, 'f', 1), "blue");
        
        // Compare with expected results if provided
        if (expectedRA != 0.0 || expectedDec != 0.0) {
            double raError = qAbs(calculatedRA - expectedRA);
            double decError = qAbs(calculatedDec - expectedDec);
            
            logMessage(QString("Expected: RA=%1°, Dec=%2°")
                          .arg(expectedRA, 0, 'f', 6)
                          .arg(expectedDec, 0, 'f', 6), "orange");
            logMessage(QString("Error:    RA=%1° (%2 arcmin), Dec=%3° (%4 arcmin)")
                          .arg(raError, 0, 'f', 6).arg(raError * 60, 0, 'f', 1)
                          .arg(decError, 0, 'f', 6).arg(decError * 60, 0, 'f', 1), 
                      (raError < 1.0 && decError < 1.0) ? "green" : "red");
            
            // Check if we're close to expected (within 1 degree)
            if (raError < 1.0 && decError < 1.0) {
                logMessage("✓ PASS: Within 1° tolerance", "green");
            } else {
                logMessage("✗ FAIL: Outside 1° tolerance", "red");
                
                // Try to diagnose the issue
                if (raError > 10.0) {
                    logMessage("Large RA error suggests coordinate system issue", "red");
                }
                if (qAbs(calculatedRA - currentRA) < 1.0) {
                    logMessage("Current result matches - may be systematic error", "orange");
                }
            }
        }
        
        // Compare with current result if provided
        if (currentRA != 0.0 || currentDec != 0.0) {
            logMessage(QString("Previous: RA=%1°, Dec=%2°")
                          .arg(currentRA, 0, 'f', 6)
                          .arg(currentDec, 0, 'f', 6), "gray");
        }
        
    } else {
        logMessage("✗ CONVERSION FAILED", "red");
    }
    
    // Restore original observer location
    m_observerLocation = savedLocation;
    
    logMessage("", "gray");  // Blank line for separation
}

bool StellinaProcessor::performCFABinning(const std::vector<float> &inputPixels, std::vector<float> &binnedPixels,
                                         long width, long height, long &binnedWidth, long &binnedHeight) {
    // Ensure dimensions are even for 2x2 binning
    if (width % 2 != 0 || height % 2 != 0) {
        logMessage("Image dimensions must be even for 2x2 CFA binning", "red");
        return false;
    }
    
    binnedWidth = width / 2;
    binnedHeight = height / 2;
    binnedPixels.resize(binnedWidth * binnedHeight);
    
    if (m_debugMode) {
        logMessage(QString("CFA binning: %1x%2 → %3x%4").arg(width).arg(height).arg(binnedWidth).arg(binnedHeight), "gray");
    }
    
    // Perform 2x2 CFA binning - add all 4 pixels in each 2x2 block and scale down
    // For Stellina RGGB pattern: R G
    //                            G B
    for (long y = 0; y < binnedHeight; ++y) {
        for (long x = 0; x < binnedWidth; ++x) {
            long srcY = y * 2;
            long srcX = x * 2;
            
            // Get the 4 pixels from the 2x2 block
            long idx00 = srcY * width + srcX;         // R (top-left)
            long idx01 = srcY * width + (srcX + 1);   // G (top-right)  
            long idx10 = (srcY + 1) * width + srcX;   // G (bottom-left)
            long idx11 = (srcY + 1) * width + (srcX + 1); // B (bottom-right)
            
            // Add all 4 pixels and scale down by 4 to prevent overflow
            // This preserves the proper intensity levels for 16-bit FITS
            float binnedValue = (inputPixels[idx00] + inputPixels[idx01] + 
                               inputPixels[idx10] + inputPixels[idx11]) / 4.0f;
            
            long binnedIdx = y * binnedWidth + x;
            binnedPixels[binnedIdx] = binnedValue;
        }
    }
    
    return true;
}

bool StellinaProcessor::createBinnedImageForPlatesolving(const QString &inputPath, const QString &binnedPath) {
    logMessage(QString("Creating 2x2 binned image for plate solving: %1").arg(QFileInfo(binnedPath).fileName()), "blue");
    
    fitsfile *inputFits = nullptr;
    fitsfile *outputFits = nullptr;
    int status = 0;
    
    // Open input FITS file
    QByteArray inputPathBytes = inputPath.toLocal8Bit();
    if (fits_open_file(&inputFits, inputPathBytes.data(), READONLY, &status)) {
        logMessage(QString("Failed to open input FITS: %1 (error: %2)").arg(inputPath).arg(status), "red");
        return false;
    }
    
    // Get image dimensions
    long naxes[2];
    if (fits_get_img_size(inputFits, 2, naxes, &status)) {
        logMessage(QString("Failed to get image dimensions (error: %1)").arg(status), "red");
        fits_close_file(inputFits, &status);
        return false;
    }
    
    long width = naxes[0];
    long height = naxes[1];
    long totalPixels = width * height;
    
    // Read input pixel data
    std::vector<float> inputPixels(totalPixels);
    long firstPixel = 1;
    if (fits_read_img(inputFits, TFLOAT, firstPixel, totalPixels, nullptr, 
                     inputPixels.data(), nullptr, &status)) {
        logMessage(QString("Failed to read input pixel data (error: %1)").arg(status), "red");
        fits_close_file(inputFits, &status);
        return false;
    }
    
    // Perform CFA-aware 2x2 binning
    std::vector<float> binnedPixels;
    long binnedWidth, binnedHeight;
    
    if (!performCFABinning(inputPixels, binnedPixels, width, height, binnedWidth, binnedHeight)) {
        fits_close_file(inputFits, &status);
        return false;
    }
    
    // Create output FITS file
    QByteArray outputPathBytes = QString("!%1").arg(binnedPath).toLocal8Bit(); // ! forces overwrite
    if (fits_create_file(&outputFits, outputPathBytes.data(), &status)) {
        logMessage(QString("Failed to create binned FITS: %1 (error: %2)").arg(binnedPath).arg(status), "red");
        fits_close_file(inputFits, &status);
        return false;
    }
    
    // Copy header from input, but update dimensions
    if (fits_copy_header(inputFits, outputFits, &status)) {
        logMessage(QString("Failed to copy FITS header (error: %1)").arg(status), "red");
        fits_close_file(inputFits, &status);
        fits_close_file(outputFits, &status);
        return false;
    }
    
    // Update NAXIS1 and NAXIS2 for binned dimensions
    if (fits_update_key(outputFits, TLONG, "NAXIS1", &binnedWidth, "Binned image width", &status) ||
        fits_update_key(outputFits, TLONG, "NAXIS2", &binnedHeight, "Binned image height", &status)) {
        logMessage("Warning: Could not update image dimensions in header", "orange");
        status = 0; // Continue anyway
    }
    
    // Write binned pixel data
    long binnedTotalPixels = binnedWidth * binnedHeight;
    if (fits_write_img(outputFits, TFLOAT, firstPixel, binnedTotalPixels, 
                      binnedPixels.data(), &status)) {
        logMessage(QString("Failed to write binned pixel data (error: %1)").arg(status), "red");
        fits_close_file(inputFits, &status);
        fits_close_file(outputFits, &status);
        QFile::remove(binnedPath);
        return false;
    }
    
    // Add processing history
    QString historyComment = QString("CFA-aware 2x2 binning with /4 scaling for plate solving (%1x%2 -> %3x%4)")
                                .arg(width).arg(height).arg(binnedWidth).arg(binnedHeight);
    QByteArray historyBytes = historyComment.toLocal8Bit();
    if (fits_write_history(outputFits, historyBytes.data(), &status)) {
        logMessage("Warning: Could not write binning history", "orange");
        status = 0;
    }
    
    // Add custom keywords
    QString purposeStr = "PLATESOLVE";
    QByteArray purposeBytes = purposeStr.toLocal8Bit();
    char* purposePtr = purposeBytes.data();
    if (fits_write_key(outputFits, TSTRING, "PURPOSE", &purposePtr, "Binned for plate solving", &status)) {
        logMessage("Warning: Could not write PURPOSE keyword", "orange");
        status = 0;
    }
    
    int binFactor = 2;
    if (fits_write_key(outputFits, TINT, "BINNING", &binFactor, "CFA binning factor", &status)) {
        logMessage("Warning: Could not write BINNING keyword", "orange");
        status = 0;
    }
    
    fits_close_file(inputFits, &status);
    fits_close_file(outputFits, &status);
    
    if (status != 0) {
        logMessage(QString("FITS error during binned image creation (error: %1)").arg(status), "red");
        QFile::remove(binnedPath);
        return false;
    }
    
    logMessage(QString("Created binned image: %1 (%2x%3 pixels)")
                  .arg(QFileInfo(binnedPath).fileName()).arg(binnedWidth).arg(binnedHeight), "green");
    return true;
}


// Enhanced coordinate extraction that can handle both JSON and FITS sources
bool StellinaProcessor::extractCoordinatesFromAnySource(const QString &fitsPath, 
                                                       const QString &jsonPath, 
                                                       CoordinateSource &coords) {
    // First priority: Try to extract from JSON if it exists (original method)
    if (!jsonPath.isEmpty() && QFile::exists(jsonPath)) {
        QJsonObject metadata = loadStellinaJson(jsonPath);
        if (!metadata.isEmpty()) {
            double alt, az;
            if (extractCoordinates(metadata, alt, az)) {
                // Convert Alt/Az to RA/Dec using existing method
                QString dateObs = extractDateObs(fitsPath);
                if (!dateObs.isEmpty()) {
                    double ra, dec;
                    if (convertAltAzToRaDec(alt, az, dateObs, ra, dec)) {
                        coords.type = CoordinateSource::FROM_JSON_ALTAZ;
                        coords.ra = ra;
                        coords.dec = dec;
                        coords.alt = alt;
                        coords.az = az;
                        coords.source_info = QString("JSON Alt/Az: %1°,%2° -> RA/Dec: %3°,%4°")
                                           .arg(alt, 0, 'f', 4).arg(az, 0, 'f', 4)
                                           .arg(ra, 0, 'f', 6).arg(dec, 0, 'f', 6);
                        
                        if (m_debugMode) {
                            logMessage(QString("Using JSON coordinates: %1").arg(coords.source_info), "green");
                        }
                        return true;
                    }
                }
            }
        }
    }
    
    // Second priority: Try to extract WCS coordinates directly from FITS
    if (extractWCSCoordinatesFromFITS(fitsPath, coords)) {
        if (m_debugMode) {
            logMessage(QString("Using FITS WCS coordinates: %1").arg(coords.source_info), "green");
        }
        return true;
    }
    
    // No valid coordinates found from either source
    coords.type = CoordinateSource::INVALID;
    return false;
}

// New function to extract WCS coordinates directly from FITS headers
bool StellinaProcessor::extractWCSCoordinatesFromFITS(const QString &fitsPath, CoordinateSource &coords) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Try to read CRVAL1 and CRVAL2 (RA/Dec reference point)
    double crval1, crval2;
    bool hasWCS = true;
    
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &crval1, nullptr, &status) != 0) {
        hasWCS = false;
    }
    
    status = 0; // Reset status for next read
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &crval2, nullptr, &status) != 0) {
        hasWCS = false;
    }
    
    if (hasWCS) {
        // Successfully read WCS coordinates
        coords.type = CoordinateSource::FROM_FITS_WCS;
        coords.ra = crval1;
        coords.dec = crval2;
        coords.alt = 0;  // Not available from WCS
        coords.az = 0;   // Not available from WCS
        coords.source_info = QString("FITS WCS: CRVAL1=%1°, CRVAL2=%2°")
                           .arg(crval1, 0, 'f', 6).arg(crval2, 0, 'f', 6);
        
        fits_close_file(fptr, &status);
        return true;
    }
    
    fits_close_file(fptr, &status);
    return false;
}

// Modified findStellinaImages function to handle both coordinate sources
bool StellinaProcessor::findStellinaImages() {
    if (m_sourceDirectory.isEmpty()) {
        logMessage("No source directory selected", "red");
        return false;
    }
    
    QDir sourceDir(m_sourceDirectory);
    if (!sourceDir.exists()) {
        logMessage("Source directory does not exist", "red");
        return false;
    }
    
    m_stellinaImageData.clear();
    m_imagesToProcess.clear();
    
    // Find all FITS files
    QStringList fitsFiles = sourceDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    fitsFiles.sort();
    
    logMessage(QString("Found %1 FITS files").arg(fitsFiles.size()), "blue");
    
    int validPairs = 0;
    int jsonMissing = 0;
    int wcsAvailable = 0;
    int totalInvalid = 0;
    
    for (const QString &fitsFile : fitsFiles) {
        QString fitsPath = sourceDir.absoluteFilePath(fitsFile);
        QString baseName = QFileInfo(fitsFile).baseName();
        
        StellinaImageData imageData;
        imageData.originalFitsPath = fitsPath;
        imageData.currentFitsPath = fitsPath;
        
        // Try to find corresponding JSON file (original behavior)
        QStringList jsonCandidates = {
            baseName + ".json",
            baseName + ".JSON",
            baseName + "-stacking.json",
            baseName + "-stacking.JSON"
        };
        
        QString jsonPath;
        for (const QString &candidate : jsonCandidates) {
            QString candidatePath = sourceDir.absoluteFilePath(candidate);
            if (QFile::exists(candidatePath)) {
                jsonPath = candidatePath;
                imageData.originalJsonPath = candidatePath;
                break;
            }
        }
        
        // Extract coordinates from any available source
        CoordinateSource coords;
        if (extractCoordinatesFromAnySource(fitsPath, jsonPath, coords)) {
            // Populate imageData based on coordinate source
            if (coords.type == CoordinateSource::FROM_JSON_ALTAZ) {
                // Original JSON-based workflow
                imageData.metadata = loadStellinaJson(jsonPath);
                imageData.altitude = coords.alt;
                imageData.azimuth = coords.az;
                imageData.hasValidCoordinates = true;
                
                // Quality filtering for JSON-based coordinates
                if (m_qualityFilter && !checkStellinaQuality(imageData.metadata)) {
                    if (m_debugMode) {
                        logMessage(QString("Rejected due to quality filter: %1").arg(fitsFile), "orange");
                    }
                    continue;
                }
                
                if (jsonPath.isEmpty()) jsonMissing++; // This shouldn't happen in this branch
                
            } else if (coords.type == CoordinateSource::FROM_FITS_WCS) {
                // New WCS-based workflow
                imageData.altitude = coords.alt;   // Will be 0
                imageData.azimuth = coords.az;     // Will be 0
                imageData.hasValidCoordinates = true;
                
                // Create minimal metadata object for compatibility
                QJsonObject minimalMetadata;
                minimalMetadata["coordinate_source"] = "fits_wcs";
                minimalMetadata["crval1"] = coords.ra;
                minimalMetadata["crval2"] = coords.dec;
                imageData.metadata = minimalMetadata;
                
                wcsAvailable++;
            }
            
            // Extract additional FITS metadata (common to both workflows)
            imageData.exposureSeconds = extractExposureTime(fitsPath);
            imageData.temperatureKelvin = extractTemperature(fitsPath);
            imageData.binning = extractBinning(fitsPath);
            imageData.dateObs = extractDateObs(fitsPath);
            imageData.bayerPattern = detectBayerPattern(fitsPath);
            
            // Store calculated coordinates for both sources
            imageData.calculatedRA = coords.ra;
            imageData.calculatedDec = coords.dec;
            imageData.hasCalculatedCoords = true;
            
            m_stellinaImageData.append(imageData);
            m_imagesToProcess.append(fitsPath);
            validPairs++;
            
            if (m_debugMode) {
                logMessage(QString("Valid image: %1 - %2").arg(fitsFile).arg(coords.source_info), "gray");
            }
        } else {
            totalInvalid++;
            if (m_debugMode) {
                logMessage(QString("No valid coordinates found for: %1").arg(fitsFile), "orange");
            }
        }
    }
    
    // Summary report
    logMessage("", "gray");
    logMessage("=== COORDINATE SOURCE SUMMARY ===", "blue");
    logMessage(QString("Total FITS files: %1").arg(fitsFiles.size()), "gray");
    logMessage(QString("Valid images found: %1").arg(validPairs), "green");
    if (validPairs - wcsAvailable > 0) {
        logMessage(QString("  - From JSON Alt/Az: %1").arg(validPairs - wcsAvailable), "gray");
    }
    if (wcsAvailable > 0) {
        logMessage(QString("  - From FITS WCS: %1").arg(wcsAvailable), "gray");
    }
    if (totalInvalid > 0) {
        logMessage(QString("Invalid/skipped: %1").arg(totalInvalid), "orange");
    }
    logMessage("", "gray");
    
    if (validPairs == 0) {
        logMessage("No valid image pairs found. Check that files have either:", "red");
        logMessage("  1. Corresponding JSON files with Alt/Az coordinates, OR", "red");
        logMessage("  2. FITS headers with CRVAL1/CRVAL2 WCS coordinates", "red");
        return false;
    }
    
    return true;
}

// Modified writeStellinaMetadataWithCoordinates to handle WCS sources
bool StellinaProcessor::writeStellinaMetadataWithCoordinates(const QString &fitsPath, const StellinaImageData &imageData) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READWRITE, &status)) {
        logMessage(QString("Failed to open FITS file for metadata writing: %1 (status: %2)").arg(fitsPath).arg(status), "red");
        return false;
    }
    
    // Check if coordinates came from WCS (in which case they're already in the file)
    bool coordsFromWCS = imageData.metadata.contains("coordinate_source") && 
                        imageData.metadata["coordinate_source"].toString() == "fits_wcs";
    
    if (!coordsFromWCS) {
        // Original workflow: write Alt/Az and calculated RA/Dec
        if (imageData.hasValidCoordinates && imageData.hasCalculatedCoords) {
            // Write Stellina Alt/Az coordinates
            double alt = imageData.altitude;
            double az = imageData.azimuth;
            
            if (fits_write_key(fptr, TDOUBLE, "STELLALT", &alt, "Stellina mount altitude (degrees)", &status) ||
                fits_write_key(fptr, TDOUBLE, "STELLAZ", &az, "Stellina mount azimuth (degrees)", &status)) {
                logMessage("Warning: Could not write Stellina Alt/Az coordinates", "orange");
                status = 0;
            }
            
            // Write calculated RA/Dec coordinates
            double ra = imageData.calculatedRA;
            double dec = imageData.calculatedDec;
            
            if (fits_write_key(fptr, TDOUBLE, "STELLRA", &ra, "Calculated RA from mount position (degrees)", &status) ||
                fits_write_key(fptr, TDOUBLE, "STELLDEC", &dec, "Calculated Dec from mount position (degrees)", &status)) {
                logMessage("Warning: Could not write calculated RA/Dec coordinates", "orange");
                status = 0;
            }
            
            // Write coordinate conversion metadata
            QString conversion_method = "Alt/Az to RA/Dec from Stellina mount position";
            QByteArray methodBytes = conversion_method.toLocal8Bit();
            char* methodPtr = methodBytes.data();
            if (fits_write_key(fptr, TSTRING, "COORDMET", &methodPtr, "Mount coordinate conversion method", &status)) {
                logMessage("Warning: Could not write conversion method", "orange");
                status = 0;
            }
        }
        
        // Write original file paths if they exist
        if (!imageData.originalJsonPath.isEmpty()) {
            QString origJsonRelative = QFileInfo(imageData.originalJsonPath).fileName();
            QByteArray origJsonBytes = origJsonRelative.toLocal8Bit();
            char* origJsonPtr = origJsonBytes.data();
            
            if (fits_write_key(fptr, TSTRING, "STELLJSON", &origJsonPtr, "Original Stellina JSON file", &status)) {
                logMessage("Warning: Could not write original JSON file reference", "orange");
                status = 0;
            }
        }
    } else {
        // WCS workflow: coordinates are already in CRVAL1/CRVAL2, just mark the source
        QString wcs_source = "WCS coordinates from FITS headers (CRVAL1/CRVAL2)";
        QByteArray sourceBytes = wcs_source.toLocal8Bit();
        char* sourcePtr = sourceBytes.data();
        if (fits_write_key(fptr, TSTRING, "COORDMET", &sourcePtr, "Coordinate source method", &status)) {
            logMessage("Warning: Could not write coordinate source method", "orange");
            status = 0;
        }
    }
    
    // Write common metadata (exposure, temperature, etc.) regardless of coordinate source
    if (imageData.exposureSeconds > 0) {
        int exposure = imageData.exposureSeconds;
        if (fits_write_key(fptr, TINT, "STELLEXP", &exposure, "Exposure time (seconds)", &status)) {
            logMessage("Warning: Could not write exposure time", "orange");
            status = 0;
        }
    }
    
    if (imageData.temperatureKelvin > 0) {
        int temperature = imageData.temperatureKelvin;
        if (fits_write_key(fptr, TINT, "STELLTMP", &temperature, "Sensor temperature (Kelvin)", &status)) {
            logMessage("Warning: Could not write temperature", "orange");
            status = 0;
        }
    }
    
    if (!imageData.binning.isEmpty()) {
        QByteArray binningBytes = imageData.binning.toLocal8Bit();
        char* binningPtr = binningBytes.data();
        if (fits_write_key(fptr, TSTRING, "STELLBIN", &binningPtr, "Binning mode", &status)) {
            logMessage("Warning: Could not write binning mode", "orange");
            status = 0;
        }
    }
    
    fits_close_file(fptr, &status);
    return true;
}
