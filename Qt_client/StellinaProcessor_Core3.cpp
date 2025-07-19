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
// Modified findStellinaImages function to support reversed images
bool StellinaProcessor::findStellinaImages() {
    m_imagesToProcess.clear();
    m_stellinaImageData.clear();
    
    QDir sourceDir(m_sourceDirectory);
    if (!sourceDir.exists()) {
        logMessage(QString("Source directory does not exist: '%1'").arg(m_sourceDirectory), "red");
        return false;
    }
    
    logMessage(QString("Scanning directory: %1").arg(sourceDir.absolutePath()), "blue");
    
    QStringList fitsFiles = sourceDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    logMessage(QString("Found %1 FITS files").arg(fitsFiles.count()), "blue");
    
    int validPairs = 0;
    int jsonMissing = 0;
    int qualityRejected = 0;
    int reversedImages = 0;
    
    for (const QString &fitsFile : fitsFiles) {
        StellinaImageData imageData;
        imageData.originalFitsPath = sourceDir.absoluteFilePath(fitsFile);
        imageData.currentFitsPath = imageData.originalFitsPath;
        
        // Check if this is a reversed stellina image
        imageData.isReversedImage = isReversedStellinaImage(fitsFile);
        if (imageData.isReversedImage) {
            reversedImages++;
            logMessage(QString("Found reversed Stellina image: %1").arg(fitsFile), "blue");
        }
        
        // Get base name (without 'r' suffix for reversed images)
        imageData.baseName = getBaseName(fitsFile);
        
        // Detect bayer pattern
        imageData.bayerPattern = detectBayerPattern(imageData.originalFitsPath);
        
        QString baseName = QFileInfo(fitsFile).baseName();
        
        // Try to find corresponding JSON file
        // For reversed images, look for JSON files both with and without 'r' suffix
        QStringList jsonCandidates;
        
        if (imageData.isReversedImage) {
            // For img-0001r.fits, try both img-0001r.json and img-0001.json
            jsonCandidates << baseName + ".json"
                          << baseName + ".JSON"
                          << imageData.baseName + ".json"
                          << imageData.baseName + ".JSON"
                          << baseName + "-stacking.json"
                          << baseName + "-stacking.JSON"
                          << imageData.baseName + "-stacking.json"
                          << imageData.baseName + "-stacking.JSON";
        } else {
            // Regular candidates for normal images
            jsonCandidates << baseName + ".json"
                          << baseName + ".JSON"
                          << baseName + "-stacking.json"
                          << baseName + "-stacking.JSON";
        }
        
        // Add complete base name candidates
        jsonCandidates << QFileInfo(fitsFile).completeBaseName() + ".json"
                      << QFileInfo(fitsFile).completeBaseName() + ".JSON"
                      << QFileInfo(fitsFile).completeBaseName() + "-stacking.json"
                      << QFileInfo(fitsFile).completeBaseName() + "-stacking.JSON";
        
        bool jsonFound = false;
        for (const QString &candidate : jsonCandidates) {
            QString candidatePath = sourceDir.absoluteFilePath(candidate);
            qDebug() << candidatePath;
            if (QFile::exists(candidatePath)) {
                imageData.originalJsonPath = candidatePath;
                jsonFound = true;
                break;
            }
        }
        
        if (!jsonFound) {
            if (m_debugMode) {
                logMessage(QString("No JSON metadata found for: %1").arg(fitsFile), "orange");
            }
            jsonMissing++;
            continue;
        }
        
        // Load and parse JSON metadata
        imageData.metadata = loadStellinaJson(imageData.originalJsonPath);
        if (imageData.metadata.isEmpty()) {
            logMessage(QString("Failed to parse JSON for %1").arg(fitsFile), "red");
            continue;
        }

        // Extract coordinates from JSON
        if (!extractCoordinates(imageData.metadata, imageData.altitude, imageData.azimuth)) {
            logMessage(QString("Failed to extract metadata from: %1").arg(imageData.originalJsonPath), "red");
            continue;
        }

        imageData.hasValidCoordinates = true;
        
        // Extract FITS metadata
        imageData.exposureSeconds = extractExposureTime(imageData.originalFitsPath);
        imageData.temperatureKelvin = extractTemperature(imageData.originalFitsPath);
        imageData.binning = extractBinning(imageData.originalFitsPath);
        imageData.dateObs = extractDateObs(imageData.originalFitsPath);
        
        // Quality filtering
        if (m_qualityFilter && !checkStellinaQuality(imageData.metadata)) {
            qualityRejected++;
            if (m_debugMode) {
                logMessage(QString("Rejected due to invalid coordinates: %1").arg(fitsFile), "orange");
            }
            qualityRejected++;
            continue;
        }
        
        m_stellinaImageData.append(imageData);
        m_imagesToProcess.append(imageData.originalFitsPath);
        validPairs++;
        
        if (m_debugMode) {
            QString typeInfo = imageData.isReversedImage ? " (REVERSED)" : "";
            logMessage(QString("Added: %1 -> %2%3 (Bayer: %4)")
                          .arg(QFileInfo(fitsFile).fileName())
                          .arg(QFileInfo(imageData.originalJsonPath).fileName())
                          .arg(typeInfo)
                          .arg(imageData.bayerPattern), "gray");
        }
    }
    
    logMessage(QString("Results: %1 valid pairs, %2 missing JSON, %3 quality rejected, %4 reversed images")
                  .arg(validPairs).arg(jsonMissing).arg(qualityRejected).arg(reversedImages), "blue");
    
    return !m_stellinaImageData.isEmpty();
}
