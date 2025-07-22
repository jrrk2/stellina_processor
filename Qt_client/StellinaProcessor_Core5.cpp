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

// Direct coordinate data extraction without plate solving
// This will help us analyze the time drift issue efficiently

void StellinaProcessor::dumpCoordinateData() {
    logMessage("=== DUMPING COORDINATE DATA FROM JSON/FITS PAIRS ===", "blue");
    
    if (m_sourceDirectory.isEmpty()) {
        logMessage("Please select source directory first", "red");
        return;
    }
    
    QDir sourceDir(m_sourceDirectory);
    if (!sourceDir.exists()) {
        logMessage(QString("Source directory does not exist: %1").arg(m_sourceDirectory), "red");
        return;
    }
    
    // Find all FITS files
    QStringList fitsFiles = sourceDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    fitsFiles.sort(); // Ensure chronological order
    
    logMessage(QString("Found %1 FITS files").arg(fitsFiles.size()), "blue");
    logMessage("", "gray");
    
    // Header for the data dump
    logMessage("Data Format: Image | Time | Alt | Az | Calculated_RA | Calculated_Dec | DATE-OBS", "green");
    logMessage("====================================================================================================", "gray");
    
    int validCount = 0;
    QDateTime startTime;
    
    for (const QString &fitsFile : fitsFiles) {
        QString fitsPath = sourceDir.absoluteFilePath(fitsFile);
        QString baseName = QFileInfo(fitsFile).baseName();
        
        // Try to find corresponding JSON file
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
                break;
            }
        }
        
        if (jsonPath.isEmpty()) {
            continue; // Skip if no JSON found
        }
        
        // Load JSON metadata
        QJsonObject metadata = loadStellinaJson(jsonPath);
        if (metadata.isEmpty()) {
            continue;
        }
        
        // Extract coordinates from JSON
        double alt, az;
        if (!extractCoordinates(metadata, alt, az)) {
            continue;
        }
        
        // Extract DATE-OBS from FITS
        QString dateObs = extractDateObs(fitsPath);
        if (dateObs.isEmpty()) {
            continue;
        }
        
        // Calculate RA/Dec using current conversion
        double calculatedRA, calculatedDec;
        if (!convertAltAzToRaDec(alt, az, dateObs, calculatedRA, calculatedDec)) {
            continue;
        }
        
        // Parse time for elapsed calculation
        QDateTime obsTime = QDateTime::fromString(dateObs, "yyyy-MM-ddThh:mm:ss");
        obsTime.setTimeSpec(Qt::UTC);
        
        if (validCount == 0) {
            startTime = obsTime;
        }
        
        double minutesElapsed = startTime.msecsTo(obsTime) / 60000.0;
        
        // Output the data
        logMessage(QString("%1 | %2 | %3 | %4 | %5 | %6 | %7")
                      .arg(QFileInfo(fitsFile).baseName(), -15)
                      .arg(QString::number(minutesElapsed, 'f', 2), 6)
                      .arg(QString::number(alt, 'f', 4), 8)
                      .arg(QString::number(az, 'f', 4), 9)
                      .arg(QString::number(calculatedRA, 'f', 6), 12)
                      .arg(QString::number(calculatedDec, 'f', 6), 13)
                      .arg(dateObs), "gray");
        
        validCount++;
        
        // Stop after reasonable number for initial analysis
        if (validCount >= 50) {
            logMessage(QString("... (showing first 50 of %1 total files)").arg(fitsFiles.size()), "blue");
            break;
        }
    }
    
    logMessage("====================================================================================================", "gray");
    logMessage(QString("Processed %1 valid image pairs").arg(validCount), "green");
    logMessage("", "gray");
    logMessage("ANALYSIS INSTRUCTIONS:", "blue");
    logMessage("1. Copy this data to a spreadsheet or analysis tool", "gray");
    logMessage("2. Plot Calculated_RA vs Time to see drift pattern", "gray");
    logMessage("3. Look for systematic increase/decrease over time", "gray");
    logMessage("4. Compare with expected RA values if available", "gray");
}

// Alternative version that saves to CSV file for easier analysis
void StellinaProcessor::dumpCoordinateDataToCSV() {
    logMessage("=== DUMPING COORDINATE DATA TO CSV FILE ===", "blue");
    
    if (m_sourceDirectory.isEmpty()) {
        logMessage("Please select source directory first", "red");
        return;
    }
    
    QDir sourceDir(m_sourceDirectory);
    if (!sourceDir.exists()) {
        logMessage(QString("Source directory does not exist: %1").arg(m_sourceDirectory), "red");
        return;
    }
    
    // Create output CSV file
    QString csvPath = QDir(m_sourceDirectory).absoluteFilePath("stellina_coordinates.csv");
    QFile csvFile(csvPath);
    
    if (!csvFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
        logMessage(QString("Failed to create CSV file: %1").arg(csvPath), "red");
        return;
    }
    
    QTextStream csv(&csvFile);
    
    // Write CSV header
    csv << "image_name,minutes_elapsed,altitude,azimuth,calculated_ra,calculated_dec,date_obs,julian_day,lst_hours\n";
    
    // Find all FITS files
    QStringList fitsFiles = sourceDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    fitsFiles.sort();
    
    int validCount = 0;
    QDateTime startTime;
    
    for (const QString &fitsFile : fitsFiles) {
        QString fitsPath = sourceDir.absoluteFilePath(fitsFile);
        QString baseName = QFileInfo(fitsFile).baseName();
        
        // Find JSON file
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
                break;
            }
        }
        
        if (jsonPath.isEmpty()) continue;
        
        // Load JSON and extract coordinates
        QJsonObject metadata = loadStellinaJson(jsonPath);
        if (metadata.isEmpty()) continue;
        
        double alt, az;
        if (!extractCoordinates(metadata, alt, az)) continue;
        
        // Extract DATE-OBS
        QString dateObs = extractDateObs(fitsPath);
        if (dateObs.isEmpty()) continue;
        
        // Calculate coordinates
        double calculatedRA, calculatedDec;
        if (!convertAltAzToRaDec(alt, az, dateObs, calculatedRA, calculatedDec)) continue;
        
        // Parse time and calculate elapsed minutes
        QDateTime obsTime = QDateTime::fromString(dateObs, "yyyy-MM-ddThh:mm:ss");
        obsTime.setTimeSpec(Qt::UTC);
        
        if (validCount == 0) {
            startTime = obsTime;
        }
        
        double minutesElapsed = startTime.msecsTo(obsTime) / 60000.0;
        
        // Calculate Julian Day and LST for debugging
        double jd = CoordinateUtils::computeJulianDay(obsTime.date().year(), obsTime.date().month(), obsTime.date().day(),
                               obsTime.time().hour(), obsTime.time().minute(), obsTime.time().second());
        double lst = calculateLST_HighPrecision(jd, -0.1278);
        
        // Write to CSV
        csv << QString("%1,%2,%3,%4,%5,%6,%7,%8,%9\n")
                   .arg(baseName)
                   .arg(minutesElapsed, 0, 'f', 3)
                   .arg(alt, 0, 'f', 6)
                   .arg(az, 0, 'f', 6)
                   .arg(calculatedRA, 0, 'f', 8)
                   .arg(calculatedDec, 0, 'f', 8)
                   .arg(dateObs)
                   .arg(jd, 0, 'f', 8)
                   .arg(lst, 0, 'f', 8);
        
        validCount++;
    }
    
    csvFile.close();
    
    logMessage(QString("Exported %1 coordinate records to: %2").arg(validCount).arg(csvPath), "green");
    logMessage("", "gray");
    logMessage("NEXT STEPS:", "blue");
    logMessage("1. Open the CSV file in Excel/Google Sheets", "gray");
    logMessage("2. Create a plot of calculated_ra vs minutes_elapsed", "gray");
    logMessage("3. Look for linear drift pattern over time", "gray");
    logMessage("4. Compare LST progression to expected sidereal rate", "gray");
}

// Quick analysis function to identify drift in the current session
void StellinaProcessor::analyzeCoordinateDrift() {
    logMessage("=== ANALYZING COORDINATE DRIFT IN CURRENT DATA ===", "blue");
    
    // Test with fixed Alt/Az over time span to isolate time drift
    double fixedAlt = 42.0410;
    double fixedAz = 286.8526;
    
    logMessage(QString("Testing fixed Alt/Az coordinates: %1°, %2°").arg(fixedAlt).arg(fixedAz), "gray");
    logMessage("If RA drifts significantly, the error is in time-dependent conversion", "gray");
    logMessage("", "gray");
    
    // Test over typical Stellina session timespan
    QStringList testTimes = {
        "2024-01-09T22:13:29",  // 0 min
        "2024-01-09T22:23:29",  // +10 min
        "2024-01-09T22:33:29",  // +20 min
        "2024-01-09T22:43:29",  // +30 min
        "2024-01-09T22:53:29"   // +40 min
    };
    
    double firstRA = 0.0;
    double maxDrift = 0.0;
    
    logMessage("Time        | RA      | Dec     | RA Drift | Rate (°/hr)", "green");
    logMessage("=======================================================", "gray");
    
    for (int i = 0; i < testTimes.size(); ++i) {
        double ra, dec;
        if (convertAltAzToRaDec(fixedAlt, fixedAz, testTimes[i], ra, dec)) {
            if (i == 0) {
                firstRA = ra;
                logMessage(QString("%1 | %2 | %3 | %4 | %5")
                              .arg(testTimes[i], -19)
                              .arg(QString::number(ra, 'f', 3), 7)
                              .arg(QString::number(dec, 'f', 3), 7)
                              .arg("0.000", 8)
                              .arg("0.00", 10), "blue");
            } else {
                double raDrift = ra - firstRA;
                double minutes = i * 10.0;
                double ratePerHour = (raDrift / minutes) * 60.0;
                
                maxDrift = qMax(maxDrift, qAbs(raDrift));
                
                logMessage(QString("%1 | %2 | %3 | %4 | %5")
                              .arg(testTimes[i], -19)
                              .arg(QString::number(ra, 'f', 3), 7)
                              .arg(QString::number(dec, 'f', 3), 7)
                              .arg(QString::number(raDrift, 'f', 3), 8)
                              .arg(QString::number(ratePerHour, 'f', 2), 10),
                          (qAbs(raDrift) > 1.0) ? "red" : "blue");
            }
        }
    }
    
    logMessage("=======================================================", "gray");
    logMessage(QString("Maximum RA drift over 40 minutes: %1°").arg(maxDrift, 0, 'f', 3), 
              (maxDrift > 1.0) ? "red" : "green");
    
    if (maxDrift > 1.0) {
        logMessage("", "gray");
        logMessage("DIAGNOSIS: SIGNIFICANT TIME DRIFT DETECTED", "red");
        logMessage("Root cause is in coordinate conversion algorithm", "red");
        logMessage("Likely issues:", "orange");
        logMessage("- Incorrect sidereal time calculation rate", "gray");
        logMessage("- Wrong epoch or time reference frame", "gray");
        logMessage("- Accumulated rounding errors in time calculations", "gray");
    } else {
        logMessage("", "gray");
        logMessage("DIAGNOSIS: Time drift within acceptable limits", "green");
    }
}

// Add these function implementations to StellinaProcessor_Core.cpp

void StellinaProcessor::applyMountTiltCorrection(double &alt, double &az, double inputAlt, double inputAz) {
    if (!m_mountTilt.enableCorrection) {
        alt = inputAlt;
        az = inputAz;
        return;
    }
    
    // Convert tilt parameters to radians
    const double theta_N = m_mountTilt.northTilt * DEG_TO_RAD;
    const double theta_E = m_mountTilt.eastTilt * DEG_TO_RAD;
    
    // Convert input coordinates to radians
    const double input_az_rad = inputAz * DEG_TO_RAD;
    
    // Apply tilt correction model (based on tilt.py algorithm)
    // The correction assumes small angle approximations for mount tilt
    
    // Calculate tilt corrections in Alt/Az space
    const double delta_alt = theta_N * cos(input_az_rad) + theta_E * sin(input_az_rad);
    const double delta_az = (theta_E * cos(input_az_rad) - theta_N * sin(input_az_rad)) / 
                           tan(inputAlt * DEG_TO_RAD);
    
    // Apply corrections (convert delta_az back to degrees)
    alt = inputAlt + (delta_alt * RAD_TO_DEG);
    az = inputAz + (delta_az * RAD_TO_DEG);
    
    // Normalize azimuth to [0, 360) range
    while (az < 0.0) az += 360.0;
    while (az >= 360.0) az -= 360.0;
    
    // Clamp altitude to valid range [-90, 90]
    if (alt > 90.0) alt = 90.0;
    if (alt < -90.0) alt = -90.0;
    
    if (m_debugMode) {
        logMessage(QString("Mount tilt correction applied:"), "blue");
        logMessage(QString("  Input Alt/Az: %1°, %2°")
                      .arg(inputAlt, 0, 'f', 4)
                      .arg(inputAz, 0, 'f', 4), "gray");
        logMessage(QString("  Tilt params: θ_N=%1°, θ_E=%2°")
                      .arg(m_mountTilt.northTilt, 0, 'f', 4)
                      .arg(m_mountTilt.eastTilt, 0, 'f', 4), "gray");
        logMessage(QString("  Corrections: Δalt=%1°, Δaz=%2°")
                      .arg(delta_alt * RAD_TO_DEG, 0, 'f', 4)
                      .arg(delta_az * RAD_TO_DEG, 0, 'f', 4), "gray");
        logMessage(QString("  Corrected Alt/Az: %1°, %2°")
                      .arg(alt, 0, 'f', 4)
                      .arg(az, 0, 'f', 4), "blue");
    }
}

void StellinaProcessor::calibrateMountTilt() {
    logMessage("=== MOUNT TILT CALIBRATION ===", "blue");
    
    // This would typically be run after collecting plate solving results
    // For now, set the values from your tilt.py analysis
    
    // Values from your tilt analysis: θ_N = 1.0832°, θ_E = 2.4314°
    m_mountTilt.northTilt = 1.0832;  // degrees
    m_mountTilt.eastTilt = 2.4314;   // degrees
    m_mountTilt.enableCorrection = true;
    
    logMessage(QString("Mount tilt parameters set:"), "green");
    logMessage(QString("  North tilt θ_N = %1°").arg(m_mountTilt.northTilt, 0, 'f', 4), "blue");
    logMessage(QString("  East tilt θ_E = %1°").arg(m_mountTilt.eastTilt, 0, 'f', 4), "blue");
    logMessage(QString("  Correction enabled: %1").arg(m_mountTilt.enableCorrection ? "Yes" : "No"), "blue");
    
    // Update UI
    if (m_enableTiltCorrectionCheck) {
        m_enableTiltCorrectionCheck->setChecked(m_mountTilt.enableCorrection);
        m_northTiltSpin->setValue(m_mountTilt.northTilt);
        m_eastTiltSpin->setValue(m_mountTilt.eastTilt);
        updateTiltUI();
    }
    
    // Save to settings
    saveMountTiltToSettings();
    
    logMessage("Mount tilt calibration complete. Correction will be applied to all subsequent processing.", "green");
}

void StellinaProcessor::testMountTiltCorrection() {
    logMessage("=== TESTING MOUNT TILT CORRECTION ===", "blue");
    
    // Test data from your log (known problematic cases)
    struct TestCase {
        QString name;
        double inputAlt, inputAz;
        QString time;
        double expectedRA, expectedDec;  // From solve-field
    };
    
    QList<TestCase> testCases = {
        {"img-0001", 42.0410, 286.8526, "2024-01-09T22:13:29", 10.6760, 41.2734},
        {"img-0004", 41.9400, 286.9612, "2024-01-09T22:14:11", 10.4917, 41.2887},
        {"img-0005", 41.9145, 286.9887, "2024-01-09T22:14:21", 10.4929, 41.2904},
        {"img-0006", 41.8891, 287.0162, "2024-01-09T22:14:32", 10.4935, 41.2916}
    };
    
    logMessage("Testing with and without tilt correction:", "gray");
    logMessage("", "gray");
    
    // Save current state
    bool originalState = m_mountTilt.enableCorrection;
    
    for (const TestCase &test : testCases) {
        logMessage(QString("=== %1 ===").arg(test.name), "green");
        
        // Test without correction
        m_mountTilt.enableCorrection = false;
        double ra_uncorrected, dec_uncorrected;
        convertAltAzToRaDec(test.inputAlt, test.inputAz, test.time, ra_uncorrected, dec_uncorrected);
        
        // Test with correction
        m_mountTilt.enableCorrection = true;
        double ra_corrected, dec_corrected;
        convertAltAzToRaDec(test.inputAlt, test.inputAz, test.time, ra_corrected, dec_corrected);
        
        // Calculate errors
        double error_uncorrected_ra = qAbs(ra_uncorrected - test.expectedRA);
        double error_uncorrected_dec = qAbs(dec_uncorrected - test.expectedDec);
        double error_corrected_ra = qAbs(ra_corrected - test.expectedRA);
        double error_corrected_dec = qAbs(dec_corrected - test.expectedDec);
        
        // Handle RA wraparound
        if (error_uncorrected_ra > 180) error_uncorrected_ra = 360 - error_uncorrected_ra;
        if (error_corrected_ra > 180) error_corrected_ra = 360 - error_corrected_ra;
        
        double total_error_uncorrected = sqrt(error_uncorrected_ra*error_uncorrected_ra + 
                                            error_uncorrected_dec*error_uncorrected_dec);
        double total_error_corrected = sqrt(error_corrected_ra*error_corrected_ra + 
                                          error_corrected_dec*error_corrected_dec);
        
        logMessage(QString("Input Alt/Az: %1°, %2°")
                      .arg(test.inputAlt, 0, 'f', 4).arg(test.inputAz, 0, 'f', 4), "gray");
        logMessage(QString("Expected RA/Dec: %1°, %2°")
                      .arg(test.expectedRA, 0, 'f', 4).arg(test.expectedDec, 0, 'f', 4), "orange");
        
        logMessage(QString("Without correction: RA=%1°, Dec=%2° (error: %3°)")
                      .arg(ra_uncorrected, 0, 'f', 4)
                      .arg(dec_uncorrected, 0, 'f', 4)
                      .arg(total_error_uncorrected, 0, 'f', 3), "red");
        
        logMessage(QString("With correction:    RA=%1°, Dec=%2° (error: %3°)")
                      .arg(ra_corrected, 0, 'f', 4)
                      .arg(dec_corrected, 0, 'f', 4)
                      .arg(total_error_corrected, 0, 'f', 3), 
                  (total_error_corrected < total_error_uncorrected) ? "green" : "orange");
        
        double improvement = total_error_uncorrected - total_error_corrected;
        logMessage(QString("Improvement: %1° (%2%)")
                      .arg(improvement, 0, 'f', 3)
                      .arg(improvement/total_error_uncorrected*100, 0, 'f', 1),
                  (improvement > 0) ? "green" : "red");
        logMessage("", "gray");
    }
    
    // Restore original state
    m_mountTilt.enableCorrection = originalState;
    
    logMessage("=== END MOUNT TILT CORRECTION TEST ===", "blue");
}

// Dynamic calibration that reads data directly from processed FITS files

void StellinaProcessor::calibrateFromProcessedFiles() {
    logMessage("=== DYNAMIC CALIBRATION FROM PROCESSED FITS FILES ===", "blue");
    
    QList<ProcessedImageData> imageData;
    
    // Scan calibrated directory for files with Stellina metadata
    QDir calibratedDir(m_calibratedDirectory);
    if (!calibratedDir.exists()) {
        logMessage("Calibrated directory not found. Run dark calibration first.", "red");
        return;
    }
    
    // Scan solved directory for solve-field results
    QDir solvedDir(m_plateSolvedDirectory);
    if (!solvedDir.exists()) {
        logMessage("Plate-solved directory not found. Run plate solving first.", "red");
        return;
    }
    
    // Get list of solved images
    QStringList solvedFiles = solvedDir.entryList(
        QStringList() << "plate_solved_*.fits", QDir::Files);
    
    if (solvedFiles.isEmpty()) {
        logMessage("No plate-solved images found for calibration", "red");
        return;
    }
    
    logMessage(QString("Found %1 plate-solved images to analyze").arg(solvedFiles.size()), "blue");
    
    // Process each solved image
    double sessionStart;
    bool sessionStartSet = false;
    
    for (const QString &solvedFile : solvedFiles) {
        ProcessedImageData data;
        data.filename = solvedFile;
        data.isValid = false;
        
        QString solvedPath = solvedDir.absoluteFilePath(solvedFile);
        
        // Extract image number from filename (e.g., "plate_solved_img-0042.fits" -> 42)
        QRegularExpression imgNumRegex(R"(img-(\d+))");
        QRegularExpressionMatch match = imgNumRegex.match(solvedFile);
        if (match.hasMatch()) {
            data.imageNumber = match.captured(1).toInt();
        } else {
            continue; // Skip if can't extract image number
        }
        
        // Read Stellina metadata from solved FITS file
        if (!readStellinaDataFromSolvedFits(solvedPath, data)) {
            continue; // Skip if can't read metadata
        }
        
        // Read solve-field results from WCS headers
        if (!readSolveFieldResults(solvedPath, data)) {
            continue; // Skip if can't read solve results
        }
        
        // Parse observation time
        data.obsTime = QDateTime::fromString(data.dateObs, "yyyy-MM-ddThh:mm:ss");
        if (!data.obsTime.isValid()) {
            // Try other formats
            QStringList timeFormats = {
                "yyyy-MM-ddThh:mm:ss.zzz",
                "yyyy-MM-dd hh:mm:ss"
            };
            for (const QString &format : timeFormats) {
                data.obsTime = QDateTime::fromString(data.dateObs, format);
                if (data.obsTime.isValid()) break;
            }
        }
        
        if (!data.obsTime.isValid()) {
            logMessage(QString("Invalid time format in %1: %2").arg(solvedFile).arg(data.dateObs), "orange");
            continue;
        }
        
        data.obsTime.setTimeSpec(Qt::UTC);
        
        // Set session start time from first valid image
        if (!sessionStartSet) {
            sessionStart = data.obsTime.toMSecsSinceEpoch() / 60000.0 ;
            sessionStartSet = true;
            data.minutesFromStart = 0.0;
        } else {
            data.minutesFromStart = data.obsTime.toMSecsSinceEpoch() - sessionStart;
        }
        
        data.isValid = true;
        imageData.append(data);
        
        if (m_debugMode) {
            logMessage(QString("Loaded: %1 (img-%2, t=%.1fmin)")
                          .arg(solvedFile).arg(data.imageNumber).arg(data.minutesFromStart), "gray");
        }
    }
    
    if (imageData.isEmpty()) {
        logMessage("No valid calibration data found in FITS files", "red");
        return;
    }
    
    // Sort by image number to ensure chronological order
    std::sort(imageData.begin(), imageData.end(), 
              [](const ProcessedImageData &a, const ProcessedImageData &b) {
                  return a.imageNumber < b.imageNumber;
              });
    
    logMessage(QString("Successfully loaded %1 calibration points spanning %.1f minutes")
                  .arg(imageData.size())
                  .arg(imageData.last().minutesFromStart), "green");
    
    // Analyze errors and calculate drift
    analyzeAndCalibrateFromData(imageData, sessionStart);
}

// Helper function to read Stellina data from solved FITS
bool StellinaProcessor::readStellinaDataFromSolvedFits(const QString &fitsPath, ProcessedImageData &data) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Read Stellina Alt/Az coordinates
    if (fits_read_key(fptr, TDOUBLE, "STELLALT", &data.stellinaAlt, nullptr, &status) != 0 ||
        fits_read_key(fptr, TDOUBLE, "STELLAZ", &data.stellinaAz, nullptr, &status) != 0) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Read predicted RA/Dec (from mount coordinate conversion)
    status = 0;
    if (fits_read_key(fptr, TDOUBLE, "STELLRA", &data.predictedRA, nullptr, &status) != 0 ||
        fits_read_key(fptr, TDOUBLE, "STELLDEC", &data.predictedDec, nullptr, &status) != 0) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Read observation time
    char dateobs[FLEN_VALUE];
    status = 0;
    if (fits_read_key(fptr, TSTRING, "DATE-OBS", dateobs, nullptr, &status) == 0) {
        data.dateObs = QString::fromLatin1(dateobs).trimmed().remove('\'').remove('"');
    } else {
        fits_close_file(fptr, &status);
        return false;
    }
    
    fits_close_file(fptr, &status);
    return true;
}

bool StellinaProcessor::readSolveFieldResults(const QString &fitsPath, ProcessedImageData &data) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Read WCS parameters
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &data.solvedRA, nullptr, &status) != 0) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    status = 0;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &data.solvedDec, nullptr, &status) != 0) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Read CD matrix for pixel scale
    double cd11, cd12, cd21, cd22;
    status = 0;
    if (fits_read_key(fptr, TDOUBLE, "CD1_1", &cd11, nullptr, &status) == 0 &&
        fits_read_key(fptr, TDOUBLE, "CD1_2", &cd12, nullptr, &status) == 0 &&
        fits_read_key(fptr, TDOUBLE, "CD2_1", &cd21, nullptr, &status) == 0 &&
        fits_read_key(fptr, TDOUBLE, "CD2_2", &cd22, nullptr, &status) == 0) {
        
        // Calculate pixel scale from CD matrix
        // FIXED: Add pixelScale member to ProcessedImageData
        data.pixelScale = sqrt(cd11*cd11 + cd12*cd12) * 3600.0; // arcsec/pixel
        data.hasValidWCS = true;
    }
    
    fits_close_file(fptr, &status);
    return data.hasValidWCS;
}
