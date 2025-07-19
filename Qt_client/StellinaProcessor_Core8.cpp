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

void StellinaProcessor::analyzeStackingCorrections(const QList<StackingCorrectionData> &stackingData, 
                                                  const double &sessionStart) {
    logMessage("", "gray");
    logMessage("ANALYZING STACKING CORRECTIONS (MOSAIC-AWARE):", "blue");
    
    const double PIXEL_SCALE_ARCSEC = 1.25; 
    const double PIXEL_SCALE_DEG = PIXEL_SCALE_ARCSEC / 3600.0;
    
    // Fix time span display bug
    double timeSpan = stackingData.last().minutesFromStart - stackingData.first().minutesFromStart;
    logMessage(QString("Time span: %.1f minutes").arg(timeSpan), "blue");
    
    // Use first image as reference
    double refCorrectionX = stackingData.first().correctionX;
    double refCorrectionY = stackingData.first().correctionY;
    
    logMessage(QString("Reference correction (img-%1): X=%2px, Y=%3px")
                  .arg(stackingData.first().imageNumber, 4, 10, QChar('0'))
                  .arg(refCorrectionX, 0, 'f', 1)
                  .arg(refCorrectionY, 0, 'f', 1), "blue");
    
    // Analyze correction patterns to detect mosaic behavior
    QList<double> xCorrections, yCorrections, timePoints;
    QMap<QString, int> patternCount;
    
    double maxCorrectionMagnitude = 0.0;
    double totalCorrectionMagnitude = 0.0;
    int largeJumps = 0;
    
    logMessage("", "gray");
    logMessage("Detecting correction patterns:", "blue");
    
    for (int i = 0; i < stackingData.size(); ++i) {
        const StackingCorrectionData &data = stackingData[i];
        
        double deltaX = data.correctionX - refCorrectionX;
        double deltaY = data.correctionY - refCorrectionY;
        double magnitude = sqrt(deltaX*deltaX + deltaY*deltaY);
        
        xCorrections.append(deltaX);
        yCorrections.append(deltaY);
        timePoints.append(data.minutesFromStart);
        
        totalCorrectionMagnitude += magnitude;
        maxCorrectionMagnitude = qMax(maxCorrectionMagnitude, magnitude);
        
        // Detect large jumps (likely mosaic position changes)
        if (magnitude > 100.0) { // 100 pixel threshold
            largeJumps++;
        }
        
        // Classify correction patterns
        QString pattern;
        if (qAbs(deltaX) < 50 && qAbs(deltaY) < 50) {
            pattern = "fine_tracking";
        } else if (qAbs(deltaY) > 500) {
            pattern = "mosaic_y_jump";
        } else if (qAbs(deltaX) > 200) {
            pattern = "mosaic_x_jump";
        } else {
            pattern = "medium_correction";
        }
        
        patternCount[pattern]++;
        
        // Show sample corrections
        if (i % 50 == 0 || i == stackingData.size() - 1) {
            logMessage(QString("img-%1 (t=%2min): dX=%3px dY=%4px mag=%5px pattern=%6")
                          .arg(data.imageNumber, 4, 10, QChar('0'))
                          .arg(data.minutesFromStart, 0, 'f', 1)
                          .arg(deltaX, 0, 'f', 1)
                          .arg(deltaY, 0, 'f', 1)
                          .arg(magnitude, 0, 'f', 1)
                          .arg(pattern), "gray");
        }
    }
    
    double avgCorrectionMagnitude = totalCorrectionMagnitude / stackingData.size();
    
    logMessage("", "gray");
    logMessage("CORRECTION PATTERN ANALYSIS:", "blue");
    logMessage(QString("Total images analyzed: %1").arg(stackingData.size()), "blue");
    logMessage(QString("Average correction magnitude: %1 pixels (%.2f arcsec)")
                  .arg(avgCorrectionMagnitude, 0, 'f', 1)
                  .arg(avgCorrectionMagnitude * PIXEL_SCALE_ARCSEC, 0, 'f', 1), "blue");
    logMessage(QString("Maximum correction magnitude: %1 pixels (%.2f arcsec)")
                  .arg(maxCorrectionMagnitude, 0, 'f', 1)
                  .arg(maxCorrectionMagnitude * PIXEL_SCALE_ARCSEC, 0, 'f', 1), "blue");
    logMessage(QString("Large jumps detected: %1").arg(largeJumps), "blue");
    
    logMessage("", "gray");
    logMessage("Pattern breakdown:", "blue");
    for (auto it = patternCount.begin(); it != patternCount.end(); ++it) {
        double percentage = (it.value() * 100.0) / stackingData.size();
        logMessage(QString("  %1: %2 images (%3%)")
                      .arg(it.key(), -15)
                      .arg(it.value(), 3)
                      .arg(percentage, 0, 'f', 1), "gray");
    }
    
    // Determine if this is mosaic or tracking data
    bool isMosaicData = (largeJumps > stackingData.size() * 0.1); // >10% large jumps
    bool hasMosaicYJumps = patternCount.contains("mosaic_y_jump") && patternCount["mosaic_y_jump"] > 0;
    
    logMessage("", "gray");
    if (isMosaicData || hasMosaicYJumps) {
        logMessage("DETECTED: MOSAIC OBSERVATION PATTERN", "orange");
        logMessage("This appears to be a mosaic observation with multiple pointings.", "orange");
        logMessage("Standard drift analysis is not applicable to mosaic data.", "orange");
        logMessage("", "gray");
        
        analyzeMosaicCorrections(stackingData, patternCount);
        
    } else {
        logMessage("DETECTED: SINGLE FIELD TRACKING", "green");
        logMessage("Proceeding with drift analysis for single field observation.", "green");
        logMessage("", "gray");
        
        performDriftAnalysis(stackingData, timePoints, xCorrections, yCorrections, sessionStart);
    }
}

void StellinaProcessor::analyzeMosaicCorrections(const QList<StackingCorrectionData> &stackingData,
                                                const QMap<QString, int> &patternCount) {
    logMessage("MOSAIC ANALYSIS:", "blue");
    
    // For mosaic data, analyze the fine tracking components only
    QList<double> fineTrackingX, fineTrackingY, fineTrackingTimes;
    
    double refX = stackingData.first().correctionX;
    double refY = stackingData.first().correctionY;
    
    for (const StackingCorrectionData &data : stackingData) {
        double deltaX = data.correctionX - refX;
        double deltaY = data.correctionY - refY;
        double magnitude = sqrt(deltaX*deltaX + deltaY*deltaY);
        
        // Only include fine tracking corrections (< 50 pixels)
        if (magnitude < 50.0) {
            fineTrackingX.append(deltaX);
            fineTrackingY.append(deltaY);
            fineTrackingTimes.append(data.minutesFromStart);
        }
    }
    
    if (fineTrackingX.size() > 10) {
        logMessage(QString("Analyzing fine tracking from %1 images (excluding mosaic jumps)")
                      .arg(fineTrackingX.size()), "blue");
        
        // Calculate RMS of fine tracking errors
        double rmsX = 0.0, rmsY = 0.0;
        for (int i = 0; i < fineTrackingX.size(); ++i) {
            rmsX += fineTrackingX[i] * fineTrackingX[i];
            rmsY += fineTrackingY[i] * fineTrackingY[i];
        }
        rmsX = sqrt(rmsX / fineTrackingX.size());
        rmsY = sqrt(rmsY / fineTrackingY.size());
        
        const double PIXEL_SCALE_ARCSEC = 1.25;
        logMessage(QString("Fine tracking RMS: X=%.1fpx (%.1f\"), Y=%.1fpx (%.1f\")")
                      .arg(rmsX, 0, 'f', 1).arg(rmsX * PIXEL_SCALE_ARCSEC, 0, 'f', 1)
                      .arg(rmsY, 0, 'f', 1).arg(rmsY * PIXEL_SCALE_ARCSEC, 0, 'f', 1), "green");
        
        if (rmsX < 10.0 && rmsY < 10.0) {
            logMessage("✓ EXCELLENT: Fine tracking errors < 10 pixels", "green");
            logMessage("Stellina's mount tracking is working well within mosaic pointings", "green");
        } else if (rmsX < 20.0 && rmsY < 20.0) {
            logMessage("✓ GOOD: Fine tracking errors < 20 pixels", "green");
        } else {
            logMessage("⚠ MODERATE: Fine tracking errors are larger than expected", "orange");
        }
    }
    
    logMessage("", "gray");
    logMessage("MOSAIC OBSERVATION RECOMMENDATIONS:", "blue");
    logMessage("• Mount drift correction is not applicable to mosaic data", "gray");
    logMessage("• Fine tracking appears adequate for mosaic observations", "gray");
    logMessage("• Large corrections are normal mosaic positioning adjustments", "gray");
    logMessage("• Consider analyzing individual mosaic tiles separately", "gray");
    
    // Don't apply drift correction for mosaic data
    logMessage("", "gray");
    logMessage("Drift correction NOT applied - mosaic pattern detected", "orange");
}

void StellinaProcessor::performDriftAnalysis(const QList<StackingCorrectionData> &stackingData,
                                           const QList<double> &timePoints,
                                           const QList<double> &xCorrections,
                                           const QList<double> &yCorrections,
                                           const double &sessionStart) {
    // Original drift analysis code (for single field observations)
    
    const double PIXEL_SCALE_DEG = 1.25 / 3600.0; // degrees per pixel
    
    // Convert corrections to angular errors
    QList<double> raErrors, decErrors;
    for (int i = 0; i < xCorrections.size(); ++i) {
        raErrors.append(xCorrections[i] * PIXEL_SCALE_DEG);
        decErrors.append(yCorrections[i] * PIXEL_SCALE_DEG);
    }
    
    // Perform linear regression
    int n = timePoints.size();
    double sumTime = 0.0, sumRA = 0.0, sumDec = 0.0;
    double sumTimeSquared = 0.0, sumTimeRA = 0.0, sumTimeDec = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double time = timePoints[i];
        double raErr = raErrors[i];
        double decErr = decErrors[i];
        
        sumTime += time;
        sumRA += raErr;
        sumDec += decErr;
        sumTimeSquared += time * time;
        sumTimeRA += time * raErr;
        sumTimeDec += time * decErr;
    }
    
    double raSlope = (n * sumTimeRA - sumTime * sumRA) / (n * sumTimeSquared - sumTime * sumTime);
    double raIntercept = (sumRA - raSlope * sumTime) / n;
    
    double decSlope = (n * sumTimeDec - sumTime * sumDec) / (n * sumTimeSquared - sumTime * sumTime);
    double decIntercept = (sumDec - decSlope * sumTime) / n;
    
    double raSlope_hourly = raSlope * 60.0;
    double decSlope_hourly = decSlope * 60.0;
    
    // Calculate correlation
    double raVariance = 0.0, timeVariance = 0.0, raTimeCovariance = 0.0;
    double raMean = sumRA / n;
    double timeMean = sumTime / n;
    
    for (int i = 0; i < n; ++i) {
        double deltaRA = raErrors[i] - raMean;
        double deltaTime = timePoints[i] - timeMean;
        raVariance += deltaRA * deltaRA;
        timeVariance += deltaTime * deltaTime;
        raTimeCovariance += deltaRA * deltaTime;
    }
    
    double raCorrelation = raTimeCovariance / sqrt(raVariance * timeVariance);
    double raR_squared = raCorrelation * raCorrelation;
    
    logMessage("LINEAR REGRESSION RESULTS:", "green");
    logMessage(QString("RA drift: %1° + %2°/hour * t")
                  .arg(raIntercept, 0, 'f', 4)
                  .arg(raSlope_hourly, 0, 'f', 3), "green");
    logMessage(QString("Dec drift: %1° + %2°/hour * t")
                  .arg(decIntercept, 0, 'f', 4)
                  .arg(decSlope_hourly, 0, 'f', 3), "green");
    logMessage(QString("RA regression R² = %1").arg(raR_squared, 0, 'f', 3), "blue");
    
    if (raR_squared > 0.5) {
        // Apply calibration for good correlations
        logMessage("", "gray");
        logMessage("APPLYING DRIFT CALIBRATION:", "green");
        
        m_mountTilt.initialRAOffset = raIntercept;
        m_mountTilt.initialDecOffset = decIntercept;
        m_mountTilt.driftRA = raSlope_hourly;
        m_mountTilt.driftDec = decSlope_hourly;
        m_mountTilt.sessionStart = sessionStart;
        m_mountTilt.enableCorrection = true;
        m_mountTilt.enableDriftCorrection = true;
        m_mountTilt.northTilt = 0.0;
        m_mountTilt.eastTilt = 0.0;
        
        updateTiltUI();
        saveMountTiltToSettings();
        
        logMessage("✓ Drift correction applied and active", "green");
    } else {
        logMessage("", "gray");
        logMessage("⚠ Poor correlation - drift correction not applied", "orange");
        logMessage("Consider using plate solving calibration method instead", "orange");
    }
}

// Replace the testSystematicOffsetCorrection function with this real data version
// Add to StellinaProcessor_Core.cpp

void StellinaProcessor::testSystematicOffsetCorrection() {
    logMessage("=== TESTING DRIFT CORRECTION WITH REAL FITS DATA ===", "blue");
    
    // Check if we have plate-solved files to work with
    QDir plateSolvedDir(m_plateSolvedDirectory);
    if (m_plateSolvedDirectory.isEmpty() || !plateSolvedDir.exists()) {
        logMessage("No plate-solved directory found. Run plate solving first.", "red");
        return;
    }
    
    // Get list of plate-solved FITS files
    QStringList plateSolvedFiles = plateSolvedDir.entryList(
        QStringList() << "plate_solved_*.fits" << "*solved*.fits", 
        QDir::Files);
    
    if (plateSolvedFiles.isEmpty()) {
        logMessage("No plate-solved FITS files found for testing", "red");
        return;
    }
    
    // Sort files and take a sample for testing (every 10th file to cover time span)
    plateSolvedFiles.sort();
    QStringList testFiles;
    for (int i = 0; i < plateSolvedFiles.size(); i += qMax(1, plateSolvedFiles.size() / 20)) {
        testFiles.append(plateSolvedFiles[i]);
        if (testFiles.size() >= 20) break;  // Limit to 20 test points
    }
    
    logMessage(QString("Testing with %1 files spanning the observation session").arg(testFiles.size()), "blue");
    logMessage("", "gray");
    
    // Display current correction parameters
    logMessage("Current drift correction parameters:", "blue");
    logMessage(QString("  Correction enabled: %1").arg(m_mountTilt.enableCorrection ? "Yes" : "No"), "blue");
    logMessage(QString("  Drift correction enabled: %1").arg(m_mountTilt.enableDriftCorrection ? "Yes" : "No"), "blue");
    
    if (m_mountTilt.enableCorrection && m_mountTilt.enableDriftCorrection) {
        logMessage(QString("  Initial RA offset: %1°").arg(m_mountTilt.initialRAOffset, 0, 'f', 4), "blue");
        logMessage(QString("  Initial Dec offset: %1°").arg(m_mountTilt.initialDecOffset, 0, 'f', 4), "blue");
        logMessage(QString("  RA drift rate: %1°/hour").arg(m_mountTilt.driftRA, 0, 'f', 3), "blue");
        logMessage(QString("  Dec drift rate: %1°/hour").arg(m_mountTilt.driftDec, 0, 'f', 3), "blue");
        logMessage(QString("  Session start: %1").arg(m_mountTilt.sessionStart), "blue");
    } else {
        logMessage("  Drift correction is disabled - enable it first", "red");
        return;
    }
    logMessage("", "gray");
    
    // Test each file
    QList<double> errorsBeforeCorrection, errorsAfterCorrection;
    QList<double> timePoints;
    QDateTime sessionStart;
    bool sessionStartSet = false;
    
    for (const QString &fileName : testFiles) {
        QString fitsPath = plateSolvedDir.absoluteFilePath(fileName);
        
        // Read coordinates from FITS headers
        fitsfile *fptr = nullptr;
        int status = 0;
        
        QByteArray pathBytes = fitsPath.toLocal8Bit();
        if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
            logMessage(QString("Failed to open FITS file: %1").arg(fileName), "orange");
            continue;
        }
        
        // Read Stellina Alt/Az coordinates
        double stellinaAlt, stellinaAz;
        if (fits_read_key(fptr, TDOUBLE, "STELLALT", &stellinaAlt, nullptr, &status) != 0 ||
            fits_read_key(fptr, TDOUBLE, "STELLAZ", &stellinaAz, nullptr, &status) != 0) {
            logMessage(QString("No Stellina Alt/Az in %1").arg(fileName), "orange");
            fits_close_file(fptr, &status);
            continue;
        }
        
        // Read observation time
        char dateobs[FLEN_VALUE];
        status = 0;
        if (fits_read_key(fptr, TSTRING, "DATE-OBS", dateobs, nullptr, &status) != 0) {
            logMessage(QString("No DATE-OBS in %1").arg(fileName), "orange");
            fits_close_file(fptr, &status);
            continue;
        }
        
        QString dateObsStr = QString::fromLatin1(dateobs).trimmed().remove('\'').remove('"');
        
        // Read solved RA/Dec (from WCS headers - this is the "ground truth")
        double solvedRA, solvedDec;
        status = 0;
        if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &solvedRA, nullptr, &status) != 0 ||
            fits_read_key(fptr, TDOUBLE, "CRVAL2", &solvedDec, nullptr, &status) != 0) {
            logMessage(QString("No WCS RA/Dec in %1").arg(fileName), "orange");
            fits_close_file(fptr, &status);
            continue;
        }
        
        fits_close_file(fptr, &status);
        
        // Parse observation time for elapsed time calculation
        QDateTime obsTime = QDateTime::fromString(dateObsStr, "yyyy-MM-ddThh:mm:ss");
        if (!obsTime.isValid()) {
            // Try other formats
            QStringList timeFormats = {
                "yyyy-MM-ddThh:mm:ss.zzz",
                "yyyy-MM-dd hh:mm:ss"
            };
            for (const QString &format : timeFormats) {
                obsTime = QDateTime::fromString(dateObsStr, format);
                if (obsTime.isValid()) break;
            }
        }
        
        if (!obsTime.isValid()) {
            logMessage(QString("Invalid time format in %1: %2").arg(fileName).arg(dateObsStr), "orange");
            continue;
        }
        
        obsTime.setTimeSpec(Qt::UTC);
        
        // Set session start from first valid file
        if (!sessionStartSet) {
            sessionStart = obsTime;
            sessionStartSet = true;
        }
        
        double minutesFromStart = sessionStart.msecsTo(obsTime) / 60000.0;
        timePoints.append(minutesFromStart);
        
        // Test WITHOUT correction (temporarily disable)
        bool originalEnabled = m_mountTilt.enableCorrection;
        bool originalDriftEnabled = m_mountTilt.enableDriftCorrection;
        m_mountTilt.enableCorrection = false;
        m_mountTilt.enableDriftCorrection = false;
        
        double ra_uncorrected, dec_uncorrected;
        if (!convertAltAzToRaDec(stellinaAlt, stellinaAz, dateObsStr, ra_uncorrected, dec_uncorrected)) {
            logMessage(QString("Failed to convert coordinates for %1").arg(fileName), "orange");
            continue;
        }
        
        // Test WITH correction (re-enable)
        m_mountTilt.enableCorrection = originalEnabled;
        m_mountTilt.enableDriftCorrection = originalDriftEnabled;
        
        double ra_corrected, dec_corrected;
        if (!convertAltAzToRaDec(stellinaAlt, stellinaAz, dateObsStr, ra_corrected, dec_corrected)) {
            logMessage(QString("Failed to convert corrected coordinates for %1").arg(fileName), "orange");
            continue;
        }
        
        // Calculate errors against solve-field solution (ground truth)
        double error_uncorrected_ra = ra_uncorrected - solvedRA;
        double error_uncorrected_dec = dec_uncorrected - solvedDec;
        double error_corrected_ra = ra_corrected - solvedRA;
        double error_corrected_dec = dec_corrected - solvedDec;
        
        // Handle RA wraparound
        if (error_uncorrected_ra > 180) error_uncorrected_ra -= 360;
        if (error_uncorrected_ra < -180) error_uncorrected_ra += 360;
        if (error_corrected_ra > 180) error_corrected_ra -= 360;
        if (error_corrected_ra < -180) error_corrected_ra += 360;
        
        double total_error_uncorrected = sqrt(error_uncorrected_ra*error_uncorrected_ra + 
                                            error_uncorrected_dec*error_uncorrected_dec);
        double total_error_corrected = sqrt(error_corrected_ra*error_corrected_ra + 
                                          error_corrected_dec*error_corrected_dec);
        
        errorsBeforeCorrection.append(total_error_uncorrected);
        errorsAfterCorrection.append(total_error_corrected);
        
        // Extract image number for display
        QString baseName = QFileInfo(fileName).baseName();
        QRegularExpression imgNumRegex(R"(img-(\d+))");
        QRegularExpressionMatch match = imgNumRegex.match(baseName);
        QString imgNum = match.hasMatch() ? match.captured(1) : "???";
        
        // Log results
        logMessage(QString("Img-%1 (t=%2min):")
                      .arg(imgNum, 3)
                      .arg(minutesFromStart, 0, 'f', 1), "green");
        logMessage(QString("  Alt/Az: %1°, %2°")
                      .arg(stellinaAlt, 0, 'f', 4).arg(stellinaAz, 0, 'f', 4), "gray");
        logMessage(QString("  Solved: RA=%1°, Dec=%2° (ground truth)")
                      .arg(solvedRA, 0, 'f', 4).arg(solvedDec, 0, 'f', 4), "orange");
        logMessage(QString("  Before correction: RA=%1°, Dec=%2° (error: %3°)")
                      .arg(ra_uncorrected, 0, 'f', 4)
                      .arg(dec_uncorrected, 0, 'f', 4)
                      .arg(total_error_uncorrected, 0, 'f', 3), "red");
        logMessage(QString("  After correction:  RA=%1°, Dec=%2° (error: %3°)")
                      .arg(ra_corrected, 0, 'f', 4)
                      .arg(dec_corrected, 0, 'f', 4)
                      .arg(total_error_corrected, 0, 'f', 3), 
                  (total_error_corrected < total_error_uncorrected) ? "green" : "orange");
        
        double improvement = total_error_uncorrected - total_error_corrected;
        logMessage(QString("  Improvement: %1° (%2%)")
                      .arg(improvement, 0, 'f', 3)
                      .arg(improvement/total_error_uncorrected*100, 0, 'f', 1),
                  (improvement > 0) ? "green" : "red");
        logMessage("", "gray");
    }
    
    if (errorsBeforeCorrection.isEmpty()) {
        logMessage("No valid test data found", "red");
        return;
    }
    
    // Calculate overall statistics
    double avgErrorBefore = 0.0, avgErrorAfter = 0.0;
    double maxErrorBefore = 0.0, maxErrorAfter = 0.0;
    double minErrorBefore = 1000.0, minErrorAfter = 1000.0;
    
    for (int i = 0; i < errorsBeforeCorrection.size(); ++i) {
        avgErrorBefore += errorsBeforeCorrection[i];
        avgErrorAfter += errorsAfterCorrection[i];
        maxErrorBefore = qMax(maxErrorBefore, errorsBeforeCorrection[i]);
        maxErrorAfter = qMax(maxErrorAfter, errorsAfterCorrection[i]);
        minErrorBefore = qMin(minErrorBefore, errorsBeforeCorrection[i]);
        minErrorAfter = qMin(minErrorAfter, errorsAfterCorrection[i]);
    }
    
    avgErrorBefore /= errorsBeforeCorrection.size();
    avgErrorAfter /= errorsAfterCorrection.size();
    
    double overallImprovement = avgErrorBefore - avgErrorAfter;
    double improvementPercent = (overallImprovement / avgErrorBefore) * 100.0;
    
    // Calculate RMS error
    double rmsBefore = 0.0, rmsAfter = 0.0;
    for (int i = 0; i < errorsBeforeCorrection.size(); ++i) {
        rmsBefore += errorsBeforeCorrection[i] * errorsBeforeCorrection[i];
        rmsAfter += errorsAfterCorrection[i] * errorsAfterCorrection[i];
    }
    rmsBefore = sqrt(rmsBefore / errorsBeforeCorrection.size());
    rmsAfter = sqrt(rmsAfter / errorsAfterCorrection.size());
    
    logMessage("=== OVERALL TEST RESULTS ===", "blue");
    logMessage(QString("Test data points: %1").arg(errorsBeforeCorrection.size()), "blue");
    logMessage(QString("Time span tested: %.1f minutes").arg(timePoints.last() - timePoints.first()), "blue");
    logMessage("", "gray");
    
    logMessage("Error Statistics:", "blue");
    logMessage(QString("  Before correction: avg=%.3f°, rms=%.3f°, range=%.3f°-%.3f°")
                  .arg(avgErrorBefore).arg(rmsBefore).arg(minErrorBefore).arg(maxErrorBefore), "orange");
    logMessage(QString("  After correction:  avg=%.3f°, rms=%.3f°, range=%.3f°-%.3f°")
                  .arg(avgErrorAfter).arg(rmsAfter).arg(minErrorAfter).arg(maxErrorAfter), "green");
    logMessage(QString("  Improvement: %.3f° (%.1f%% reduction)")
                  .arg(overallImprovement).arg(improvementPercent), 
              (overallImprovement > 0) ? "green" : "red");
    
    // Performance assessment
    logMessage("", "gray");
    if (avgErrorAfter < 0.5) {
        logMessage("✓ EXCELLENT: Average error now < 0.5° - drift correction working very well!", "green");
    } else if (avgErrorAfter < 1.0) {
        logMessage("✓ GOOD: Average error now < 1.0° - significant improvement", "green");
    } else if (avgErrorAfter < 2.0) {
        logMessage("✓ MODERATE: Average error now < 2.0° - some improvement", "orange");
    } else if (overallImprovement > 0) {
        logMessage("⚠ PARTIAL: Some improvement but errors still large", "orange");
    } else {
        logMessage("✗ PROBLEM: Drift correction not working - check parameters", "red");
    }
    
    // Check for residual drift
    if (errorsAfterCorrection.size() >= 5) {
        // Simple check: compare first few and last few errors
        double avgEarlyError = 0.0, avgLateError = 0.0;
        int nEarly = qMin(3, errorsAfterCorrection.size() / 3);
        int nLate = qMin(3, errorsAfterCorrection.size() / 3);
        
        for (int i = 0; i < nEarly; ++i) {
            avgEarlyError += errorsAfterCorrection[i];
        }
        avgEarlyError /= nEarly;
        
        for (int i = errorsAfterCorrection.size() - nLate; i < errorsAfterCorrection.size(); ++i) {
            avgLateError += errorsAfterCorrection[i];
        }
        avgLateError /= nLate;
        
        double residualDrift = avgLateError - avgEarlyError;
        logMessage(QString("Residual drift check: early=%.3f°, late=%.3f°, drift=%.3f°")
                      .arg(avgEarlyError).arg(avgLateError).arg(residualDrift), "blue");
        
        if (qAbs(residualDrift) < 0.5) {
            logMessage("✓ Residual drift is small - linear model fits well", "green");
        } else {
            logMessage("⚠ Significant residual drift - may need higher-order correction", "orange");
        }
    }
    
    logMessage("=== END DRIFT CORRECTION TEST ===", "blue");
}
