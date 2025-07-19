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

// Helper function to read solve-field results from WCS headers
bool StellinaProcessor::readSolveFieldResults(const QString &fitsPath, ProcessedImageData &data) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Read WCS reference point (CRVAL1/CRVAL2 are RA/Dec of reference pixel)
    double crval1, crval2;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &crval1, nullptr, &status) == 0 &&
        fits_read_key(fptr, TDOUBLE, "CRVAL2", &crval2, nullptr, &status) == 0) {
        data.solvedRA = crval1;
        data.solvedDec = crval2;
        fits_close_file(fptr, &status);
        return true;
    }
    
    // If WCS not available, try to parse from solve-field comments in header
    // Look for "Field center: (RA,Dec) = (xxx, yyy)" in HISTORY cards
    status = 0;
    int nkeys;
    fits_get_hdrspace(fptr, &nkeys, nullptr, &status);
    
    for (int i = 1; i <= nkeys; ++i) {
        char card[FLEN_CARD];
        if (fits_read_record(fptr, i, card, &status) == 0) {
            QString cardStr = QString::fromLatin1(card);
            
            // Look for solve-field result pattern
            QRegularExpression fieldCenterRegex(R"(Field center.*\(RA,Dec\)\s*=\s*\(([\d.]+),\s*([\d.]+)\))");
            QRegularExpressionMatch match = fieldCenterRegex.match(cardStr);
            if (match.hasMatch()) {
                data.solvedRA = match.captured(1).toDouble();
                data.solvedDec = match.captured(2).toDouble();
                fits_close_file(fptr, &status);
                return true;
            }
        }
        status = 0; // Reset for next iteration
    }
    
    fits_close_file(fptr, &status);
    return false;
}
// Update the analyzeAndCalibrateFromData function to use linear regression
// Replace the existing function in StellinaProcessor_Core.cpp

void StellinaProcessor::analyzeAndCalibrateFromData(const QList<ProcessedImageData> &imageData, 
                                                   const double &sessionStart) {
    logMessage("", "gray");
    logMessage("ANALYZING MOUNT ERRORS WITH LINEAR REGRESSION:", "blue");
    
    // Calculate errors for each image
    QList<double> timePoints, raErrors, decErrors, totalErrors;
    double maxError = 0.0;
    double totalErrorSum = 0.0;
    
    for (const ProcessedImageData &data : imageData) {
        double raError = data.predictedRA - data.solvedRA;
        double decError = data.predictedDec - data.solvedDec;
        
        // Handle RA wraparound
        if (raError > 180) raError -= 360;
        if (raError < -180) raError += 360;
        
        double totalError = sqrt(raError*raError + decError*decError);
        
        timePoints.append(data.minutesFromStart);
        raErrors.append(raError);
        decErrors.append(decError);
        totalErrors.append(totalError);
        
        maxError = qMax(maxError, totalError);
        totalErrorSum += totalError;
        
        QString errorLevel = (totalError < 1.0) ? "green" : 
                           (totalError < 2.0) ? "orange" : "red";
        
        logMessage(QString("Img-%1 (t=%2min): Error RA=%3° Dec=%4° Total=%5°")
                      .arg(data.imageNumber, 3)
                      .arg(data.minutesFromStart, 0, 'f', 1)
                      .arg(raError, 0, 'f', 3)
                      .arg(decError, 0, 'f', 3)
                      .arg(totalError, 0, 'f', 3), errorLevel);
    }
    
    int n = imageData.size();
    double avgError = totalErrorSum / n;
    
    logMessage("", "gray");
    logMessage("ERROR STATISTICS:", "blue");
    logMessage(QString("Number of data points: %1").arg(n), "blue");
    logMessage(QString("Time span: %.1f minutes").arg(timePoints.last() - timePoints.first()), "blue");
    logMessage(QString("Average total error: %1°").arg(avgError, 0, 'f', 3), "orange");
    logMessage(QString("Maximum error: %1°").arg(maxError, 0, 'f', 3), "orange");
    
    // Perform linear regression on RA and Dec errors vs time
    logMessage("", "gray");
    logMessage("LINEAR REGRESSION ANALYSIS:", "blue");
    
    // Calculate linear regression: error = intercept + slope * time
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
    
    // RA regression: raError = raIntercept + raSlope * time
    double raSlope = (n * sumTimeRA - sumTime * sumRA) / (n * sumTimeSquared - sumTime * sumTime);
    double raIntercept = (sumRA - raSlope * sumTime) / n;
    
    // Dec regression: decError = decIntercept + decSlope * time  
    double decSlope = (n * sumTimeDec - sumTime * sumDec) / (n * sumTimeSquared - sumTime * sumTime);
    double decIntercept = (sumDec - decSlope * sumTime) / n;
    
    // Convert slopes from degrees/minute to degrees/hour
    double raSlope_hourly = raSlope * 60.0;
    double decSlope_hourly = decSlope * 60.0;
    
    logMessage(QString("RA error regression: %1° + %2°/min * t")
                  .arg(raIntercept, 0, 'f', 3)
                  .arg(raSlope, 0, 'f', 4), "blue");
    logMessage(QString("                     %1° + %2°/hour * t")
                  .arg(raIntercept, 0, 'f', 3)
                  .arg(raSlope_hourly, 0, 'f', 3), "blue");
    
    logMessage(QString("Dec error regression: %1° + %2°/min * t")
                  .arg(decIntercept, 0, 'f', 3)
                  .arg(decSlope, 0, 'f', 4), "blue");
    logMessage(QString("                      %1° + %2°/hour * t")
                  .arg(decIntercept, 0, 'f', 3)
                  .arg(decSlope_hourly, 0, 'f', 3), "blue");
    
    // Calculate correlation coefficients to assess fit quality
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
    
    logMessage(QString("RA regression R² = %1 (correlation = %2)")
                  .arg(raR_squared, 0, 'f', 3)
                  .arg(raCorrelation, 0, 'f', 3), "blue");
    
    // Calculate similar for Dec
    double decVariance = 0.0, decTimeCovariance = 0.0;
    double decMean = sumDec / n;
    
    for (int i = 0; i < n; ++i) {
        double deltaDec = decErrors[i] - decMean;
        double deltaTime = timePoints[i] - timeMean;
        decVariance += deltaDec * deltaDec;
        decTimeCovariance += deltaDec * deltaTime;
    }
    
    double decCorrelation = decTimeCovariance / sqrt(decVariance * timeVariance);
    double decR_squared = decCorrelation * decCorrelation;
    
    logMessage(QString("Dec regression R² = %1 (correlation = %2)")
                  .arg(decR_squared, 0, 'f', 3)
                  .arg(decCorrelation, 0, 'f', 3), "blue");
    
    // Store the linear drift parameters
    logMessage("", "gray");
    logMessage("APPLYING TIME-DEPENDENT DRIFT CORRECTION:", "green");
    
    m_mountTilt.initialRAOffset = raIntercept;     // RA error at t=0
    m_mountTilt.initialDecOffset = decIntercept;   // Dec error at t=0
    m_mountTilt.driftRA = raSlope_hourly;          // RA drift rate (degrees/hour)
    m_mountTilt.driftDec = decSlope_hourly;        // Dec drift rate (degrees/hour)
    m_mountTilt.sessionStart = sessionStart;      // Session start time for t=0 reference
    m_mountTilt.enableCorrection = true;
    m_mountTilt.enableDriftCorrection = true;
    
    // Clear geometric tilt parameters (not using them)
    m_mountTilt.northTilt = 0.0;
    m_mountTilt.eastTilt = 0.0;
    m_mountTilt.systematicRAOffset = 0.0;   // Not using constant offset
    m_mountTilt.systematicDecOffset = 0.0;
    
    logMessage(QString("Initial RA offset (t=0): %1°").arg(m_mountTilt.initialRAOffset, 0, 'f', 4), "green");
    logMessage(QString("Initial Dec offset (t=0): %1°").arg(m_mountTilt.initialDecOffset, 0, 'f', 4), "green");
    logMessage(QString("RA drift rate: %1°/hour").arg(m_mountTilt.driftRA, 0, 'f', 3), "green");
    logMessage(QString("Dec drift rate: %1°/hour").arg(m_mountTilt.driftDec, 0, 'f', 3), "green");
    
    // Estimate improvement using linear model
    double totalPredictedError = 0.0;
    for (int i = 0; i < n; ++i) {
        double time_hours = timePoints[i] / 60.0;
        double predictedRAError = raIntercept + raSlope_hourly * time_hours;
        double predictedDecError = decIntercept + decSlope_hourly * time_hours;
        double correctedRAError = raErrors[i] - predictedRAError;
        double correctedDecError = decErrors[i] - predictedDecError;
        double correctedTotalError = sqrt(correctedRAError*correctedRAError + 
                                         correctedDecError*correctedDecError);
        totalPredictedError += correctedTotalError;
    }
    
    double avgPredictedError = totalPredictedError / n;
    double improvement = avgError - avgPredictedError;
    
    logMessage("", "gray");
    logMessage("EXPECTED IMPROVEMENT:", "green");
    logMessage(QString("Current average error: %1°").arg(avgError, 0, 'f', 3), "orange");
    logMessage(QString("Predicted residual error: %1°").arg(avgPredictedError, 0, 'f', 3), "green");
    logMessage(QString("Expected improvement: %1° (%2%)")
                  .arg(improvement, 0, 'f', 3)
                  .arg(improvement/avgError*100, 0, 'f', 1), "green");
    
    if (raR_squared > 0.8) {
        logMessage("✓ EXCELLENT: Strong linear correlation in RA drift", "green");
    } else if (raR_squared > 0.5) {
        logMessage("✓ GOOD: Moderate linear correlation in RA drift", "green");
    } else {
        logMessage("⚠ WARNING: Weak linear correlation - drift may be more complex", "orange");
    }
    
    if (avgPredictedError < 1.0) {
        logMessage("✓ Expected residual error < 1.0° - excellent correction!", "green");
    } else if (avgPredictedError < 2.0) {
        logMessage("✓ Expected residual error < 2.0° - good correction", "green");
    } else {
        logMessage("⚠ Expected residual error still high - may need additional correction", "orange");
    }
    
    // Update UI and save
    if (m_enableTiltCorrectionCheck && m_enableDriftCorrectionCheck) {
        m_enableTiltCorrectionCheck->setChecked(true);
        m_enableDriftCorrectionCheck->setChecked(true);
        m_northTiltSpin->setValue(0.0);
        m_eastTiltSpin->setValue(0.0);
        updateTiltUI();
    }
    
    saveMountTiltToSettings();
    
    logMessage("=== LINEAR REGRESSION CALIBRATION COMPLETE ===", "blue");
}

bool StellinaProcessor::convertAltAzToRaDec(double alt, double az, const QString &dateObs, double &ra, double &dec) {
  double observer_lat, observer_lon, jd, lst, ha;
  return convertAltAzToRaDecExt(alt, az, dateObs, ra, dec, observer_lat, observer_lon, jd, lst, ha);

}

// Update the coordinate conversion to apply time-dependent correction:
bool StellinaProcessor::convertAltAzToRaDecExt(double alt, double az, const QString &dateObs,
					       double &ra, double &dec, double &observer_lat, double &observer_lon,
					       double &jd, double &lst, double &ha) {
    // Apply mount tilt correction first (currently just passes through)
    double correctedAlt, correctedAz;
    applyMountTiltCorrection(correctedAlt, correctedAz, alt, az);
    
    // Parse observer location from settings
    QStringList locationParts = m_observerLocation.split(',');
    observer_lat = 51.5074;
    observer_lon = -0.1278;
    
    if (locationParts.size() >= 2) {
        bool ok1, ok2;
        double lat = locationParts[0].trimmed().toDouble(&ok1);
        double lon = locationParts[1].trimmed().toDouble(&ok2);
        if (ok1 && ok2) {
            observer_lat = lat;
            observer_lon = lon;
        }
    }
    
    // Parse observation time
    QDateTime obsTime;
    if (!dateObs.isEmpty()) {
        QStringList formats = {
            "yyyy-MM-ddThh:mm:ss.zzz",
            "yyyy-MM-ddThh:mm:ss.zzzZ", 
            "yyyy-MM-ddThh:mm:ss",
            "yyyy-MM-ddThh:mm:ssZ",
            "yyyy-MM-dd hh:mm:ss.zzz",
            "yyyy-MM-dd hh:mm:ss",
            "yyyy-MM-dd"
        };
        
        for (const QString &format : formats) {
            obsTime = QDateTime::fromString(dateObs, format);
            if (obsTime.isValid()) {
                obsTime.setTimeSpec(Qt::UTC);
                break;
            }
        }
    }
    
    if (!obsTime.isValid()) {
        logMessage(QString("ERROR: Could not parse observation time '%1'").arg(dateObs), "red");
        return false;
    }

    // Calculate Julian Date and LST
    jd = CoordinateUtils::computeJulianDay(obsTime.date().year(),
                           obsTime.date().month(),
                           obsTime.date().day(),
                           obsTime.time().hour(),
                           obsTime.time().minute(),
                           obsTime.time().second());
    
    lst = calculateLST_HighPrecision(jd, observer_lon);
    // Calculate RA/Dec from Alt/Az using the CoordinateUtils class
    // Convert horizontal to equatorial
    auto [raNow, decNow, haNow] = CoordinateUtils::altAzToRaDec(correctedAlt, correctedAz, observer_lat, observer_lon, lst);
    
    // Convert current epoch to J2000
    auto [ra2000, dec2000] = CoordinateUtils::jNowToJ2000(raNow, decNow);

    ra = raNow; // FIXME!
    dec = decNow;
    ha = haNow;
    
    // Apply time-dependent drift correction if enabled
    if (m_mountTilt.enableCorrection && m_mountTilt.enableDriftCorrection && 
        m_mountTilt.sessionStart != 0) {
        
        double originalRA = ra, originalDec = dec;
        
        // Calculate time elapsed since session start (in hours)
        double elapsedHours = (obsTime.toMSecsSinceEpoch() - m_mountTilt.sessionStart) / 3600000.0;
        
        // Apply linear drift correction: correction = initial_offset + drift_rate * elapsed_time
        double raCorrection = -(m_mountTilt.initialRAOffset + m_mountTilt.driftRA * elapsedHours);
        double decCorrection = -(m_mountTilt.initialDecOffset + m_mountTilt.driftDec * elapsedHours);
        
        ra += raCorrection;
        dec += decCorrection;
        
        // Normalize RA to [0, 360) range
        while (ra < 0) ra += 360.0;
        while (ra >= 360.0) ra -= 360.0;
        
        // Clamp Dec to valid range
        if (dec > 90.0) dec = 90.0;
        if (dec < -90.0) dec = -90.0;
        
        if (m_debugMode) {
            logMessage(QString("Time-dependent drift correction:")
                          .arg(elapsedHours, 0, 'f', 3), "blue");
            logMessage(QString("  Elapsed time: %1 hours").arg(elapsedHours, 0, 'f', 3), "gray");
            logMessage(QString("  RA correction: %1°").arg(raCorrection, 0, 'f', 3), "gray");
            logMessage(QString("  Dec correction: %1°").arg(decCorrection, 0, 'f', 3), "gray");
            logMessage(QString("  Result: RA %1°→%2°, Dec %3°→%4°")
                          .arg(originalRA, 0, 'f', 3).arg(ra, 0, 'f', 3)
                          .arg(originalDec, 0, 'f', 3).arg(dec, 0, 'f', 3), "blue");
        }
    }
    
    if (m_debugMode) {
        logMessage(QString("Final coordinates: RA=%1°, Dec=%2°")
                      .arg(ra, 0, 'f', 6).arg(dec, 0, 'f', 6), "green");
    }
    
    return true;
}

// Update the coordinate conversion to apply time-dependent correction:
bool StellinaProcessor::convertRaDecToAltAzExt(double ra, double dec, const QString &dateObs,
					       double &alt, double &az, double &observer_lat, double &observer_lon,
					       double &jd, double &lst, double &ha) {
    // Parse observer location from settings
    QStringList locationParts = m_observerLocation.split(',');
    observer_lat = 51.5074;
    observer_lon = -0.1278;
    
    if (locationParts.size() >= 2) {
        bool ok1, ok2;
        double lat = locationParts[0].trimmed().toDouble(&ok1);
        double lon = locationParts[1].trimmed().toDouble(&ok2);
        if (ok1 && ok2) {
            observer_lat = lat;
            observer_lon = lon;
        }
    }
    
    // Parse observation time
    QDateTime obsTime;
    if (!dateObs.isEmpty()) {
        QStringList formats = {
            "yyyy-MM-ddThh:mm:ss.zzz",
            "yyyy-MM-ddThh:mm:ss.zzzZ", 
            "yyyy-MM-ddThh:mm:ss",
            "yyyy-MM-ddThh:mm:ssZ",
            "yyyy-MM-dd hh:mm:ss.zzz",
            "yyyy-MM-dd hh:mm:ss",
            "yyyy-MM-dd"
        };
        
        for (const QString &format : formats) {
            obsTime = QDateTime::fromString(dateObs, format);
            if (obsTime.isValid()) {
                obsTime.setTimeSpec(Qt::UTC);
                break;
            }
        }
    }
    
    if (!obsTime.isValid()) {
        logMessage(QString("ERROR: Could not parse observation time '%1'").arg(dateObs), "red");
        return false;
    }

    // Calculate Julian Date and LST
    jd = CoordinateUtils::computeJulianDay(obsTime.date().year(),
                           obsTime.date().month(),
                           obsTime.date().day(),
                           obsTime.time().hour(),
                           obsTime.time().minute(),
                           obsTime.time().second());
    
    lst = calculateLST_HighPrecision(jd, observer_lon);
    // Calculate RA/Dec from Alt/Az using the CoordinateUtils class
    // Convert horizontal to equatorial
    auto [altNow, azNow, haNow] = CoordinateUtils::raDecToAltAz(ra/15.0, dec, observer_lat, observer_lon, lst);

    alt = altNow;
    az = azNow;
    ha = haNow;
    
    return true;
}

// Update the save/load settings functions:
void StellinaProcessor::saveMountTiltToSettings() {
    QSettings settings;
    
    settings.setValue("mountTilt/northTilt", m_mountTilt.northTilt);
    settings.setValue("mountTilt/eastTilt", m_mountTilt.eastTilt);
    settings.setValue("mountTilt/driftRA", m_mountTilt.driftRA);
    settings.setValue("mountTilt/driftDec", m_mountTilt.driftDec);
    settings.setValue("mountTilt/systematicRAOffset", m_mountTilt.systematicRAOffset);
    settings.setValue("mountTilt/systematicDecOffset", m_mountTilt.systematicDecOffset);
    settings.setValue("mountTilt/initialRAOffset", m_mountTilt.initialRAOffset);
    settings.setValue("mountTilt/initialDecOffset", m_mountTilt.initialDecOffset);
    settings.setValue("mountTilt/enableCorrection", m_mountTilt.enableCorrection);
    settings.setValue("mountTilt/enableDriftCorrection", m_mountTilt.enableDriftCorrection);
    settings.setValue("mountTilt/sessionStart", m_mountTilt.sessionStart);
    
    if (m_debugMode) {
        logMessage("Mount drift correction parameters saved to settings", "gray");
    }
}

bool StellinaProcessor::loadMountTiltFromSettings() {
    QSettings settings;
    
    m_mountTilt.northTilt = settings.value("mountTilt/northTilt", 0.0).toDouble();
    m_mountTilt.eastTilt = settings.value("mountTilt/eastTilt", 0.0).toDouble();
    m_mountTilt.driftRA = settings.value("mountTilt/driftRA", 0.0).toDouble();
    m_mountTilt.driftDec = settings.value("mountTilt/driftDec", 0.0).toDouble();
    m_mountTilt.systematicRAOffset = settings.value("mountTilt/systematicRAOffset", 0.0).toDouble();
    m_mountTilt.systematicDecOffset = settings.value("mountTilt/systematicDecOffset", 0.0).toDouble();
    m_mountTilt.initialRAOffset = settings.value("mountTilt/initialRAOffset", 0.0).toDouble();
    m_mountTilt.initialDecOffset = settings.value("mountTilt/initialDecOffset", 0.0).toDouble();
    m_mountTilt.enableCorrection = settings.value("mountTilt/enableCorrection", false).toBool();
    m_mountTilt.enableDriftCorrection = settings.value("mountTilt/enableDriftCorrection", false).toBool();
    m_mountTilt.sessionStart = settings.value("mountTilt/sessionStart").toDateTime().toMSecsSinceEpoch();
    
    if (m_mountTilt.enableCorrection && m_mountTilt.enableDriftCorrection) {
        logMessage(QString("Loaded mount drift correction: initial RA=%1°, drift=%2°/h")
                      .arg(m_mountTilt.initialRAOffset, 0, 'f', 4)
                      .arg(m_mountTilt.driftRA, 0, 'f', 3), "blue");
    }
    
    return m_mountTilt.enableCorrection;
}

// Modified processImageDarkCalibration to handle bayer patterns
bool StellinaProcessor::processImageDarkCalibration(const QString &lightFrame) {
    m_currentTaskLabel->setText("Dark calibration...");
    
    // Find the corresponding StellinaImageData
    StellinaImageData imageData;
    bool found = false;
    
    for (int i = 0; i < m_stellinaImageData.size(); ++i) {
        if (m_stellinaImageData[i].currentFitsPath == lightFrame || 
            m_stellinaImageData[i].originalFitsPath == lightFrame) {
            imageData = m_stellinaImageData[i];
            found = true;
            break;
        }
    }
    
    if (!found) {
        logMessage(QString("No metadata found for light frame: %1").arg(QFileInfo(lightFrame).fileName()), "red");
        return false;
    }
    
    // Log reversed image information
    if (imageData.isReversedImage) {
        logMessage(QString("Processing reversed stellina image: %1 (Bayer: %2)")
                      .arg(QFileInfo(lightFrame).fileName())
                      .arg(imageData.bayerPattern), "blue");
    }
    
    // Use metadata from imageData
    int lightExposure = imageData.exposureSeconds;
    int lightTemperatureK = imageData.temperatureKelvin;
    QString lightBinning = imageData.binning;
    QString lightBayerPattern = imageData.bayerPattern;
    
    if (lightExposure <= 0) {
        logMessage("Invalid exposure time in image metadata", "red");
        return false;
    }
    
    // Find matching dark frames with bayer pattern consideration
    QStringList matchingDarks = findAllMatchingDarkFrames(lightExposure, lightTemperatureK, 
                                                         lightBinning, lightBayerPattern);
    
    if (matchingDarks.isEmpty()) {
        int temperatureC = lightTemperatureK - 273;
        logMessage(QString("No matching dark frames found for exposure=%1s, temp=%2K (%3°C), binning=%4, bayer=%5")
                      .arg(lightExposure).arg(lightTemperatureK).arg(temperatureC)
                      .arg(lightBinning).arg(lightBayerPattern), "orange");
        
        // Still process with coordinate conversion even without dark calibration
        if (!writeStellinaMetadataWithCoordinates(lightFrame, imageData)) {
            logMessage("Failed to write metadata with coordinates to original FITS file", "red");
        }
        
        m_skippedCount++;
        return true;
    }
    
    logMessage(QString("Found %1 matching dark frames (including rotated if needed)").arg(matchingDarks.size()), "blue");
    
    // Create master dark with bayer pattern in filename
    QString masterDarkName = QString("master_dark_%1s_%2K_%3_%4.fits")
                                .arg(lightExposure)
                                .arg(lightTemperatureK)
                                .arg(lightBinning)
                                .arg(lightBayerPattern);
    QString masterDarkPath = QDir(m_calibratedDirectory).absoluteFilePath(masterDarkName);
    
    if (!QFile::exists(masterDarkPath)) {
        if (!createMasterDark(matchingDarks, masterDarkPath)) {
            logMessage("Failed to create master dark", "red");
            return false;
        }
        int temperatureC = lightTemperatureK - 273;
        logMessage(QString("Created master dark: %1 (from %2K/%3°C data, %4 pattern)")
                      .arg(masterDarkName).arg(lightTemperatureK).arg(temperatureC).arg(lightBayerPattern), "green");
    } else {
        logMessage(QString("Using existing master dark: %1").arg(masterDarkName), "blue");
    }
    
    // Apply master dark
    QString outputName = QString("dark_calibrated_%1.fits")
                            .arg(QFileInfo(lightFrame).baseName());
    QString outputPath = QDir(m_calibratedDirectory).absoluteFilePath(outputName);
    
    if (applyMasterDark(lightFrame, masterDarkPath, outputPath)) {
        // Write metadata with coordinate conversion
        StellinaImageData calibratedImageData = imageData;
        calibratedImageData.currentFitsPath = outputPath;
        
        if (!writeStellinaMetadataWithCoordinates(outputPath, calibratedImageData)) {
            logMessage("Warning: Failed to write metadata with coordinates to calibrated FITS file", "orange");
        } else {
            updateProcessingStage(outputPath, "DARK_CALIBRATED_WITH_COORDS");
            logMessage(QString("Wrote dark-calibrated file with RA/DEC coordinates: %1").arg(QFileInfo(outputPath).fileName()), "green");
        }
        
        // Update tracking data
        for (int i = 0; i < m_stellinaImageData.size(); ++i) {
            if (m_stellinaImageData[i].originalFitsPath == imageData.originalFitsPath) {
                m_stellinaImageData[i].currentFitsPath = outputPath;
                break;
            }
        }
        
        m_darkCalibratedFiles.append(outputPath);
        m_darkCalibratedCount++;
        
        QString reversedInfo = imageData.isReversedImage ? " (reversed image)" : "";
        logMessage(QString("Dark calibration with coordinate conversion successful: %1%2")
                      .arg(outputName).arg(reversedInfo), "green");
        return true;
    } else {
        logMessage("Dark calibration failed", "red");
        return false;
    }
}

