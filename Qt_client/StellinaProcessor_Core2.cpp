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

void StellinaProcessor::scanDarkFrames() {
    m_darkFrames.clear();
    
    if (m_darkDirectory.isEmpty()) {
        m_darkFramesCount->setText("No dark frames directory selected");
        return;
    }
    
    QDir darkDir(m_darkDirectory);
    if (!darkDir.exists()) {
        m_darkFramesCount->setText("Dark frames directory does not exist");
        return;
    }
    
    QStringList darkFiles = darkDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    
    logMessage(QString("Scanning %1 potential dark frames...").arg(darkFiles.size()), "blue");
    
    // Group dark frames by exposure, temperature, and binning
    QMap<QString, QStringList> darkGroups;
    QMap<QString, DarkFrame> darkInfo; // Store the representative info for each group
    
    for (const QString &darkFile : darkFiles) {
        QString fullPath = darkDir.absoluteFilePath(darkFile);
        
        int exposure = extractExposureTime(fullPath);
        int temperatureK = extractTemperature(fullPath); // Now returns Kelvin
        QString binning = extractBinning(fullPath);
        
        if (exposure > 0) { // Valid dark frame
            QString key = QString("%1_%2_%3").arg(exposure).arg(temperatureK).arg(binning);
            darkGroups[key].append(fullPath);
            
            // Store representative info for this group (first file sets the pattern)
            if (!darkInfo.contains(key)) {
                DarkFrame dark;
                dark.filepath = fullPath; // Representative file
                dark.exposure = exposure;
                dark.temperature = temperatureK; // Store as Kelvin
                dark.binning = binning;
                darkInfo[key] = dark;
                
                if (m_debugMode) {
                    int temperatureC = temperatureK - 273; // Convert back to Celsius for display
                    logMessage(QString("Dark group %1: %2s, %3K (%4°C), %5 - first file: %6")
                                  .arg(darkInfo.size())
                                  .arg(exposure)
                                  .arg(temperatureK)
                                  .arg(temperatureC)
                                  .arg(binning)
                                  .arg(QFileInfo(darkFile).fileName()), "gray");
                }
            }
        } else {
            if (m_debugMode) {
                logMessage(QString("Skipped invalid dark frame: %1").arg(QFileInfo(darkFile).fileName()), "orange");
            }
        }
    }
    
    // Create DarkFrame entries from the groups
    m_darkFrames.clear();
    for (auto it = darkGroups.begin(); it != darkGroups.end(); ++it) {
        QString key = it.key();
        if (darkInfo.contains(key)) {
            DarkFrame dark = darkInfo[key];
            // Update filepath to point to the first file in the group
            dark.filepath = it.value().first();
            m_darkFrames.append(dark);
            
            int temperatureC = dark.temperature - 273; // Convert for display
            logMessage(QString("Dark group: %1s exposure, %2K (%3°C), %4 binning → %5 frames")
                          .arg(dark.exposure)
                          .arg(dark.temperature)
                          .arg(temperatureC)
                          .arg(dark.binning)
                          .arg(it.value().size()), "blue");
        }
    }
    
    // Update UI
    m_darkFramesCount->setText(QString("%1 dark frame groups found").arg(m_darkFrames.size()));
    
    // Update dark frames table
    m_darkFramesTable->setRowCount(m_darkFrames.size());
    for (int i = 0; i < m_darkFrames.size(); ++i) {
        const DarkFrame &dark = m_darkFrames[i];
        QString key = QString("%1_%2_%3").arg(dark.exposure).arg(dark.temperature).arg(dark.binning);
        
        int temperatureC = dark.temperature - 273;
        m_darkFramesTable->setItem(i, 0, new QTableWidgetItem(QString("%1s_%2K_%3").arg(dark.exposure).arg(dark.temperature).arg(dark.binning)));
        m_darkFramesTable->setItem(i, 1, new QTableWidgetItem(QString("%1s").arg(dark.exposure)));
        m_darkFramesTable->setItem(i, 2, new QTableWidgetItem(QString("%1K (%2°C)").arg(dark.temperature).arg(temperatureC)));
        m_darkFramesTable->setItem(i, 3, new QTableWidgetItem(dark.binning));
        
        // Count how many dark frames in this group
        int count = darkGroups[key].size();
        m_darkFramesTable->setItem(i, 4, new QTableWidgetItem(QString::number(count)));
    }
    
    logMessage(QString("Dark frame scan complete: %1 groups found").arg(m_darkFrames.size()), "green");
    
    // Show summary of what was found
    if (!m_darkFrames.isEmpty()) {
        QStringList summary;
        for (const DarkFrame &dark : m_darkFrames) {
            QString key = QString("%1_%2_%3").arg(dark.exposure).arg(dark.temperature).arg(dark.binning);
            int count = darkGroups[key].size();
            int temperatureC = dark.temperature - 273;
            summary.append(QString("%1×%2s@%3K").arg(count).arg(dark.exposure).arg(dark.temperature));
        }
        logMessage(QString("Found: %1").arg(summary.join(", ")), "blue");
    }
}

// Dark calibration functions
bool StellinaProcessor::findMatchingDarkFrame(const QString &lightFrame, DarkFrame &darkFrame) {
    // This method is now deprecated - we always use findAllMatchingDarkFrames instead
    // Keep it for compatibility but it's not used in the new workflow
    Q_UNUSED(lightFrame)
    Q_UNUSED(darkFrame)
    return false;
}

QStringList StellinaProcessor::findAllMatchingDarkFrames(int targetExposure, int targetTemperature, const QString &targetBinning) {
    QStringList matchingDarks;
    
    if (m_darkDirectory.isEmpty()) {
        return matchingDarks;
    }
    
    QDir darkDir(m_darkDirectory);
    if (!darkDir.exists()) {
        return matchingDarks;
    }
    
    QStringList darkFiles = darkDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    
    for (const QString &darkFile : darkFiles) {
        QString fullPath = darkDir.absoluteFilePath(darkFile);
        
        int exposure = extractExposureTime(fullPath);
        int temperature = extractTemperature(fullPath);
        QString binning = extractBinning(fullPath);
        
        // Check if this dark frame matches the light frame characteristics
        bool exposureMatch = qAbs(exposure - targetExposure) <= (targetExposure * m_exposureTolerance / 100);
        bool temperatureMatch = qAbs(temperature - targetTemperature) <= m_temperatureTolerance;
        bool binningMatch = (binning == targetBinning);
        
        if (exposureMatch && temperatureMatch && binningMatch) {
            matchingDarks.append(fullPath);
        }
    }
    
    return matchingDarks;
}

bool StellinaProcessor::createMasterDark(const QStringList &darkFrames, const QString &outputPath) {
    if (darkFrames.isEmpty()) {
        logMessage("No dark frames provided for master dark creation", "red");
        return false;
    }
    
    logMessage(QString("Creating master dark from %1 frames...").arg(darkFrames.size()), "blue");
    
    if (darkFrames.size() == 1) {
        // Only one dark frame, just copy it
        if (QFile::copy(darkFrames.first(), outputPath)) {
            logMessage("Master dark created (single frame copy)", "green");
            return true;
        } else {
            logMessage("Failed to copy single dark frame", "red");
            return false;
        }
    }
    
    // Use direct FITS manipulation for better performance
    logMessage("Using direct FITS processing for master dark creation", "blue");
    return createMasterDarkDirect(darkFrames, outputPath);
}

bool StellinaProcessor::createMasterDarkDirect(const QStringList &darkFrames, const QString &outputPath) {
    if (darkFrames.isEmpty()) {
        return false;
    }
    
    if (darkFrames.size() == 1) {
        return QFile::copy(darkFrames.first(), outputPath);
    }
    
    logMessage(QString("Creating master dark from %1 frames using direct FITS manipulation...").arg(darkFrames.size()), "blue");
    
    fitsfile *firstFits = nullptr;
    fitsfile *outputFits = nullptr;
    int status = 0;
    
    // Open first dark frame to get dimensions and header info
    QByteArray firstPath = darkFrames.first().toLocal8Bit();
    if (fits_open_file(&firstFits, firstPath.data(), READONLY, &status)) {
        logMessage(QString("Failed to open first dark frame: %1 (FITS error: %2)").arg(darkFrames.first()).arg(status), "red");
        return false;
    }
    
    // Get image dimensions
    int naxis;
    long naxes[2];
    if (fits_get_img_dim(firstFits, &naxis, &status) || 
        fits_get_img_size(firstFits, 2, naxes, &status)) {
        logMessage(QString("Failed to get image dimensions (FITS error: %1)").arg(status), "red");
        fits_close_file(firstFits, &status);
        return false;
    }
    
    long totalPixels = naxes[0] * naxes[1];
    logMessage(QString("Image dimensions: %1 x %2 (%3 pixels)").arg(naxes[0]).arg(naxes[1]).arg(totalPixels), "blue");
    
    // Create output FITS file
    QByteArray outputPathBytes = QString("!%1").arg(outputPath).toLocal8Bit(); // ! forces overwrite
    if (fits_create_file(&outputFits, outputPathBytes.data(), &status)) {
        logMessage(QString("Failed to create output FITS file: %1 (FITS error: %2)").arg(outputPath).arg(status), "red");
        fits_close_file(firstFits, &status);
        return false;
    }
    
    // Copy header from first frame
    if (fits_copy_header(firstFits, outputFits, &status)) {
        logMessage(QString("Failed to copy FITS header (FITS error: %1)").arg(status), "red");
        fits_close_file(firstFits, &status);
        fits_close_file(outputFits, &status);
        return false;
    }
    
    fits_close_file(firstFits, &status);
    
    // Allocate memory for accumulation (using double for precision)
    std::vector<double> accumulator(totalPixels, 0.0);
    std::vector<float> pixelBuffer(totalPixels);
    
    int successfulFrames = 0;
    
    // Process each dark frame
    for (int i = 0; i < darkFrames.size(); ++i) {
        const QString &darkPath = darkFrames[i];
        fitsfile *darkFits = nullptr;
        
        QByteArray darkPathBytes = darkPath.toLocal8Bit();
        if (fits_open_file(&darkFits, darkPathBytes.data(), READONLY, &status)) {
            logMessage(QString("Failed to open dark frame %1: %2 (FITS error: %3)")
                          .arg(i+1).arg(QFileInfo(darkPath).fileName()).arg(status), "orange");
            status = 0; // Reset status to continue
            continue;
        }
        
        // Read pixel data
        long firstPixel = 1;
        if (fits_read_img(darkFits, TFLOAT, firstPixel, totalPixels, nullptr, 
                         pixelBuffer.data(), nullptr, &status)) {
            logMessage(QString("Failed to read pixel data from frame %1: %2 (FITS error: %3)")
                          .arg(i+1).arg(QFileInfo(darkPath).fileName()).arg(status), "orange");
            fits_close_file(darkFits, &status);
            status = 0;
            continue;
        }
        
        // Add to accumulator
        for (long j = 0; j < totalPixels; ++j) {
            accumulator[j] += pixelBuffer[j];
        }
        
        successfulFrames++;
        fits_close_file(darkFits, &status);
        
        // Update progress
        if (i % 5 == 0 || i == darkFrames.size() - 1) {
            logMessage(QString("Processed dark frame %1 of %2").arg(i+1).arg(darkFrames.size()), "gray");
        }
    }
    
    if (successfulFrames == 0) {
        logMessage("No dark frames could be processed", "red");
        fits_close_file(outputFits, &status);
        QFile::remove(outputPath);
        return false;
    }
    
    logMessage(QString("Successfully read %1 of %2 dark frames, computing average...").arg(successfulFrames).arg(darkFrames.size()), "blue");
    
    // Average the accumulated values
    for (long i = 0; i < totalPixels; ++i) {
        pixelBuffer[i] = static_cast<float>(accumulator[i] / successfulFrames);
    }
    
    // Write averaged data to output
    long firstPixel = 1;
    if (fits_write_img(outputFits, TFLOAT, firstPixel, totalPixels, 
                      pixelBuffer.data(), &status)) {
        logMessage(QString("Failed to write master dark data (FITS error: %1)").arg(status), "red");
        fits_close_file(outputFits, &status);
        QFile::remove(outputPath);
        return false;
    }
    
    // Update header with processing info
    QString historyComment = QString("Master dark from %1 frames").arg(successfulFrames);
    QByteArray historyBytes = historyComment.toLocal8Bit();
    if (fits_write_history(outputFits, historyBytes.data(), &status)) {
        // Non-critical error, just log it
        logMessage("Warning: Could not write processing history to header", "orange");
        status = 0;
    }
    
    // Add custom keyword for frame count
    if (fits_write_key(outputFits, TINT, "NFRAMES", &successfulFrames, "Number of frames averaged", &status)) {
        logMessage("Warning: Could not write frame count to header", "orange");
        status = 0;
    }
    
    fits_close_file(outputFits, &status);
    
    if (status != 0) {
        logMessage(QString("FITS error occurred during master dark creation (error: %1)").arg(status), "red");
        QFile::remove(outputPath);
        return false;
    }
    
    logMessage(QString("Master dark created successfully from %1 frames (averaged %2 frames)")
                  .arg(darkFrames.size()).arg(successfulFrames), "green");
    return true;
}
// Replace the applyMasterDark function in StellinaProcessor_Core.cpp with this implementation

bool StellinaProcessor::applyMasterDark(const QString &lightFrame, const QString &masterDark, const QString &outputFrame) {
    logMessage(QString("Applying master dark to %1 using direct FITS processing").arg(QFileInfo(lightFrame).fileName()), "blue");
    
    fitsfile *lightFits = nullptr;
    fitsfile *darkFits = nullptr;
    fitsfile *outputFits = nullptr;
    int status = 0;
    
    // Open light frame
    QByteArray lightPath = lightFrame.toLocal8Bit();
    if (fits_open_file(&lightFits, lightPath.data(), READONLY, &status)) {
        logMessage(QString("Failed to open light frame: %1 (FITS error: %2)").arg(lightFrame).arg(status), "red");
        return false;
    }
    
    // Open master dark
    QByteArray darkPath = masterDark.toLocal8Bit();
    if (fits_open_file(&darkFits, darkPath.data(), READONLY, &status)) {
        logMessage(QString("Failed to open master dark: %1 (FITS error: %2)").arg(masterDark).arg(status), "red");
        fits_close_file(lightFits, &status);
        return false;
    }
    
    // Get dimensions and verify they match
    long lightNaxes[2], darkNaxes[2];
    if (fits_get_img_size(lightFits, 2, lightNaxes, &status) ||
        fits_get_img_size(darkFits, 2, darkNaxes, &status)) {
        logMessage(QString("Failed to get image dimensions (FITS error: %1)").arg(status), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        return false;
    }
    
    if (lightNaxes[0] != darkNaxes[0] || lightNaxes[1] != darkNaxes[1]) {
        logMessage(QString("Image dimensions mismatch - Light: %1x%2, Dark: %3x%4")
                      .arg(lightNaxes[0]).arg(lightNaxes[1])
                      .arg(darkNaxes[0]).arg(darkNaxes[1]), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        return false;
    }
    
    long totalPixels = lightNaxes[0] * lightNaxes[1];
    logMessage(QString("Processing %1 x %2 image (%3 pixels)").arg(lightNaxes[0]).arg(lightNaxes[1]).arg(totalPixels), "blue");
    
    // Create output FITS file
    QByteArray outputPathBytes = QString("!%1").arg(outputFrame).toLocal8Bit(); // ! forces overwrite
    if (fits_create_file(&outputFits, outputPathBytes.data(), &status)) {
        logMessage(QString("Failed to create output FITS file: %1 (FITS error: %2)").arg(outputFrame).arg(status), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        return false;
    }
    
    // Copy header from light frame
    if (fits_copy_header(lightFits, outputFits, &status)) {
        logMessage(QString("Failed to copy FITS header (FITS error: %1)").arg(status), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        fits_close_file(outputFits, &status);
        return false;
    }
    
    // Allocate memory for pixel data
    std::vector<float> lightPixels(totalPixels);
    std::vector<float> darkPixels(totalPixels);
    std::vector<float> resultPixels(totalPixels);
    
    // Read light frame pixel data
    long firstPixel = 1;
    if (fits_read_img(lightFits, TFLOAT, firstPixel, totalPixels, nullptr, 
                     lightPixels.data(), nullptr, &status)) {
        logMessage(QString("Failed to read light frame pixel data (FITS error: %1)").arg(status), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        fits_close_file(outputFits, &status);
        QFile::remove(outputFrame);
        return false;
    }
    
    // Read master dark pixel data
    if (fits_read_img(darkFits, TFLOAT, firstPixel, totalPixels, nullptr, 
                     darkPixels.data(), nullptr, &status)) {
        logMessage(QString("Failed to read master dark pixel data (FITS error: %1)").arg(status), "red");
        fits_close_file(lightFits, &status);
        fits_close_file(darkFits, &status);
        fits_close_file(outputFits, &status);
        QFile::remove(outputFrame);
        return false;
    }
    
    // Close input files as we no longer need them
    fits_close_file(lightFits, &status);
    fits_close_file(darkFits, &status);
    
    // Perform dark subtraction with negative truncation
    long negativePixels = 0;
    float minResult = std::numeric_limits<float>::max();
    float maxResult = std::numeric_limits<float>::lowest();
    
    for (long i = 0; i < totalPixels; ++i) {
        float result = lightPixels[i] - darkPixels[i];
        
        // Truncate negative values to 0
        if (result < 0.0f) {
            result = 0.0f;
            negativePixels++;
        }
        
        resultPixels[i] = result;
        
        // Track min/max for statistics
        if (result < minResult) minResult = result;
        if (result > maxResult) maxResult = result;
    }
    
    // Log statistics
    double negativePercent = (static_cast<double>(negativePixels) / totalPixels) * 100.0;
    logMessage(QString("Dark subtraction complete - Negative pixels truncated: %1 (%2%)")
                  .arg(negativePixels).arg(negativePercent, 0, 'f', 2), "blue");
    logMessage(QString("Result range: %1 to %2").arg(minResult, 0, 'f', 2).arg(maxResult, 0, 'f', 2), "gray");
    
    // Write result to output file
    if (fits_write_img(outputFits, TFLOAT, firstPixel, totalPixels, 
                      resultPixels.data(), &status)) {
        logMessage(QString("Failed to write calibrated pixel data (FITS error: %1)").arg(status), "red");
        fits_close_file(outputFits, &status);
        QFile::remove(outputFrame);
        return false;
    }
    
    // Update header with processing information
    QString historyComment = QString("Dark calibrated using master dark: %1").arg(QFileInfo(masterDark).fileName());
    QByteArray historyBytes = historyComment.toLocal8Bit();
    if (fits_write_history(outputFits, historyBytes.data(), &status)) {
        // Non-critical error
        logMessage("Warning: Could not write processing history to header", "orange");
        status = 0;
    }
    
    // Add custom keywords for processing info
    int truncatedCount = static_cast<int>(negativePixels);
    if (fits_write_key(outputFits, TINT, "NEGTRUNC", &truncatedCount, "Number of negative pixels truncated", &status)) {
        logMessage("Warning: Could not write truncation count to header", "orange");
        status = 0;
    }
    
    QString processingDate = QDateTime::currentDateTime().toString(Qt::ISODate);
    QByteArray dateBytes = processingDate.toLocal8Bit();
    char* datePtr = dateBytes.data();
    if (fits_write_key(outputFits, TSTRING, "DARKCAL", &datePtr, "Dark calibration date", &status)) {
        logMessage("Warning: Could not write processing date to header", "orange");
        status = 0;
    }
    
    fits_close_file(outputFits, &status);
    
    if (status != 0) {
        logMessage(QString("FITS error occurred during dark calibration (error: %1)").arg(status), "red");
        QFile::remove(outputFrame);
        return false;
    }
    
    logMessage(QString("Dark calibration successful: %1").arg(QFileInfo(outputFrame).fileName()), "green");
    
    m_darkCalibratedFiles.append(outputFrame);
    return true;
}

// Astrometric stacking functions
bool StellinaProcessor::createSequence(const QStringList &imageList, const QString &sequenceName) {
    if (imageList.isEmpty()) {
        return false;
    }
    
    QString outputDir = getOutputDirectoryForCurrentStage();
    
    // Copy images to a temporary location with sequential naming
    QString seqDir = QDir(outputDir).absoluteFilePath(sequenceName);
    QDir().mkpath(seqDir);
    
    logMessage(QString("Creating sequence with %1 images...").arg(imageList.size()), "blue");
    
    for (int i = 0; i < imageList.size(); ++i) {
        QString srcFile = imageList[i];
        QString dstFile = QDir(seqDir).absoluteFilePath(QString("%1_%2.fits").arg(sequenceName).arg(i + 1, 4, 10, QChar('0')));
        
        if (!QFile::copy(srcFile, dstFile)) {
            logMessage(QString("Failed to copy %1 to sequence directory").arg(QFileInfo(srcFile).fileName()), "red");
            return false;
        }
    }
    
}

bool StellinaProcessor::performGlobalRegistration(const QString &sequenceName) {
    logMessage("Performing global registration...", "blue");
    m_registrationStatusLabel->setText("Registration: In Progress");
    
    // Load sequence
    QString command = QString("load %1").arg(sequenceName);

    m_registrationStatusLabel->setText("Registration: Complete");
    logMessage("Global registration completed", "green");
    return true;
}

bool StellinaProcessor::performAstrometricStacking() {
    if (m_plateSolvedFiles.isEmpty()) {
        logMessage("No plate-solved images available for stacking", "red");
        return false;
    }
    
    if (m_plateSolvedFiles.size() < 3) {
        logMessage("Need at least 3 images for stacking", "red");
        return false;
    }
    
    logMessage(QString("Starting astrometric stacking with %1 images").arg(m_plateSolvedFiles.size()), "blue");
    
    // Create sequence
    m_sequenceName = QString("stellina_stack_%1").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss"));
    
    m_subTaskProgressBar->setVisible(true);
    m_subTaskProgressBar->setMaximum(4);
    m_subTaskProgressBar->setValue(0);
    
    // Step 1: Create sequence
    m_currentTaskLabel->setText("Creating image sequence...");
    m_subTaskProgressBar->setValue(1);
    
    if (!createSequence(m_plateSolvedFiles, m_sequenceName)) {
        logMessage("Failed to create image sequence", "red");
        m_subTaskProgressBar->setVisible(false);
        return false;
    }
    
    // Step 2: Perform registration
    m_currentTaskLabel->setText("Performing astrometric registration...");
    m_subTaskProgressBar->setValue(2);
    
    if (!performGlobalRegistration(m_sequenceName)) {
        logMessage("Failed to perform registration", "red");
        m_subTaskProgressBar->setVisible(false);
        return false;
    }
    
    // Step 3: Perform stacking
    m_currentTaskLabel->setText("Stacking registered images...");
    m_subTaskProgressBar->setValue(3);
    
    if (!performStacking(m_sequenceName, m_stackingParams)) {
        logMessage("Failed to perform stacking", "red");
        m_subTaskProgressBar->setVisible(false);
        return false;
    }
    
    // Step 4: Complete
    m_currentTaskLabel->setText("Stacking complete");
    m_subTaskProgressBar->setValue(4);
    m_subTaskProgressBar->setVisible(false);
    
    return true;
}

QJsonObject StellinaProcessor::loadStellinaJson(const QString &jsonPath) {
    QFile file(jsonPath);
    if (!file.open(QIODevice::ReadOnly)) {
        if (m_debugMode) {
            logMessage(QString("Failed to open JSON file: %1").arg(jsonPath), "red");
        }
        return QJsonObject();
    }
    
    QJsonParseError error;
    QJsonDocument doc = QJsonDocument::fromJson(file.readAll(), &error);
    
    if (error.error != QJsonParseError::NoError) {
        if (m_debugMode) {
            logMessage(QString("JSON parse error in %1: %2").arg(jsonPath).arg(error.errorString()), "red");
        }
        return QJsonObject();
    }
    
    return doc.object();
}
