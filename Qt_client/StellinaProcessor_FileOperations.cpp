// StellinaProcessor_FileOperations.cpp - File management and operations
// Extracted from StellinaProcessor_Core.cpp

#include "StellinaProcessor.h"
#include <QDir>
#include <QFileInfo>
#include <QRegularExpression>
#include <QStandardPaths>

// ============================================================================
// File Discovery and Validation
// ============================================================================

bool StellinaProcessor::findStellinaImages() {
    m_imagesToProcess.clear();
    
    if (m_sourceDirectory.isEmpty()) {
        logMessage("No source directory specified", "red");
        return false;
    }
    
    QDir sourceDir(m_sourceDirectory);
    if (!sourceDir.exists()) {
        logMessage(QString("Source directory does not exist: %1").arg(m_sourceDirectory), "red");
        return false;
    }
    
    // Look for FITS files
    QStringList fitsFiles = sourceDir.entryList(
        QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT",
        QDir::Files, QDir::Name);
    
    if (fitsFiles.isEmpty()) {
        logMessage("No FITS files found in source directory", "red");
        return false;
    }
    
    // Validate and filter files
    int validFiles = 0;
    for (const QString &fitsFile : fitsFiles) {
        QString fullPath = sourceDir.absoluteFilePath(fitsFile);
        
        if (validateImageFile(fullPath)) {
            m_imagesToProcess.append(fullPath);
            validFiles++;
        } else {
            logMessage(QString("Skipping invalid file: %1").arg(fitsFile), "orange");
        }
    }
    
    if (validFiles == 0) {
        logMessage("No valid Stellina images found", "red");
        return false;
    }
    
    logMessage(QString("Found %1 valid FITS files").arg(validFiles), "blue");
    
    // Sort files by name for consistent processing order
    m_imagesToProcess.sort();
    
    return true;
}

bool StellinaProcessor::validateImageFile(const QString &filePath) {
    // Check file exists and is readable
    QFileInfo fileInfo(filePath);
    if (!fileInfo.exists() || !fileInfo.isReadable()) {
        return false;
    }
    
    // Check file size (too small likely corrupt)
    if (fileInfo.size() < 1024) {
        if (m_debugMode) {
            logMessage(QString("File too small: %1").arg(filePath), "gray");
        }
        return false;
    }
    
    // Validate FITS format
    if (!validateFITSFile(filePath)) {
        if (m_debugMode) {
            logMessage(QString("Invalid FITS format: %1").arg(filePath), "gray");
        }
        return false;
    }
    
    // Check for Stellina naming pattern (optional - some files may not follow pattern)
    QString filename = fileInfo.baseName();
    QRegularExpression stellinaPattern(R"(img-\d{4}[r]?)");
    
    if (!stellinaPattern.match(filename).hasMatch()) {
        // Not a strict Stellina pattern, but still allow if it's a valid FITS
        if (m_debugMode) {
            logMessage(QString("Non-standard filename pattern: %1").arg(filename), "gray");
        }
    }
    
    return true;
}

QStringList StellinaProcessor::findMatchingFiles(const QString &directory, const QStringList &patterns) {
    QStringList matchingFiles;
    
    QDir dir(directory);
    if (!dir.exists()) {
        return matchingFiles;
    }
    
    QStringList files = dir.entryList(patterns, QDir::Files, QDir::Name);
    
    for (const QString &file : files) {
        QString fullPath = dir.absoluteFilePath(file);
        matchingFiles.append(fullPath);
    }
    
    return matchingFiles;
}

QStringList StellinaProcessor::findFilesWithMetadata(const QString &directory) {
    QStringList filesWithMetadata;
    
    QStringList fitsFiles = findMatchingFiles(directory, QStringList() << "*.fits");
    
    for (const QString &fitsPath : fitsFiles) {
        // Check if file has Stellina metadata or corresponding JSON
        if (hasStellinaMetadata(fitsPath)) {
            filesWithMetadata.append(fitsPath);
        } else {
            // Check for corresponding JSON file
            QString jsonPath = fitsPath;
            jsonPath.replace(".fits", ".json");
            if (QFile::exists(jsonPath)) {
                filesWithMetadata.append(fitsPath);
            }
        }
    }
    
    return filesWithMetadata;
}

// ============================================================================
// File Path Utilities
// ============================================================================

QString StellinaProcessor::generateOutputPath(const QString &inputPath, const QString &prefix, const QString &directory) {
    QFileInfo inputInfo(inputPath);
    QString baseName = inputInfo.baseName();
    
    // Remove existing prefixes to avoid duplication
    const QStringList existingPrefixes = {"dark_calibrated_", "plate_solved_", "registered_"};
    for (const QString &existingPrefix : existingPrefixes) {
        if (baseName.startsWith(existingPrefix)) {
            baseName = baseName.mid(existingPrefix.length());
            break;
        }
    }
    
    QString outputName = QString("%1%2.fits").arg(prefix).arg(baseName);
    return QDir(directory).absoluteFilePath(outputName);
}

QString StellinaProcessor::generateUniqueOutputPath(const QString &inputPath, const QString &prefix, const QString &directory) {
    QString basePath = generateOutputPath(inputPath, prefix, directory);
    
    if (!QFile::exists(basePath)) {
        return basePath;
    }
    
    // Generate unique name with counter
    QFileInfo baseInfo(basePath);
    QString baseName = baseInfo.baseName();
    QString extension = baseInfo.suffix();
    QString dirPath = baseInfo.absolutePath();
    
    int counter = 1;
    QString uniquePath;
    
    do {
        uniquePath = QDir(dirPath).absoluteFilePath(
            QString("%1_%2.%3").arg(baseName).arg(counter).arg(extension));
        counter++;
    } while (QFile::exists(uniquePath) && counter < 1000);
    
    return uniquePath;
}

bool StellinaProcessor::ensureDirectoryExists(const QString &directory) {
    if (directory.isEmpty()) {
        return false;
    }
    
    QDir dir(directory);
    if (dir.exists()) {
        return true;
    }
    
    if (dir.mkpath(".")) {
        logMessage(QString("Created directory: %1").arg(directory), "blue");
        return true;
    } else {
        logMessage(QString("Failed to create directory: %1").arg(directory), "red");
        return false;
    }
}

bool StellinaProcessor::createOutputDirectories() {
    bool success = true;
    
    if (!m_calibratedDirectory.isEmpty()) {
        success &= ensureDirectoryExists(m_calibratedDirectory);
    }
    
    if (!m_plateSolvedDirectory.isEmpty()) {
        success &= ensureDirectoryExists(m_plateSolvedDirectory);
    }
    
    if (!m_stackedDirectory.isEmpty()) {
        success &= ensureDirectoryExists(m_stackedDirectory);
    }
    
    return success;
}

// ============================================================================
// File Backup and Safety
// ============================================================================

bool StellinaProcessor::createBackup(const QString &filePath) {
    if (!QFile::exists(filePath)) {
        return false;
    }
    
    QString backupPath = filePath + ".backup";
    
    // Remove existing backup
    if (QFile::exists(backupPath)) {
        QFile::remove(backupPath);
    }
    
    return QFile::copy(filePath, backupPath);
}

bool StellinaProcessor::restoreBackup(const QString &filePath) {
    QString backupPath = filePath + ".backup";
    
    if (!QFile::exists(backupPath)) {
        return false;
    }
    
    // Remove current file
    if (QFile::exists(filePath)) {
        QFile::remove(filePath);
    }
    
    // Restore from backup
    return QFile::copy(backupPath, filePath);
}

bool StellinaProcessor::cleanupBackups(const QString &directory) {
    QDir dir(directory);
    if (!dir.exists()) {
        return false;
    }
    
    QStringList backupFiles = dir.entryList(QStringList() << "*.backup", QDir::Files);
    
    int removedCount = 0;
    for (const QString &backupFile : backupFiles) {
        QString backupPath = dir.absoluteFilePath(backupFile);
        if (QFile::remove(backupPath)) {
            removedCount++;
        }
    }
    
    if (removedCount > 0) {
        logMessage(QString("Cleaned up %1 backup files").arg(removedCount), "blue");
    }
    
    return true;
}

// ============================================================================
// Temporary File Management
// ============================================================================

bool StellinaProcessor::cleanupTemporaryFiles() {
    bool success = true;
    
    // Clean up solve-field temporary files
    for (const QString &directory : {m_calibratedDirectory, m_plateSolvedDirectory}) {
        if (!directory.isEmpty()) {
            success &= cleanupSolveFieldTemporaryFiles(directory);
        }
    }
    
    // Clean up processing temporary files
    QString tempDir = QStandardPaths::writableLocation(QStandardPaths::TempLocation);
    QDir temp(tempDir);
    
    QStringList tempPatterns = {
        "stellina_temp_*",
        "solve_temp_*",
        "stack_temp_*"
    };
    
    for (const QString &pattern : tempPatterns) {
        QStringList tempFiles = temp.entryList(QStringList() << pattern, QDir::Files);
        for (const QString &tempFile : tempFiles) {
            QString tempPath = temp.absoluteFilePath(tempFile);
            if (QFile::remove(tempPath)) {
                if (m_debugMode) {
                    logMessage(QString("Removed temp file: %1").arg(tempFile), "gray");
                }
            }
        }
    }
    
    return success;
}

bool StellinaProcessor::cleanupSolveFieldTemporaryFiles(const QString &directory) {
    QDir dir(directory);
    if (!dir.exists()) {
        return false;
    }
    
    // solve-field creates these temporary files
    QStringList tempExtensions = {
        "*.axy", "*.corr", "*.match", "*.rdls", "*.solved", "*.wcs",
        "*.xyls", "*.new", "*-indx.png", "*-ngc.png"
    };
    
    int removedCount = 0;
    for (const QString &pattern : tempExtensions) {
        QStringList tempFiles = dir.entryList(QStringList() << pattern, QDir::Files);
        for (const QString &tempFile : tempFiles) {
            QString tempPath = dir.absoluteFilePath(tempFile);
            if (QFile::remove(tempPath)) {
                removedCount++;
            }
        }
    }
    
    if (removedCount > 0 && m_debugMode) {
        logMessage(QString("Cleaned up %1 solve-field temp files in %2")
                  .arg(removedCount).arg(directory), "gray");
    }
    
    return true;
}

QString StellinaProcessor::createTemporaryFile(const QString &prefix, const QString &suffix) {
    QString tempDir = QStandardPaths::writableLocation(QStandardPaths::TempLocation);
    QString timestamp = QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss_zzz");
    
    QString tempPath = QDir(tempDir).absoluteFilePath(
        QString("%1_%2%3").arg(prefix).arg(timestamp).arg(suffix));
    
    return tempPath;
}

// ============================================================================
// File Metadata and Tracking
// ============================================================================

StellinaImageData* StellinaProcessor::findImageDataByPath(const QString &path) {
    // Search by current path first
    for (auto &data : m_stellinaImageData) {
        if (data.currentFitsPath == path) {
            return &data;
        }
    }
    
    // Search by original path
    for (auto &data : m_stellinaImageData) {
        if (data.originalFitsPath == path) {
            return &data;
        }
    }
    
    return nullptr;
}

StellinaImageData* StellinaProcessor::findImageDataByOriginalPath(const QString &originalPath) {
    for (auto &data : m_stellinaImageData) {
        if (data.originalFitsPath == originalPath) {
            return &data;
        }
    }
    
    return nullptr;
}

bool StellinaProcessor::updateImageDataPath(const QString &oldPath, const QString &newPath) {
    StellinaImageData* data = findImageDataByPath(oldPath);
    if (data) {
        data->currentFitsPath = newPath;
        return true;
    }
    
    return false;
}

QStringList StellinaProcessor::getProcessedFilePaths(const QString &stage) const {
    QStringList paths;
    
    for (const auto &data : m_stellinaImageData) {
        if (data.processingStage == stage) {
            paths.append(data.currentFitsPath);
        }
    }
    
    return paths;
}

// ============================================================================
// File Statistics and Analysis
// ============================================================================

void StellinaProcessor::analyzeFileStatistics() {
    logMessage("=== FILE STATISTICS ===", "blue");
    
    // Analyze source files
    if (!m_sourceDirectory.isEmpty()) {
        QDir sourceDir(m_sourceDirectory);
        QStringList sourceFiles = sourceDir.entryList(QStringList() << "*.fits", QDir::Files);
        
        qint64 totalSize = 0;
        for (const QString &file : sourceFiles) {
            QFileInfo info(sourceDir.absoluteFilePath(file));
            totalSize += info.size();
        }
        
        logMessage(QString("Source files: %1 files, %2 MB total")
                  .arg(sourceFiles.size())
                  .arg(totalSize / (1024 * 1024)), "blue");
    }
    
    // Analyze processing stages
    const QStringList stages = {"DARK_CALIBRATED", "PLATE_SOLVED", "REGISTERED", "STACKED"};
    for (const QString &stage : stages) {
        QStringList stageFiles = getProcessedFilePaths(stage);
        logMessage(QString("%1: %2 files").arg(stage).arg(stageFiles.size()), "blue");
    }
    
    // Analyze disk usage
    qint64 totalProcessedSize = 0;
    for (const QString &dir : {m_calibratedDirectory, m_plateSolvedDirectory, m_stackedDirectory}) {
        if (!dir.isEmpty() && QDir(dir).exists()) {
            QDir procDir(dir);
            QStringList procFiles = procDir.entryList(QStringList() << "*.fits", QDir::Files);
            
            for (const QString &file : procFiles) {
                QFileInfo info(procDir.absoluteFilePath(file));
                totalProcessedSize += info.size();
            }
        }
    }
    
    logMessage(QString("Processed files total: %1 MB").arg(totalProcessedSize / (1024 * 1024)), "blue");
}

bool StellinaProcessor::checkDiskSpace(const QString &directory, qint64 requiredMB) {
    QFileInfo dirInfo(directory);
    
    if (!dirInfo.exists()) {
        return false;
    }
    
    // Note: QStorageInfo would be more accurate but requires Qt 5.4+
    // For now, we'll do a basic check
    QString testFile = QDir(directory).absoluteFilePath("disk_test_temp");
    QFile test(testFile);
    
    if (test.open(QIODevice::WriteOnly)) {
        test.close();
        QFile::remove(testFile);
        return true; // Basic write test passed
    }
    
    return false;
}

void StellinaProcessor::generateFileReport() {
    QStringList report;
    
    report << "=== STELLINA PROCESSOR FILE REPORT ===";
    report << QString("Generated: %1").arg(QDateTime::currentDateTime().toString());
    report << "";
    
    report << "Directories:";
    report << QString("  Source: %1").arg(m_sourceDirectory);
    report << QString("  Dark frames: %1").arg(m_darkDirectory);
    report << QString("  Calibrated: %1").arg(m_calibratedDirectory);
    report << QString("  Plate solved: %1").arg(m_plateSolvedDirectory);
    report << QString("  Stacked: %1").arg(m_stackedDirectory);
    report << "";
    
    report << "Processing Statistics:";
    report << QString("  Images processed: %1").arg(m_processedCount);
    report << QString("  Errors encountered: %1").arg(m_errorCount);
    report << QString("  Dark calibrated: %1").arg(m_darkCalibratedCount);
    report << QString("  Files tracked: %1").arg(m_stellinaImageData.size());
    report << "";
    
    // Write report to file
    QString reportPath = QDir(m_stackedDirectory).absoluteFilePath("processing_report.txt");
    QFile reportFile(reportPath);
    
    if (reportFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream stream(&reportFile);
        for (const QString &line : report) {
            stream << line << "\n";
            logMessage(line, "blue");
        }
        
        logMessage(QString("File report saved: %1").arg(reportPath), "green");
    }
}