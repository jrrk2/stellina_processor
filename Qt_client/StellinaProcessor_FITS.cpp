// StellinaProcessor_FITS.cpp - FITS file operations extracted from Core
// This demonstrates how to extract FITS-related functions into their own module

#include "StellinaProcessor.h"
#include <QFileInfo>
#include <QDir>
#include <QDebug>
#include <opencv2/opencv.hpp>
#include <fitsio.h>

// ============================================================================
// FITS Metadata Operations
// ============================================================================

bool StellinaProcessor::readStellinaMetadataFromFits(const QString &fitsPath, StellinaImageData &imageData) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    // Open FITS file
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        if (m_debugMode) {
            logMessage(QString("Could not open FITS file for metadata reading: %1").arg(fitsPath), "orange");
        }
        return false;
    }
    
    // Read Stellina-specific keywords
    char keyword_value[FLEN_VALUE];
    char comment[FLEN_COMMENT];
    
    // Read altitude
    if (fits_read_key(fptr, TDOUBLE, "ST_ALT", &imageData.altitude, comment, &status) == 0) {
        imageData.hasValidCoordinates = true;
    } else {
        status = 0; // Reset status for next read
    }
    
    // Read azimuth
    if (fits_read_key(fptr, TDOUBLE, "ST_AZ", &imageData.azimuth, comment, &status) == 0) {
        imageData.hasValidCoordinates = true;
    } else {
        status = 0;
    }
    
    // Read calculated RA/Dec if present
    if (fits_read_key(fptr, TDOUBLE, "ST_RA", &imageData.calculatedRA, comment, &status) == 0) {
        imageData.hasCalculatedCoords = true;
    } else {
        status = 0;
    }
    
    if (fits_read_key(fptr, TDOUBLE, "ST_DEC", &imageData.calculatedDec, comment, &status) == 0) {
        imageData.hasCalculatedCoords = true;
    } else {
        status = 0;
    }
    
    // Read observation time
    if (fits_read_key(fptr, TSTRING, "DATE-OBS", keyword_value, comment, &status) == 0) {
        imageData.dateObs = QString(keyword_value);
    } else {
        status = 0;
    }
    
    // Read exposure time (use exposureSeconds field)
    double exptime;
    if (fits_read_key(fptr, TDOUBLE, "EXPTIME", &exptime, comment, &status) == 0) {
        imageData.exposureSeconds = static_cast<int>(exptime);
    } else {
        status = 0;
    }
    
    // Read Bayer pattern
    if (fits_read_key(fptr, TSTRING, "BAYERPAT", keyword_value, comment, &status) == 0) {
        imageData.bayerPattern = QString(keyword_value);
    } else {
        status = 0;
    }
    
    fits_close_file(fptr, &status);
    
    if (m_debugMode) {
        logMessage(QString("Read FITS metadata: Alt=%1°, Az=%2°")
                  .arg(imageData.altitude, 0, 'f', 2)
                  .arg(imageData.azimuth, 0, 'f', 2), "gray");
    }
    
    return imageData.hasValidCoordinates;
}

bool StellinaProcessor::writeStellinaMetadataWithCoordinates(const QString &fitsPath, const StellinaImageData &imageData) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    // Open FITS file for writing
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READWRITE, &status)) {
        logMessage(QString("Could not open FITS file for metadata writing: %1").arg(fitsPath), "red");
        return false;
    }
    
    // Clean existing Stellina keywords first
    cleanExistingStellinaKeywordsFITS(fptr);
    
    // Write Stellina coordinates
    if (fits_write_key(fptr, TDOUBLE, "ST_ALT", const_cast<double*>(&imageData.altitude), 
                      "Stellina altitude (degrees)", &status)) {
        logMessage("Failed to write ST_ALT keyword", "orange");
        status = 0;
    }
    
    if (fits_write_key(fptr, TDOUBLE, "ST_AZ", const_cast<double*>(&imageData.azimuth), 
                      "Stellina azimuth (degrees)", &status)) {
        logMessage("Failed to write ST_AZ keyword", "orange");
        status = 0;
    }
    
    // Write calculated RA/Dec if available
    if (imageData.hasCalculatedCoords) {
        if (fits_write_key(fptr, TDOUBLE, "ST_RA", const_cast<double*>(&imageData.calculatedRA), 
                          "Calculated RA (degrees)", &status)) {
            logMessage("Failed to write ST_RA keyword", "orange");
            status = 0;
        }
        
        if (fits_write_key(fptr, TDOUBLE, "ST_DEC", const_cast<double*>(&imageData.calculatedDec), 
                          "Calculated Dec (degrees)", &status)) {
            logMessage("Failed to write ST_DEC keyword", "orange");
            status = 0;
        }
    }
    
    // Write exposure time if not already present
    if (imageData.exposureSeconds > 0) {
        double exptime = static_cast<double>(imageData.exposureSeconds);
        if (fits_write_key(fptr, TDOUBLE, "EXPTIME", &exptime, 
                          "Exposure time (seconds)", &status)) {
            logMessage("Failed to write EXPTIME keyword", "orange");
            status = 0;
        }
    }
    
    // Write Bayer pattern
    if (!imageData.bayerPattern.isEmpty()) {
        QByteArray bayerBytes = imageData.bayerPattern.toLocal8Bit();
        if (fits_write_key(fptr, TSTRING, "BAYERPAT", bayerBytes.data(), 
                          "Bayer pattern", &status)) {
            logMessage("Failed to write BAYERPAT keyword", "orange");
            status = 0;
        }
    }
    
    // Write timestamp
    QString timestamp = QDateTime::currentDateTime().toString(Qt::ISODate);
    QByteArray timestampBytes = timestamp.toLocal8Bit();
    if (fits_write_key(fptr, TSTRING, "ST_PROC", timestampBytes.data(), 
                      "Processing timestamp", &status)) {
        logMessage("Failed to write ST_PROC keyword", "orange");
        status = 0;
    }
    
    fits_close_file(fptr, &status);
    
    if (status != 0) {
        logMessage(QString("FITS errors occurred while writing metadata (status: %1)").arg(status), "orange");
        return false;
    }
    
    return true;
}

bool StellinaProcessor::updateProcessingStage(const QString &fitsPath, const QString &stage) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READWRITE, &status)) {
        return false;
    }
    
    // Update timestamp
    QString timestamp = QDateTime::currentDateTime().toString(Qt::ISODate);
    QByteArray timestampBytes = timestamp.toLocal8Bit();
    fits_update_key(fptr, TSTRING, "ST_PROC", timestampBytes.data(), 
                   "Processing timestamp", &status);
    
    fits_close_file(fptr, &status);
    return (status == 0);
}

bool StellinaProcessor::cleanExistingStellinaKeywords(const QString &fitsPath) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READWRITE, &status)) {
        return false;
    }
    
    bool result = cleanExistingStellinaKeywordsFITS(fptr);
    fits_close_file(fptr, &status);
    
    return result;
}

// ============================================================================
// FITS Image Operations
// ============================================================================

bool StellinaProcessor::loadFITSImage(const QString &fitsPath, cv::Mat &image) {
    fitsfile *fptr = nullptr;
    int status = 0;
    int naxis, bitpix;
    long naxes[2];
    
    // Open FITS file
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        logMessage(QString("Could not open FITS file: %1").arg(fitsPath), "red");
        return false;
    }
    
    // Get image dimensions
    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status)) {
        logMessage("Could not read FITS image parameters", "red");
        fits_close_file(fptr, &status);
        return false;
    }
    
    if (naxis != 2) {
        logMessage(QString("Unsupported FITS image dimensions: %1").arg(naxis), "red");
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Create OpenCV matrix
    int width = static_cast<int>(naxes[0]);
    int height = static_cast<int>(naxes[1]);
    
    if (bitpix == 16 || bitpix == -16) {
        // 16-bit image
        image = cv::Mat::zeros(height, width, CV_16U);
        
        long fpixel[2] = {1, 1};
        if (fits_read_pix(fptr, TUSHORT, fpixel, width * height, nullptr, 
                         image.data, nullptr, &status)) {
            logMessage("Could not read FITS image data", "red");
            fits_close_file(fptr, &status);
            return false;
        }
    } else if (bitpix == 32 || bitpix == -32) {
        // 32-bit image
        image = cv::Mat::zeros(height, width, CV_32F);
        
        long fpixel[2] = {1, 1};
        if (fits_read_pix(fptr, TFLOAT, fpixel, width * height, nullptr, 
                         image.data, nullptr, &status)) {
            logMessage("Could not read FITS image data", "red");
            fits_close_file(fptr, &status);
            return false;
        }
    } else {
        logMessage(QString("Unsupported FITS bit depth: %1").arg(bitpix), "red");
        fits_close_file(fptr, &status);
        return false;
    }
    
    fits_close_file(fptr, &status);
    
    if (m_debugMode) {
        logMessage(QString("Loaded FITS image: %1x%2, %3-bit")
                  .arg(width).arg(height).arg(bitpix), "gray");
    }
    
    return true;
}

bool StellinaProcessor::saveFITSImage(const QString &fitsPath, const cv::Mat &image) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    // Determine FITS parameters from OpenCV image
    long naxes[2] = {image.cols, image.rows};
    int bitpix;
    int datatype;
    
    if (image.type() == CV_16U) {
        bitpix = 16;
        datatype = TUSHORT;
    } else if (image.type() == CV_32F) {
        bitpix = -32;
        datatype = TFLOAT;
    } else {
        logMessage("Unsupported OpenCV image type for FITS export", "red");
        return false;
    }
    
    // Create FITS file
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_create_file(&fptr, pathBytes.data(), &status)) {
        logMessage(QString("Could not create FITS file: %1").arg(fitsPath), "red");
        return false;
    }
    
    // Create image HDU
    if (fits_create_img(fptr, bitpix, 2, naxes, &status)) {
        logMessage("Could not create FITS image HDU", "red");
        fits_close_file(fptr, &status);
        QFile::remove(fitsPath);
        return false;
    }
    
    // Write image data
    long fpixel[2] = {1, 1};
    if (fits_write_pix(fptr, datatype, fpixel, image.total(), 
                      const_cast<void*>(static_cast<const void*>(image.data)), &status)) {
        logMessage("Could not write FITS image data", "red");
        fits_close_file(fptr, &status);
        QFile::remove(fitsPath);
        return false;
    }
    
    fits_close_file(fptr, &status);
    
    if (status != 0) {
        logMessage(QString("FITS errors occurred while saving (status: %1)").arg(status), "red");
        QFile::remove(fitsPath);
        return false;
    }
    
    return true;
}

bool StellinaProcessor::copyFITSHeaders(const QString &sourcePath, const QString &destPath) {
    fitsfile *sourcePtr = nullptr, *destPtr = nullptr;
    int status = 0;
    
    // Open source file
    QByteArray sourceBytes = sourcePath.toLocal8Bit();
    if (fits_open_file(&sourcePtr, sourceBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Open destination file
    QByteArray destBytes = destPath.toLocal8Bit();
    if (fits_open_file(&destPtr, destBytes.data(), READWRITE, &status)) {
        fits_close_file(sourcePtr, &status);
        return false;
    }
    
    // Copy header keywords (excluding structural keywords)
    int nkeys;
    if (fits_get_hdrspace(sourcePtr, &nkeys, nullptr, &status) == 0) {
        char keyname[FLEN_KEYWORD], value[FLEN_VALUE], comment[FLEN_COMMENT];
        
        for (int i = 1; i <= nkeys; i++) {
            if (fits_read_keyn(sourcePtr, i, keyname, value, comment, &status) == 0) {
                // Skip structural keywords
                QString key(keyname);
                if (key.startsWith("SIMPLE") || key.startsWith("BITPIX") || 
                    key.startsWith("NAXIS") || key.startsWith("EXTEND")) {
                    continue;
                }
                
                // Copy the keyword
                fits_write_key_str(destPtr, keyname, value, comment, &status);
                status = 0; // Reset status for next iteration
            }
        }
    }
    
    fits_close_file(sourcePtr, &status);
    fits_close_file(destPtr, &status);
    
    return true;
}

// ============================================================================
// FITS Validation and Utilities
// ============================================================================

bool StellinaProcessor::validateFITSFile(const QString &fitsPath) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Check if it's a valid image
    int naxis;
    if (fits_get_img_dim(fptr, &naxis, &status) || naxis != 2) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    fits_close_file(fptr, &status);
    return (status == 0);
}

bool StellinaProcessor::hasStellinaMetadata(const QString &fitsPath) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Check for Stellina-specific keywords
    char value[FLEN_VALUE];
    bool hasStellina = (fits_read_key(fptr, TSTRING, "ST_ALT", value, nullptr, &status) == 0) ||
                      (fits_read_key(fptr, TSTRING, "ST_AZ", value, nullptr, &status) == 0);
    
    fits_close_file(fptr, &status);
    return hasStellina;
}

QString StellinaProcessor::getBayerPatternFromFits(const QString &fitsPath) {
    fitsfile *fptr = nullptr;
    int status = 0;
    char pattern[FLEN_VALUE];
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return "RGGB"; // Default
    }
    
    if (fits_read_key(fptr, TSTRING, "BAYERPAT", pattern, nullptr, &status) == 0) {
        fits_close_file(fptr, &status);
        return QString(pattern);
    }
    
    fits_close_file(fptr, &status);
    
    // Default heuristic based on filename
    QString filename = QFileInfo(fitsPath).fileName().toLower();
    if (filename.contains("r.fits") || filename.endsWith("r")) {
        return "BGGR";  // Reversed images
    }
    
    return "RGGB";  // Default
}

// ============================================================================
// Private Helper Functions
// ============================================================================

bool StellinaProcessor::cleanExistingStellinaKeywordsFITS(fitsfile *fptr) {
    int status = 0;
    
    // List of Stellina-specific keywords to remove
    const QStringList stellinaKeywords = {
        "ST_ALT", "ST_AZ", "ST_RA", "ST_DEC", "ST_PROC"
    };
    
    for (const QString &keyword : stellinaKeywords) {
        QByteArray keyBytes = keyword.toLocal8Bit();
        fits_delete_key(fptr, keyBytes.data(), &status);
        status = 0; // Reset status - deletion failure is not critical
    }
    
    return true;
}

// ============================================================================
// JSON Metadata Operations
// ============================================================================

/*
bool StellinaProcessor::readStellinaDataFromJSON(const QString &jsonPath, StellinaImageData &data) {
    QFile jsonFile(jsonPath);
    if (!jsonFile.open(QIODevice::ReadOnly)) {
        if (m_debugMode) {
            logMessage(QString("Could not open JSON file: %1").arg(jsonPath), "orange");
        }
        return false;
    }
    
    QJsonParseError parseError;
    QJsonDocument doc = QJsonDocument::fromJson(jsonFile.readAll(), &parseError);
    
    if (parseError.error != QJsonParseError::NoError) {
        logMessage(QString("JSON parse error in %1: %2").arg(jsonPath).arg(parseError.errorString()), "red");
        return false;
    }
    
    QJsonObject rootObj = doc.object();
    data.metadata = rootObj;
    
    // Extract essential coordinates
    if (rootObj.contains("altitude") && rootObj.contains("azimuth")) {
        data.altitude = rootObj["altitude"].toDouble();
        data.azimuth = rootObj["azimuth"].toDouble();
        data.hasValidCoordinates = true;
    }
    
    // Extract timing information
    if (rootObj.contains("dateObs")) {
        data.dateObs = rootObj["dateObs"].toString();
    }
    
    return true;
}
    int status = 0;
    int naxis, bitpix;
    long naxes[2];
    
    // Open FITS file
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        logMessage(QString("Could not open FITS file: %1").arg(fitsPath), "red");
        return false;
    }
    
    // Get image dimensions
    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status)) {
        logMessage("Could not read FITS image parameters", "red");
        fits_close_file(fptr, &status);
        return false;
    }
    
    if (naxis != 2) {
        logMessage(QString("Unsupported FITS image dimensions: %1").arg(naxis), "red");
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Create OpenCV matrix
    int width = static_cast<int>(naxes[0]);
    int height = static_cast<int>(naxes[1]);
    
    if (bitpix == 16 || bitpix == -16) {
        // 16-bit image
        image = cv::Mat::zeros(height, width, CV_16U);
        
        long fpixel[2] = {1, 1};
        if (fits_read_pix(fptr, TUSHORT, fpixel, width * height, nullptr, 
                         image.data, nullptr, &status)) {
            logMessage("Could not read FITS image data", "red");
            fits_close_file(fptr, &status);
            return false;
        }
    } else if (bitpix == 32 || bitpix == -32) {
        // 32-bit image
        image = cv::Mat::zeros(height, width, CV_32F);
        
        long fpixel[2] = {1, 1};
        if (fits_read_pix(fptr, TFLOAT, fpixel, width * height, nullptr, 
                         image.data, nullptr, &status)) {
            logMessage("Could not read FITS image data", "red");
            fits_close_file(fptr, &status);
            return false;
        }
    } else {
        logMessage(QString("Unsupported FITS bit depth: %1").arg(bitpix), "red");
        fits_close_file(fptr, &status);
        return false;
    }
    
    fits_close_file(fptr, &status);
    
    if (m_debugMode) {
        logMessage(QString("Loaded FITS image: %1x%2, %3-bit")
                  .arg(width).arg(height).arg(bitpix), "gray");
    }
    
    return true;
}

bool StellinaProcessor::saveFITSImage(const QString &fitsPath, const cv::Mat &image) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    // Determine FITS parameters from OpenCV image
    long naxes[2] = {image.cols, image.rows};
    int bitpix;
    int datatype;
    
    if (image.type() == CV_16U) {
        bitpix = 16;
        datatype = TUSHORT;
    } else if (image.type() == CV_32F) {
        bitpix = -32;
        datatype = TFLOAT;
    } else {
        logMessage("Unsupported OpenCV image type for FITS export", "red");
        return false;
    }
    
    // Create FITS file
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_create_file(&fptr, pathBytes.data(), &status)) {
        logMessage(QString("Could not create FITS file: %1").arg(fitsPath), "red");
        return false;
    }
    
    // Create image HDU
    if (fits_create_img(fptr, bitpix, 2, naxes, &status)) {
        logMessage("Could not create FITS image HDU", "red");
        fits_close_file(fptr, &status);
        QFile::remove(fitsPath);
        return false;
    }
    
    // Write image data
    long fpixel[2] = {1, 1};
    if (fits_write_pix(fptr, datatype, fpixel, image.total(), 
                      const_cast<void*>(static_cast<const void*>(image.data)), &status)) {
        logMessage("Could not write FITS image data", "red");
        fits_close_file(fptr, &status);
        QFile::remove(fitsPath);
        return false;
    }
    
    fits_close_file(fptr, &status);
    
    if (status != 0) {
        logMessage(QString("FITS errors occurred while saving (status: %1)").arg(status), "red");
        QFile::remove(fitsPath);
        return false;
    }
    
    return true;
}

bool StellinaProcessor::copyFITSHeaders(const QString &sourcePath, const QString &destPath) {
    fitsfile *sourcePtr = nullptr, *destPtr = nullptr;
    int status = 0;
    
    // Open source file
    QByteArray sourceBytes = sourcePath.toLocal8Bit();
    if (fits_open_file(&sourcePtr, sourceBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Open destination file
    QByteArray destBytes = destPath.toLocal8Bit();
    if (fits_open_file(&destPtr, destBytes.data(), READWRITE, &status)) {
        fits_close_file(sourcePtr, &status);
        return false;
    }
    
    // Copy header keywords (excluding structural keywords)
    int nkeys;
    if (fits_get_hdrspace(sourcePtr, &nkeys, nullptr, &status) == 0) {
        char keyname[FLEN_KEYWORD], value[FLEN_VALUE], comment[FLEN_COMMENT];
        
        for (int i = 1; i <= nkeys; i++) {
            if (fits_read_keyn(sourcePtr, i, keyname, value, comment, &status) == 0) {
                // Skip structural keywords
                QString key(keyname);
                if (key.startsWith("SIMPLE") || key.startsWith("BITPIX") || 
                    key.startsWith("NAXIS") || key.startsWith("EXTEND")) {
                    continue;
                }
                
                // Copy the keyword
                fits_write_key_str(destPtr, keyname, value, comment, &status);
                status = 0; // Reset status for next iteration
            }
        }
    }
    
    fits_close_file(sourcePtr, &status);
    fits_close_file(destPtr, &status);
    
    return true;
}


// ============================================================================
// FITS Validation and Utilities
// ============================================================================

bool StellinaProcessor::validateFITSFile(const QString &fitsPath) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Check if it's a valid image
    int naxis;
    if (fits_get_img_dim(fptr, &naxis, &status) || naxis != 2) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    fits_close_file(fptr, &status);
    return (status == 0);
}
 
bool StellinaProcessor::hasStellinaMetadata(const QString &fitsPath) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Check for Stellina-specific keywords
    char value[FLEN_VALUE];
    bool hasStellina = (fits_read_key(fptr, TSTRING, "ST_ALT", value, nullptr, &status) == 0) ||
                      (fits_read_key(fptr, TSTRING, "ST_AZ", value, nullptr, &status) == 0);
    
    fits_close_file(fptr, &status);
    return hasStellina;
}
 
QString StellinaProcessor::getBayerPatternFromFits(const QString &fitsPath) {
    fitsfile *fptr = nullptr;
    int status = 0;
    char pattern[FLEN_VALUE];
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return "RGGB"; // Default
    }
    
    if (fits_read_key(fptr, TSTRING, "BAYERPAT", pattern, nullptr, &status) == 0) {
        fits_close_file(fptr, &status);
        return QString(pattern);
    }
    
    fits_close_file(fptr, &status);
    
    // Default heuristic based on filename
    QString filename = QFileInfo(fitsPath).fileName().toLower();
    if (filename.contains("r.fits") || filename.endsWith("r")) {
        return "BGGR";  // Reversed images
    }
    
    return "RGGB";  // Default
}
 
// ============================================================================
// Private Helper Functions
// ============================================================================

bool StellinaProcessor::cleanExistingStellinaKeywordsFITS(fitsfile *fptr) {
    int status = 0;
    
    // List of Stellina-specific keywords to remove
    const QStringList stellinaKeywords = {
        "ST_ALT", "ST_AZ", "ST_RA", "ST_DEC", "ST_STAGE", "ST_PROC"
    };
    
    for (const QString &keyword : stellinaKeywords) {
        QByteArray keyBytes = keyword.toLocal8Bit();
        fits_delete_key(fptr, keyBytes.data(), &status);
        status = 0; // Reset status - deletion failure is not critical
    }
    
    return true;
}
 */
// ============================================================================
// JSON Metadata Operations
// ============================================================================

bool StellinaProcessor::readStellinaDataFromJSON(const QString &jsonPath, StellinaImageData &data) {
    QFile jsonFile(jsonPath);
    if (!jsonFile.open(QIODevice::ReadOnly)) {
        if (m_debugMode) {
            logMessage(QString("Could not open JSON file: %1").arg(jsonPath), "orange");
        }
        return false;
    }
    
    QJsonParseError parseError;
    QJsonDocument doc = QJsonDocument::fromJson(jsonFile.readAll(), &parseError);
    
    if (parseError.error != QJsonParseError::NoError) {
        logMessage(QString("JSON parse error in %1: %2").arg(jsonPath).arg(parseError.errorString()), "red");
        return false;
    }
    
    QJsonObject rootObj = doc.object();
    data.metadata = rootObj;
    
    // Extract essential coordinates
    if (rootObj.contains("altitude") && rootObj.contains("azimuth")) {
        data.altitude = rootObj["altitude"].toDouble();
        data.azimuth = rootObj["azimuth"].toDouble();
        data.hasValidCoordinates = true;
    }
    
    // Extract timing information
    if (rootObj.contains("dateObs")) {
        data.dateObs = rootObj["dateObs"].toString();
    }
    
    return true;
}
