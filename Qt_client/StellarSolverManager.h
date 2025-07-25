#include <QCoreApplication>
#include <QObject>
#include <QQueue>
#include <QJsonObject>
#include <iostream>
#include <fitsio.h>

// StellarSolver includes
#include <stellarsolver.h>
#include <parameters.h>
#include <structuredefinitions.h>
#include "SimplifiedXISFWriter.h"

struct SolveJob {
    QString filename;
    int jobId;
    QDateTime startTime;
    QString status = "PENDING";
    double ra = -1;
    double dec = -91;
    double pixelScale = 0;
    double solveTime = 0;
    double totalTime = 0;
};

class ParallelStellarSolver : public QObject
{
    Q_OBJECT

public:
    ParallelStellarSolver(int maxConcurrent = 8, QObject *parent = nullptr) 
        : QObject(parent), m_maxConcurrent(maxConcurrent), m_completedJobs(0) {
        
        // Setup shared parameters once
        setupCommonParameters();
    }
    
    ~ParallelStellarSolver() {
        // Clean up any remaining solvers
        for (auto* solver : m_activeSolvers.keys()) {
            if (solver && solver->isRunning()) {
                solver->abortAndWait();
            }
            solver->deleteLater();
        }
    }

    void addJob(const QString& filename) {
        SolveJob job;
        job.filename = filename;
        job.jobId = m_jobQueue.size() + 1;
        m_jobQueue.enqueue(job);
    }

    void startBatchSolving() {
        if (m_jobQueue.isEmpty()) {
            std::cerr << "No jobs to process!" << std::endl;
            QCoreApplication::quit();
            return;
        }

	m_totalJobs = m_jobQueue.size();
        std::cout << "=== Parallel Plate Solving ===" << std::endl;
        std::cout << "Total jobs: " << m_totalJobs << std::endl;
        std::cout << "Max concurrent: " << m_maxConcurrent << std::endl;
        std::cout << std::endl;

        m_batchStartTime = QDateTime::currentDateTime();
	    emit statusChanged(QString("Starting batch solve: %1 jobs").arg(m_totalJobs));
       
        // Start initial batch of jobs
        startNextJobs();
    }

    // Add this method to set the output directory (call from StellinaProcessor):
    void setOutputDirectory(const QString& outputDir) {
    m_outputDirectory = outputDir;
    }

signals:
    void progressChanged(int completed, int total);  // Main progress signal
    void statusChanged(const QString& message);      // Status updates
    void batchFinished(int successful, int failed);  // Completion signal

private slots:

// Replace your existing onSolverFinished() slot with this:
void onSolverFinished() {
    StellarSolver* solver = qobject_cast<StellarSolver*>(sender());
    if (!solver || !m_activeSolvers.contains(solver)) {
        return;
    }

    SolveJob& job = m_activeSolvers[solver];
    
    // Calculate times
    qint64 totalMs = job.startTime.msecsTo(QDateTime::currentDateTime());
    job.totalTime = totalMs / 1000.0;

    // Process results
    if (solver->solvingDone() && solver->hasWCSData()) {
        FITSImage::Solution solution = solver->getSolution();
        job.status = "SUCCESS";
        job.ra = solution.ra;
        job.dec = solution.dec;
        job.pixelScale = solution.pixscale / 2.0; // Correct for 2x downsampling
        job.solveTime = job.totalTime;
        
        // Write WCS results to output FITS file
        QString inputFile = job.filename;
        QString baseName = QFileInfo(inputFile).baseName();
        
        if (m_useXISF) {
            QString outputFile = QString("%1/plate_solved_%2.xisf").arg(m_outputDirectory).arg(baseName);
            if (writeWCSToXISF(inputFile, outputFile, solution)) {
                std::cout << "    XISF written to: " << QFileInfo(outputFile).fileName().toStdString() << std::endl;
            } else {
                std::cout << "    Failed to write XISF file" << std::endl;
            }
        } else {
            QString outputFile = QString("%1/plate_solved_%2.fits").arg(m_outputDirectory).arg(baseName);
	    if (writeWCSToFITS(inputFile, outputFile, solution)) {
		std::cout << "[Job " << job.jobId << "] ✓ SUCCESS: " << baseName.toStdString() << std::endl;
		std::cout << "    RA: " << job.ra << "°, Dec: " << job.dec << "°" << std::endl;
		std::cout << "    Pixel scale: " << job.pixelScale << " arcsec/pix" << std::endl;
		std::cout << "    WCS written to: " << QFileInfo(outputFile).fileName().toStdString() << std::endl;
		std::cout << "    Total time: " << std::fixed << std::setprecision(2) << job.totalTime << "s" << std::endl;
	    } else {
		std::cout << "[Job " << job.jobId << "] ⚠ SOLVED but failed to write WCS: " << baseName.toStdString() << std::endl;
		job.status = "SOLVED_NO_WRITE";
	    }
        }
        
    } else {
        job.status = "FAILED";
        std::cout << "[Job " << job.jobId << "] ✗ FAILED: " << QFileInfo(job.filename).baseName().toStdString();
        if (solver->failed()) {
            std::cout << " (solver failed)";
        } else if (!solver->solvingDone()) {
            std::cout << " (incomplete)";
        } else {
            std::cout << " (no WCS data)";
        }
        std::cout << std::endl;
        std::cout << "    Total time: " << std::fixed << std::setprecision(2) << job.totalTime << "s" << std::endl;
    }

    // Store completed job
    m_completedJobs = m_totalJobs - m_jobQueue.size();
    emit progressChanged(m_completedJobs, m_totalJobs);
    m_results.append(job);
    
    // Clean up solver
    m_activeSolvers.remove(solver);
    solver->deleteLater();

    // Progress update
    std::cout << "Progress: " << m_completedJobs << "/" << m_totalJobs << " completed" << std::endl;
    std::cout << std::endl;

    // Check if all jobs are complete or start next job
    if (m_completedJobs >= m_totalJobs) {
        finishBatch();
    } else {
        startNextJobs();
    }
}
    
private:
    void setupCommonParameters() {
        // Get built-in profiles
        QList<Parameters> profiles = StellarSolver::getBuiltInProfiles();
        if (profiles.isEmpty()) {
            std::cerr << "ERROR: No parameter profiles available" << std::endl;
            return;
        }

        m_commonParams = profiles.at(0); // Use first profile as base

        // Configure parameters
        m_commonParams.multiAlgorithm = SSolver::MULTI_AUTO;
        m_commonParams.search_radius = 2.0;      // Match --radius=2 
        m_commonParams.minwidth = 0.1;
        m_commonParams.maxwidth = 10.0;
        m_commonParams.resort = true;
        m_commonParams.autoDownsample = false;
        m_commonParams.downsample = 1;           // We do our own downsampling
        m_commonParams.inParallel = true;
        m_commonParams.solverTimeLimit = 300;

        // Star extraction parameters
        m_commonParams.initialKeep = 2000;
        m_commonParams.keepNum = 500;
        m_commonParams.r_min = 1.0;
        m_commonParams.removeBrightest = 0;
        m_commonParams.removeDimmest = 50;
        m_commonParams.saturationLimit = 65000;
        m_commonParams.minarea = 5;
        m_commonParams.threshold_offset = 0;
        m_commonParams.threshold_bg_multiple = 2.0;

        // Find index files
        m_indexPaths = findIndexFiles();
        if (m_indexPaths.isEmpty()) {
            std::cerr << "ERROR: No astrometry index files found!" << std::endl;
        }
    }

    void startNextJobs() {
        while (m_activeSolvers.size() < m_maxConcurrent && !m_jobQueue.isEmpty()) {
            SolveJob job = m_jobQueue.dequeue();
            startSingleJob(job);
        }
    }

    void startSingleJob(SolveJob job) {
        // Create new solver
        StellarSolver* solver = new StellarSolver(this);
        
        // Connect signals
        connect(solver, &StellarSolver::finished, this, &ParallelStellarSolver::onSolverFinished);
        
        // Configure solver
        solver->setProperty("ProcessType", SSolver::SOLVE);
        solver->setProperty("ExtractorType", SSolver::EXTRACTOR_INTERNAL);
        solver->setProperty("SolverType", SSolver::SOLVER_STELLARSOLVER);
        solver->setIndexFolderPaths(m_indexPaths);
        solver->setParameters(m_commonParams);

        // Load FITS image
        if (!loadFITSImage(solver, job.filename)) {
            job.status = "SKIPPED_NO_COORDS";
            m_completedJobs++;
            m_results.append(job);
            solver->deleteLater();
            
            std::cout << "[Job " << job.jobId << "] ⚬ SKIPPED: " << QFileInfo(job.filename).baseName().toStdString() 
                      << " (no STELLRA/STELLDEC keywords)" << std::endl;
            
            if (m_completedJobs >= m_totalJobs) {
                finishBatch();
            } else {
                startNextJobs();
            }
            return;
        }

        // Start timing and add to active jobs
        job.startTime = QDateTime::currentDateTime();
        m_activeSolvers[solver] = job;
        
        std::cout << "[Job " << job.jobId << "] Starting: " << QFileInfo(job.filename).baseName().toStdString() << std::endl;
        
        // Start solving
        solver->start();
    }

    bool loadFITSImage(StellarSolver* solver, const QString& fitsFile) {
        fitsfile *fptr;
        int status = 0;
        int naxis, bitpix;
        long naxes[2];
        
        // Open FITS file
        if (fits_open_file(&fptr, fitsFile.toLocal8Bit().data(), READONLY, &status)) {
            return false;
        }

	// Try to find coordinate information in priority order:
	// 1. STELLRA/STELLDEC (preferred - from Stellina processing)
	// 2. CRVAL1/CRVAL2 (fallback - from existing WCS)

	double coordinateRA = -1, coordinateDec = -91;
	QString coordinateSource = "unknown";

	// First, try STELLRA/STELLDEC
	char comment[FLEN_COMMENT];
	status = 0;
	bool hasSTELLRA = (fits_read_key(fptr, TDOUBLE, "STELLRA", &coordinateRA, comment, &status) == 0);
	status = 0;
	bool hasSTELLDEC = (fits_read_key(fptr, TDOUBLE, "STELLDEC", &coordinateDec, comment, &status) == 0);

	if (hasSTELLRA && hasSTELLDEC && coordinateRA >= 0 && coordinateRA <= 360 && 
        coordinateDec >= -90 && coordinateDec <= 90) {
	    coordinateSource = "STELLRA/STELLDEC";
	} else {
	    // Fallback to WCS coordinates (CRVAL1/CRVAL2)
	    status = 0;
	    bool hasCRVAL1 = (fits_read_key(fptr, TDOUBLE, "CRVAL1", &coordinateRA, comment, &status) == 0);
	    status = 0;
	    bool hasCRVAL2 = (fits_read_key(fptr, TDOUBLE, "CRVAL2", &coordinateDec, comment, &status) == 0);

	    if (hasCRVAL1 && hasCRVAL2 && coordinateRA >= 0 && coordinateRA <= 360 && 
		coordinateDec >= -90 && coordinateDec <= 90 && (coordinateRA != 0 && coordinateDec != 0)) {
		coordinateSource = "CRVAL1/CRVAL2 (WCS)";
	    } else {
		// No valid coordinates found
		fits_close_file(fptr, &status);
		std::cout << "[SKIP] No valid coordinates in " << QFileInfo(fitsFile).baseName().toStdString() 
			  << " (checked STELLRA/STELLDEC and CRVAL1/CRVAL2)" << std::endl;
		return false;
	    }
	}

	// Set coordinate hints for the solver
	solver->setSearchPositionInDegrees(coordinateRA, coordinateDec);
	
        // Get image dimensions
        status = 0;
        if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) || naxis != 2) {
            fits_close_file(fptr, &status);
            return false;
        }
        
        int width = static_cast<int>(naxes[0]);
        int height = static_cast<int>(naxes[1]);
        long npixels = width * height;
        
        // Allocate and read image data
        std::vector<float> imageData(npixels);
        long firstPix[2] = {1, 1};
        if (fits_read_pix(fptr, TFLOAT, firstPix, npixels, nullptr, 
                         imageData.data(), nullptr, &status)) {
            fits_close_file(fptr, &status);
            return false;
        }
        
        fits_close_file(fptr, &status);
        
        // Calculate statistics
        float minVal = *std::min_element(imageData.begin(), imageData.end());
        float maxVal = *std::max_element(imageData.begin(), imageData.end());
        
        // Downsample with Bayer-aware processing
        const int downsampleFactor = 2;
        int outputWidth = width / downsampleFactor;
        int outputHeight = height / downsampleFactor;
        
        std::vector<uint8_t> buffer(outputWidth * outputHeight);
        float range = maxVal - minVal;
        
        if (range > 0) {
            for (int y = 0; y < outputHeight; y++) {
                for (int x = 0; x < outputWidth; x++) {
                    int srcY = y * downsampleFactor;
                    int srcX = x * downsampleFactor;
                    
                    if (srcY + 1 >= height || srcX + 1 >= width) continue;
                    
                    // 2x2 Bayer pattern with proper weighting
                    float r = imageData[srcY * width + srcX];
                    float g1 = imageData[srcY * width + srcX + 1];
                    float g2 = imageData[(srcY + 1) * width + srcX];
                    float b = imageData[(srcY + 1) * width + srcX + 1];
                    
                    float avg = (r + (g1 + g2) * 0.5f + b) / 3.0f;
                    float normalized = (avg - minVal) / range;
                    buffer[y * outputWidth + x] = static_cast<uint8_t>(std::clamp(normalized * 255.0f, 0.0f, 255.0f));
                }
            }
        } else {
            std::fill(buffer.begin(), buffer.end(), 128);
        }
        
        // Create statistics
        FITSImage::Statistic stats{};
        stats.width = outputWidth;
        stats.height = outputHeight;
        stats.channels = 1;
        stats.dataType = TBYTE;
        stats.bytesPerPixel = 1;
        
        for (int i = 0; i < 3; i++) {
            stats.min[i] = (i == 0) ? minVal : 0.0;
            stats.max[i] = (i == 0) ? maxVal : 0.0;
            stats.mean[i] = (i == 0) ? (minVal + maxVal) / 2.0 : 0.0;
            stats.stddev[i] = 0.0;
            stats.median[i] = (i == 0) ? (minVal + maxVal) / 2.0 : 0.0;
        }
        stats.SNR = 1.0;
        
        // Store buffer in solver-specific storage to keep it alive
        m_imageBuffers[solver] = std::move(buffer);
        
        // Load into solver
        return solver->loadNewImageBuffer(stats, m_imageBuffers[solver].data());
    }

    void finishBatch() {
        qint64 totalBatchTime = m_batchStartTime.msecsTo(QDateTime::currentDateTime());
        
        std::cout << "=== Batch Complete ===" << std::endl;
        std::cout << "Total time: " << (totalBatchTime / 1000.0) << " seconds" << std::endl;
        std::cout << std::endl;
        
        // Print summary
        int successful = 0;
        int failed = 0;
        int skipped = 0;
        double totalSolveTime = 0;
        
        for (const auto& result : m_results) {
            if (result.status == "SUCCESS") {
                successful++;
                totalSolveTime += result.totalTime;
            } else if (result.status == "SKIPPED_NO_COORDS") {
                skipped++;
            } else {
                failed++;
            }
        }
        
        std::cout << "Summary:" << std::endl;
        std::cout << "  Successful: " << successful << "/" << m_totalJobs << std::endl;
        std::cout << "  Failed: " << failed << "/" << m_totalJobs << std::endl;
        std::cout << "  Skipped (no coordinates): " << skipped << "/" << m_totalJobs << std::endl;
        if (successful > 0) {
            std::cout << "  Average solve time: " << std::fixed << std::setprecision(2) 
                      << (totalSolveTime / successful) << "s" << std::endl;
            std::cout << "  Throughput (successful): " << std::fixed << std::setprecision(2) 
                      << (successful / (totalBatchTime / 1000.0)) << " solves/second" << std::endl;
        }
	
	emit batchFinished(successful, failed);
    }

    QStringList findIndexFiles() {
        QStringList indexPaths;
        QStringList searchPaths = {
            "/usr/local/astrometry/data",
            "/opt/homebrew/share/astrometry", 
            "/usr/local/share/astrometry"
        };
        
        for (const QString &path : searchPaths) {
            QDir indexDir(path);
            if (indexDir.exists()) {
                QStringList filters;
                filters << "index-*.fits";
                QFileInfoList indexFiles = indexDir.entryInfoList(filters, QDir::Files);
                
                if (!indexFiles.isEmpty()) {
                    indexPaths.append(path);
                    break;
                }
            }
        }
        
        return indexPaths;
    }

    // Optional: Handle fileio log messages
    void handleFileIOLog(const QString& message) {
	// You can process fileio log messages here if needed
	// std::cout << "FileIO: " << message.toStdString() << std::endl;
    }

// =============================================================================
// INTEGRATION WITH YOUR PIPELINE
// =============================================================================

    // Enhanced plate solving with XISF output
    bool writeWCSToXISF(const QString& inputPath, 
					      const QString& outputPath, 
                        const FITSImage::Solution& solution)
     {

	// Read original FITS image data (you already have this code)
	fitsfile *fptr = nullptr;
	int status = 0;

	QByteArray inputPathBytes = inputPath.toLocal8Bit();
	if (fits_open_file(&fptr, inputPathBytes.data(), READONLY, &status)) {
	    return false;
	}

	// Get image dimensions
	long naxes[2];
	if (fits_get_img_size(fptr, 2, naxes, &status)) {
	    fits_close_file(fptr, &status);
	    return false;
	}

	int width = naxes[0];
	int height = naxes[1];
	long totalPixels = width * height;

	// Read pixel data
	std::vector<float> pixels(totalPixels);
	if (fits_read_img(fptr, TFLOAT, 1, totalPixels, nullptr, pixels.data(), nullptr, &status)) {
	    fits_close_file(fptr, &status);
	    return false;
	}

	fits_close_file(fptr, &status);

	// Create XISF file
	SimplifiedXISFWriter writer(outputPath);

	// Add the image
	if (!writer.addImage("main", pixels.data(), width, height, 1)) {
	    return false;
	}

	// Add WCS data
	writer.addWCSData(solution);

	// Add processing history
	QJsonObject solvingParams;
	solvingParams["solver"] = "StellarSolver";
	solvingParams["input_file"] = QFileInfo(inputPath).fileName();
	solvingParams["solve_time"] = QString::number(m_results[0].totalTime, 'f', 2);
	solvingParams["pixel_scale"] = solution.pixscale;
	writer.addProcessingHistory("plate_solving", solvingParams);

	// Add file information
	writer.addProperty("File:OriginalPath", "String", inputPath, "Original FITS file path");
	writer.addProperty("File:ProcessingDate", "String", QDateTime::currentDateTimeUtc().toString(Qt::ISODate), "Processing date");

	return writer.write();
    }
    
    bool writeWCSToFITS(const QString& inputPath,
					      const QString& outputPath, 
					      const FITSImage::Solution& solution) {
	fitsfile *inputFits = nullptr, *outputFits = nullptr;
	int status = 0;

	// Open input FITS file
	QByteArray inputPathBytes = inputPath.toLocal8Bit();
	if (fits_open_file(&inputFits, inputPathBytes.data(), READONLY, &status)) {
	    std::cout << "Failed to open input FITS: " << inputPath.toStdString() 
		      << " (status: " << status << ")" << std::endl;
	    return false;
	}

	// Create output FITS file (force overwrite with !)
	QByteArray outputPathBytes = QString("!%1").arg(outputPath).toLocal8Bit();
	if (fits_create_file(&outputFits, outputPathBytes.data(), &status)) {
	    std::cout << "Failed to create output FITS: " << outputPath.toStdString() 
		      << " (status: " << status << ")" << std::endl;
	    fits_close_file(inputFits, &status);
	    return false;
	}

	// Copy entire input file to output (header + data)
	if (fits_copy_file(inputFits, outputFits, 1, 1, 1, &status)) {
	    std::cout << "Failed to copy FITS file (status: " << status << ")" << std::endl;
	    fits_close_file(inputFits, &status);
	    fits_close_file(outputFits, &status);
	    return false;
	}

	// Close input file (no longer needed)
	fits_close_file(inputFits, &status);

	// Get image dimensions for reference pixel calculation
	long naxes[2];
	if (fits_get_img_size(outputFits, 2, naxes, &status)) {
	    std::cout << "Warning: Could not get image dimensions for WCS reference pixel" << std::endl;
	    naxes[0] = 1920; // Default fallback
	    naxes[1] = 1080;
	    status = 0; // Continue anyway
	}

	// Write WCS keywords to the output file
	// Reference point (RA/Dec of reference pixel in degrees)
	double crval1 = solution.ra;
	double crval2 = solution.dec;
	fits_write_key(outputFits, TDOUBLE, "CRVAL1", &crval1, "Reference RA (degrees)", &status);
	fits_write_key(outputFits, TDOUBLE, "CRVAL2", &crval2, "Reference Dec (degrees)", &status);

	// Reference pixel (center of image, 1-indexed)
	double crpix1 = naxes[0] / 2.0 + 0.5;
	double crpix2 = naxes[1] / 2.0 + 0.5;
	fits_write_key(outputFits, TDOUBLE, "CRPIX1", &crpix1, "Reference pixel X", &status);
	fits_write_key(outputFits, TDOUBLE, "CRPIX2", &crpix2, "Reference pixel Y", &status);

	// Pixel scale converted to degrees/pixel (solution.pixscale is in arcsec/pixel)
	double pixscale_deg = solution.pixscale / 3600.0;

	// CD matrix (assuming square pixels, no rotation initially)
	double cd11 = -pixscale_deg;  // Negative for standard orientation
	double cd12 = 0.0;
	double cd21 = 0.0;
	double cd22 = pixscale_deg;

	// Apply orientation if available
	if (solution.orientation != 0.0) {
	    double orient_rad = solution.orientation * M_PI / 180.0;
	    double cos_orient = cos(orient_rad);
	    double sin_orient = sin(orient_rad);

	    cd11 = -pixscale_deg * cos_orient;
	    cd12 = pixscale_deg * sin_orient;
	    cd21 = pixscale_deg * sin_orient;
	    cd22 = pixscale_deg * cos_orient;
	}

	fits_write_key(outputFits, TDOUBLE, "CD1_1", &cd11, "Coordinate matrix element", &status);
	fits_write_key(outputFits, TDOUBLE, "CD1_2", &cd12, "Coordinate matrix element", &status);
	fits_write_key(outputFits, TDOUBLE, "CD2_1", &cd21, "Coordinate matrix element", &status);
	fits_write_key(outputFits, TDOUBLE, "CD2_2", &cd22, "Coordinate matrix element", &status);

	// Coordinate types (TAN projection)
	char ctype1[] = "RA---TAN";
	char ctype2[] = "DEC--TAN";
	char* ctype1_ptr = ctype1;
	char* ctype2_ptr = ctype2;
	fits_write_key(outputFits, TSTRING, "CTYPE1", &ctype1_ptr, "Coordinate type", &status);
	fits_write_key(outputFits, TSTRING, "CTYPE2", &ctype2_ptr, "Coordinate type", &status);

	// Coordinate units
	char cunit1[] = "deg";
	char cunit2[] = "deg";
	char* cunit1_ptr = cunit1;
	char* cunit2_ptr = cunit2;
	fits_write_key(outputFits, TSTRING, "CUNIT1", &cunit1_ptr, "Coordinate unit", &status);
	fits_write_key(outputFits, TSTRING, "CUNIT2", &cunit2_ptr, "Coordinate unit", &status);

	// Additional WCS information
	double pixscale_arcsec = solution.pixscale;
	fits_write_key(outputFits, TDOUBLE, "PIXSCALE", &pixscale_arcsec, "Pixel scale (arcsec/pixel)", &status);

	double orientation = solution.orientation;
	fits_write_key(outputFits, TDOUBLE, "ORIENTAT", &orientation, "Position angle (degrees)", &status);

	// Field dimensions (using correct field names from structure)
	double fieldw = solution.fieldWidth;   // Field width in arcminutes
	double fieldh = solution.fieldHeight;  // Field height in arcminutes
	fits_write_key(outputFits, TDOUBLE, "FIELDW", &fieldw, "Field width (arcmin)", &status);
	fits_write_key(outputFits, TDOUBLE, "FIELDH", &fieldh, "Field height (arcmin)", &status);

	// Processing information
	char solver[] = "StellarSolver";
	char* solver_ptr = solver;
	fits_write_key(outputFits, TSTRING, "SOLVER", &solver_ptr, "Plate solving software", &status);

	// Parity information
	QString parityText = (solution.parity == FITSImage::NEGATIVE) ? "negative" : "positive";
	QByteArray parityBytes = parityText.toLocal8Bit();
	char* parityPtr = parityBytes.data();
	fits_write_key(outputFits, TSTRING, "PARITY", &parityPtr, "Image parity", &status);

	// Solution errors if available
	if (solution.raError > 0 && solution.decError > 0) {
	    double raError = solution.raError;
	    double decError = solution.decError;
	    fits_write_key(outputFits, TDOUBLE, "RAERR", &raError, "RA error (arcsec)", &status);
	    fits_write_key(outputFits, TDOUBLE, "DECERR", &decError, "Dec error (arcsec)", &status);
	}

	// Add processing history
	QString history = QString("Plate solved: RA=%1°, Dec=%2°, Scale=%3\"/pix, Orient=%4°, Parity=%5")
			    .arg(solution.ra, 0, 'f', 6)
			    .arg(solution.dec, 0, 'f', 6)
			    .arg(solution.pixscale, 0, 'f', 3)
			    .arg(solution.orientation, 0, 'f', 2)
			    .arg(parityText);
	QByteArray historyBytes = history.toLocal8Bit();
	fits_write_history(outputFits, historyBytes.data(), &status);

	// Close output file
	fits_close_file(outputFits, &status);

	if (status != 0) {
	    std::cout << "FITS error while writing WCS data (status: " << status << ")" << std::endl;
	    return false;
	}

	return true;
    }

private:
    int m_maxConcurrent;
    int m_totalJobs = 0;
    std::atomic<int> m_completedJobs;
    
    Parameters m_commonParams;
    QStringList m_indexPaths;
    
    QQueue<SolveJob> m_jobQueue;
    QHash<StellarSolver*, SolveJob> m_activeSolvers;
    QHash<StellarSolver*, std::vector<uint8_t>> m_imageBuffers; // Keep buffers alive per solver
    QList<SolveJob> m_results;
    
    QDateTime m_batchStartTime;
    QString m_outputDirectory;  // Directory where solved FITS files will be written
    bool m_useXISF = true;
};
