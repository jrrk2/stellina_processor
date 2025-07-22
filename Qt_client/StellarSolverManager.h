#include <QCoreApplication>
#include <QObject>
#include <QQueue>
#include <iostream>
#include <fitsio.h>

// StellarSolver includes
#include <stellarsolver.h>
#include <parameters.h>
#include <structuredefinitions.h>

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

signals:
    void progressChanged(int completed, int total);  // Main progress signal
    void statusChanged(const QString& message);      // Status updates
    void batchFinished(int successful, int failed);  // Completion signal

private slots:
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
            job.solveTime = job.totalTime; // Approximate since we don't have internal timing
            
            std::cout << "[Job " << job.jobId << "] ✓ SUCCESS: " << QFileInfo(job.filename).baseName().toStdString() << std::endl;
            std::cout << "    RA: " << job.ra << "°, Dec: " << job.dec << "°" << std::endl;
            std::cout << "    Pixel scale: " << job.pixelScale << " arcsec/pix" << std::endl;
            std::cout << "    Total time: " << std::fixed << std::setprecision(2) << job.totalTime << "s" << std::endl;
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
	    // progressUpdated(m_completedJobs);
        std::cout << "Progress: " << m_completedJobs << "/" << m_totalJobs << " completed" << std::endl;
        std::cout << std::endl;

        // Start next job if available
        startNextJobs();

        // Check if all jobs complete
        if (m_completedJobs >= m_totalJobs) {
            finishBatch();
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
        
        // Check for required STELLRA and STELLDEC keywords - skip if missing
        double fitsRA = -1, fitsDec = -91;
        char comment[FLEN_COMMENT];
        
        status = 0;
        bool hasSTELLRA = (fits_read_key(fptr, TDOUBLE, "STELLRA", &fitsRA, comment, &status) == 0);
        status = 0;
        bool hasSTELLDEC = (fits_read_key(fptr, TDOUBLE, "STELLDEC", &fitsDec, comment, &status) == 0);
        
        // Skip files without both STELLRA and STELLDEC
        if (!hasSTELLRA || !hasSTELLDEC || fitsRA < 0 || fitsRA > 360 || fitsDec < -90 || fitsDec > 90) {
            fits_close_file(fptr, &status);
            return false;  // This will trigger LOAD_FAILED handling
        }
        
        // Set coordinate hints
        solver->setSearchPositionInDegrees(fitsRA, fitsDec);
        
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
};
