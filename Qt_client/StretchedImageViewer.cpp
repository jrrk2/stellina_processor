#pragma once

#include <QMainWindow>
#include <QMenuBar>
#include <QLabel>
#include <QScrollArea>
#include <QSlider>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QPushButton>
#include <QFileDialog>
#include <QMessageBox>
#include <QProgressBar>
#include <QStatusBar>
#include <fitsio.h>
#include <cmath>
#include <algorithm>

class StretchedImageViewer : public QMainWindow
{
    Q_OBJECT

public:
    explicit StretchedImageViewer(QWidget *parent = nullptr);
    ~StretchedImageViewer();

    // Public methods to integrate with your main application
    void loadStackedImage(const QString &fitsPath);
    void showViewer();

private slots:
    void onLoadImage();
    void onStretchTypeChanged();
    void onStretchParametersChanged();
    void onZoomChanged(int value);
    void onResetZoom();
    void onSaveStretchedImage();
    void onAutoStretch();

private:
    // UI Components
    QWidget *m_centralWidget;
    QScrollArea *m_scrollArea;
    QLabel *m_imageLabel;
    
    // Controls
    QGroupBox *m_stretchControlsGroup;
    QComboBox *m_stretchTypeCombo;
    QDoubleSpinBox *m_gammaSpinBox;
    QDoubleSpinBox *m_blackPointSpinBox;
    QDoubleSpinBox *m_whitePointSpinBox;
    QSlider *m_zoomSlider;
    QLabel *m_zoomLabel;
    QPushButton *m_loadButton;
    QPushButton *m_saveButton;
    QPushButton *m_autoStretchButton;
    QPushButton *m_resetZoomButton;
    QProgressBar *m_progressBar;
    
    // Image data
    QString m_currentImagePath;
    QImage m_originalImage;
    QImage m_stretchedImage;
    std::vector<float> m_imageData;
    long m_width, m_height;
    float m_minValue, m_maxValue;
    double m_currentZoom;

    // Methods
    void setupUI();
    void setupMenuBar();
    void connectSignals();
    bool loadFITSImage(const QString &filePath);
    void calculateImageStatistics();
    void applyStretch();
    void updateImageDisplay();
    void updateZoom();
    QImage convertToQImage(const std::vector<float> &data, long width, long height);
    
    // Stretch algorithms
    std::vector<float> applyLinearStretch(const std::vector<float> &data);
    std::vector<float> applyGammaStretch(const std::vector<float> &data);
    std::vector<float> applyLogStretch(const std::vector<float> &data);
    std::vector<float> applySinhStretch(const std::vector<float> &data);
    std::vector<float> applyAutoStretch(const std::vector<float> &data);
    
    // Utility methods
    void logMessage(const QString &message, const QString &color = "black");
    float percentile(const std::vector<float> &data, float percent);
};

// Implementation
StretchedImageViewer::StretchedImageViewer(QWidget *parent)
    : QMainWindow(parent)
    , m_centralWidget(nullptr)
    , m_scrollArea(nullptr)
    , m_imageLabel(nullptr)
    , m_width(0)
    , m_height(0)
    , m_minValue(0.0f)
    , m_maxValue(0.0f)
    , m_currentZoom(1.0)
{
    setWindowTitle("Stretched Stacked Image Viewer");
    setMinimumSize(800, 600);
    setupUI();
    setupMenuBar();
    connectSignals();
}

StretchedImageViewer::~StretchedImageViewer()
{
}

void StretchedImageViewer::setupUI()
{
    m_centralWidget = new QWidget(this);
    setCentralWidget(m_centralWidget);
    
    // Main layout
    QHBoxLayout *mainLayout = new QHBoxLayout(m_centralWidget);
    
    // Left panel for controls
    QVBoxLayout *controlsLayout = new QVBoxLayout();
    controlsLayout->setSpacing(10);
    
    // Load button
    m_loadButton = new QPushButton("Load Stacked Image");
    m_loadButton->setMinimumHeight(30);
    controlsLayout->addWidget(m_loadButton);
    
    // Stretch controls group
    m_stretchControlsGroup = new QGroupBox("Stretch Parameters");
    QVBoxLayout *stretchLayout = new QVBoxLayout(m_stretchControlsGroup);
    
    // Stretch type selection
    QLabel *stretchTypeLabel = new QLabel("Stretch Type:");
    m_stretchTypeCombo = new QComboBox();
    m_stretchTypeCombo->addItems({"Linear", "Gamma", "Logarithmic", "Sinh", "Auto"});
    m_stretchTypeCombo->setCurrentText("Auto");
    
    stretchLayout->addWidget(stretchTypeLabel);
    stretchLayout->addWidget(m_stretchTypeCombo);
    
    // Gamma parameter
    QLabel *gammaLabel = new QLabel("Gamma:");
    m_gammaSpinBox = new QDoubleSpinBox();
    m_gammaSpinBox->setRange(0.1, 5.0);
    m_gammaSpinBox->setValue(2.2);
    m_gammaSpinBox->setSingleStep(0.1);
    m_gammaSpinBox->setDecimals(2);
    
    stretchLayout->addWidget(gammaLabel);
    stretchLayout->addWidget(m_gammaSpinBox);
    
    // Black point
    QLabel *blackLabel = new QLabel("Black Point (%):");
    m_blackPointSpinBox = new QDoubleSpinBox();
    m_blackPointSpinBox->setRange(0.0, 50.0);
    m_blackPointSpinBox->setValue(0.1);
    m_blackPointSpinBox->setSingleStep(0.1);
    m_blackPointSpinBox->setDecimals(2);
    
    stretchLayout->addWidget(blackLabel);
    stretchLayout->addWidget(m_blackPointSpinBox);
    
    // White point
    QLabel *whiteLabel = new QLabel("White Point (%):");
    m_whitePointSpinBox = new QDoubleSpinBox();
    m_whitePointSpinBox->setRange(50.0, 100.0);
    m_whitePointSpinBox->setValue(99.9);
    m_whitePointSpinBox->setSingleStep(0.1);
    m_whitePointSpinBox->setDecimals(2);
    
    stretchLayout->addWidget(whiteLabel);
    stretchLayout->addWidget(m_whitePointSpinBox);
    
    controlsLayout->addWidget(m_stretchControlsGroup);
    
    // Auto stretch button
    m_autoStretchButton = new QPushButton("Auto Stretch");
    controlsLayout->addWidget(m_autoStretchButton);
    
    // Zoom controls
    QGroupBox *zoomGroup = new QGroupBox("Zoom");
    QVBoxLayout *zoomLayout = new QVBoxLayout(zoomGroup);
    
    m_zoomSlider = new QSlider(Qt::Horizontal);
    m_zoomSlider->setRange(10, 500); // 10% to 500%
    m_zoomSlider->setValue(100);
    
    m_zoomLabel = new QLabel("100%");
    m_zoomLabel->setAlignment(Qt::AlignCenter);
    
    m_resetZoomButton = new QPushButton("Reset Zoom");
    
    zoomLayout->addWidget(m_zoomLabel);
    zoomLayout->addWidget(m_zoomSlider);
    zoomLayout->addWidget(m_resetZoomButton);
    
    controlsLayout->addWidget(zoomGroup);
    
    // Save button
    m_saveButton = new QPushButton("Save Stretched Image");
    m_saveButton->setEnabled(false);
    controlsLayout->addWidget(m_saveButton);
    
    controlsLayout->addStretch();
    
    // Right panel for image display
    m_scrollArea = new QScrollArea();
    m_imageLabel = new QLabel();
    m_imageLabel->setAlignment(Qt::AlignCenter);
    m_imageLabel->setStyleSheet("border: 1px solid gray;");
    m_imageLabel->setText("Load an image to begin");
    m_scrollArea->setWidget(m_imageLabel);
    m_scrollArea->setWidgetResizable(true);
    
    // Add to main layout
    mainLayout->addLayout(controlsLayout, 0);
    mainLayout->addWidget(m_scrollArea, 1);
    
    // Status bar
    m_progressBar = new QProgressBar();
    m_progressBar->setVisible(false);
    statusBar()->addPermanentWidget(m_progressBar);
    statusBar()->showMessage("Ready");
}

void StretchedImageViewer::setupMenuBar()
{
    QMenuBar *menuBar = this->menuBar();
    
    // File menu
    QMenu *fileMenu = menuBar->addMenu("File");
    
    QAction *openAction = fileMenu->addAction("Open Image...");
    openAction->setShortcut(QKeySequence::Open);
    connect(openAction, &QAction::triggered, this, &StretchedImageViewer::onLoadImage);
    
    fileMenu->addSeparator();
    
    QAction *saveAction = fileMenu->addAction("Save Stretched Image...");
    saveAction->setShortcut(QKeySequence::Save);
    connect(saveAction, &QAction::triggered, this, &StretchedImageViewer::onSaveStretchedImage);
    
    fileMenu->addSeparator();
    
    QAction *exitAction = fileMenu->addAction("Exit");
    exitAction->setShortcut(QKeySequence::Quit);
    connect(exitAction, &QAction::triggered, this, &QWidget::close);
    
    // View menu
    QMenu *viewMenu = menuBar->addMenu("View");
    
    QAction *resetZoomAction = viewMenu->addAction("Reset Zoom");
    resetZoomAction->setShortcut(QKeySequence("Ctrl+0"));
    connect(resetZoomAction, &QAction::triggered, this, &StretchedImageViewer::onResetZoom);
    
    QAction *autoStretchAction = viewMenu->addAction("Auto Stretch");
    autoStretchAction->setShortcut(QKeySequence("Ctrl+A"));
    connect(autoStretchAction, &QAction::triggered, this, &StretchedImageViewer::onAutoStretch);
}

void StretchedImageViewer::connectSignals()
{
    connect(m_loadButton, &QPushButton::clicked, this, &StretchedImageViewer::onLoadImage);
    connect(m_saveButton, &QPushButton::clicked, this, &StretchedImageViewer::onSaveStretchedImage);
    connect(m_autoStretchButton, &QPushButton::clicked, this, &StretchedImageViewer::onAutoStretch);
    connect(m_resetZoomButton, &QPushButton::clicked, this, &StretchedImageViewer::onResetZoom);
    
    connect(m_stretchTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &StretchedImageViewer::onStretchTypeChanged);
    
    connect(m_gammaSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &StretchedImageViewer::onStretchParametersChanged);
    connect(m_blackPointSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &StretchedImageViewer::onStretchParametersChanged);
    connect(m_whitePointSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &StretchedImageViewer::onStretchParametersChanged);
    
    connect(m_zoomSlider, &QSlider::valueChanged, this, &StretchedImageViewer::onZoomChanged);
}

void StretchedImageViewer::loadStackedImage(const QString &fitsPath)
{
    m_currentImagePath = fitsPath;
    if (loadFITSImage(fitsPath)) {
        calculateImageStatistics();
        applyStretch();
        updateImageDisplay();
        m_saveButton->setEnabled(true);
        statusBar()->showMessage("Image loaded successfully");
    }
}

void StretchedImageViewer::showViewer()
{
    show();
    raise();
    activateWindow();
}

bool StretchedImageViewer::loadFITSImage(const QString &filePath)
{
    fitsfile *fptr = nullptr;
    int status = 0;
    
    m_progressBar->setVisible(true);
    m_progressBar->setValue(0);
    
    QByteArray pathBytes = filePath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        QMessageBox::warning(this, "Error", "Failed to open FITS file");
        m_progressBar->setVisible(false);
        return false;
    }
    
    m_progressBar->setValue(20);
    
    // Get image dimensions
    int naxis;
    long naxes[2];
    if (fits_get_img_param(fptr, 2, nullptr, &naxis, naxes, &status)) {
        fits_close_file(fptr, &status);
        QMessageBox::warning(this, "Error", "Failed to read image parameters");
        m_progressBar->setVisible(false);
        return false;
    }
    
    m_width = naxes[0];
    m_height = naxes[1];
    
    m_progressBar->setValue(40);
    
    // Read image data
    long totalPixels = m_width * m_height;
    m_imageData.resize(totalPixels);
    
    if (fits_read_img(fptr, TFLOAT, 1, totalPixels, nullptr, 
                      m_imageData.data(), nullptr, &status)) {
        fits_close_file(fptr, &status);
        QMessageBox::warning(this, "Error", "Failed to read image data");
        m_progressBar->setVisible(false);
        return false;
    }
    
    m_progressBar->setValue(80);
    
    fits_close_file(fptr, &status);
    
    m_progressBar->setValue(100);
    m_progressBar->setVisible(false);
    
    return true;
}

void StretchedImageViewer::calculateImageStatistics()
{
    if (m_imageData.empty()) return;
    
    auto minMax = std::minmax_element(m_imageData.begin(), m_imageData.end());
    m_minValue = *minMax.first;
    m_maxValue = *minMax.second;
}

void StretchedImageViewer::applyStretch()
{
    if (m_imageData.empty()) return;
    
    std::vector<float> stretchedData;
    QString stretchType = m_stretchTypeCombo->currentText();
    
    if (stretchType == "Linear") {
        stretchedData = applyLinearStretch(m_imageData);
    } else if (stretchType == "Gamma") {
        stretchedData = applyGammaStretch(m_imageData);
    } else if (stretchType == "Logarithmic") {
        stretchedData = applyLogStretch(m_imageData);
    } else if (stretchType == "Sinh") {
        stretchedData = applySinhStretch(m_imageData);
    } else { // Auto
        stretchedData = applyAutoStretch(m_imageData);
    }
    
    m_stretchedImage = convertToQImage(stretchedData, m_width, m_height);
}

std::vector<float> StretchedImageViewer::applyLinearStretch(const std::vector<float> &data)
{
    float blackPoint = percentile(data, m_blackPointSpinBox->value());
    float whitePoint = percentile(data, m_whitePointSpinBox->value());
    
    std::vector<float> result(data.size());
    float range = whitePoint - blackPoint;
    
    for (size_t i = 0; i < data.size(); ++i) {
        float normalized = (data[i] - blackPoint) / range;
        result[i] = std::clamp(normalized, 0.0f, 1.0f);
    }
    
    return result;
}

std::vector<float> StretchedImageViewer::applyGammaStretch(const std::vector<float> &data)
{
    auto linearData = applyLinearStretch(data);
    std::vector<float> result(linearData.size());
    float gamma = m_gammaSpinBox->value();
    
    for (size_t i = 0; i < linearData.size(); ++i) {
        result[i] = std::pow(linearData[i], 1.0f / gamma);
    }
    
    return result;
}

std::vector<float> StretchedImageViewer::applyLogStretch(const std::vector<float> &data)
{
    auto linearData = applyLinearStretch(data);
    std::vector<float> result(linearData.size());
    
    for (size_t i = 0; i < linearData.size(); ++i) {
        result[i] = std::log(1.0f + linearData[i]) / std::log(2.0f);
    }
    
    return result;
}

std::vector<float> StretchedImageViewer::applySinhStretch(const std::vector<float> &data)
{
    auto linearData = applyLinearStretch(data);
    std::vector<float> result(linearData.size());
    
    for (size_t i = 0; i < linearData.size(); ++i) {
        result[i] = std::sinh(linearData[i] * 3.0f) / std::sinh(3.0f);
    }
    
    return result;
}

std::vector<float> StretchedImageViewer::applyAutoStretch(const std::vector<float> &data)
{
    // Auto-calculate optimal black and white points
    float autoBlack = percentile(data, 0.1);
    float autoWhite = percentile(data, 99.5);
    float range = autoWhite - autoBlack;
    
    std::vector<float> result(data.size());
    
    for (size_t i = 0; i < data.size(); ++i) {
        float normalized = (data[i] - autoBlack) / range;
        // Apply mild gamma correction for astronomical images
        result[i] = std::clamp(std::pow(normalized, 0.8f), 0.0f, 1.0f);
    }
    
    return result;
}

float StretchedImageViewer::percentile(const std::vector<float> &data, float percent)
{
    std::vector<float> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());
    
    size_t index = static_cast<size_t>((percent / 100.0) * (sortedData.size() - 1));
    return sortedData[index];
}

QImage StretchedImageViewer::convertToQImage(const std::vector<float> &data, long width, long height)
{
    QImage image(width, height, QImage::Format_Grayscale8);
    
    for (long y = 0; y < height; ++y) {
        uchar *scanLine = image.scanLine(y);
        for (long x = 0; x < width; ++x) {
            float value = data[y * width + x];
            scanLine[x] = static_cast<uchar>(value * 255.0f);
        }
    }
    
    return image;
}

void StretchedImageViewer::updateImageDisplay()
{
    if (m_stretchedImage.isNull()) return;
    
    updateZoom();
}

void StretchedImageViewer::updateZoom()
{
    if (m_stretchedImage.isNull()) return;
    
    m_currentZoom = m_zoomSlider->value() / 100.0;
    
    QSize scaledSize = m_stretchedImage.size() * m_currentZoom;
    QPixmap pixmap = QPixmap::fromImage(m_stretchedImage);
    QPixmap scaledPixmap = pixmap.scaled(scaledSize, Qt::KeepAspectRatio, Qt::SmoothTransformation);
    
    m_imageLabel->setPixmap(scaledPixmap);
    m_zoomLabel->setText(QString("%1%").arg(static_cast<int>(m_currentZoom * 100)));
}

// Slots implementation
void StretchedImageViewer::onLoadImage()
{
    QString fileName = QFileDialog::getOpenFileName(this,
        "Open Stacked FITS Image", "",
        "FITS Files (*.fits *.fit *.fts);;All Files (*)");
    
    if (!fileName.isEmpty()) {
        loadStackedImage(fileName);
    }
}

void StretchedImageViewer::onStretchTypeChanged()
{
    applyStretch();
    updateImageDisplay();
}

void StretchedImageViewer::onStretchParametersChanged()
{
    applyStretch();
    updateImageDisplay();
}

void StretchedImageViewer::onZoomChanged(int value)
{
    updateZoom();
}

void StretchedImageViewer::onResetZoom()
{
    m_zoomSlider->setValue(100);
}

void StretchedImageViewer::onSaveStretchedImage()
{
    if (m_stretchedImage.isNull()) return;
    
    QString fileName = QFileDialog::getSaveFileName(this,
        "Save Stretched Image", "",
        "PNG Files (*.png);;JPEG Files (*.jpg);;TIFF Files (*.tiff)");
    
    if (!fileName.isEmpty()) {
        if (m_stretchedImage.save(fileName)) {
            statusBar()->showMessage("Image saved successfully");
        } else {
            QMessageBox::warning(this, "Error", "Failed to save image");
        }
    }
}

void StretchedImageViewer::onAutoStretch()
{
    m_stretchTypeCombo->setCurrentText("Auto");
    applyStretch();
    updateImageDisplay();
}
