#include "SimplifiedXISFWriter.h"

// Now we can safely include PCL headers in the .cpp file
#include <pcl/Image.h>
#include <pcl/Property.h>
#include <pcl/Variant.h>
#include <pcl/String.h>
#include <pcl/XISF.h>

// Your project includes
#include "StellinaProcessor.h"
#include <structuredefinitions.h>

// Qt includes
#include <QDebug>
#include <QDateTime>

// Private implementation class
class SimplifiedXISFWriterPrivate
{
public:
    struct PropertyData {
        QString id;
        QString type;
        QVariant value;
        QString comment;
        bool isImageProperty;
        
        PropertyData(const QString& i, const QString& t, const QVariant& v, 
                    const QString& c, bool img = false)
            : id(i), type(t), value(v), comment(c), isImageProperty(img) {}
    };

    struct ImageData {
        QString id;
        QByteArray pixelData; // Store as raw bytes to avoid PCL in header
        int width;
        int height;
        int channels;
        QList<PropertyData> properties;
        
        ImageData(const QString& i, const float* pixels, int w, int h, int c)
            : id(i), width(w), height(h), channels(c)
        {
            // Convert float array to byte array
            qint64 pixelCount = static_cast<qint64>(w) * h * c;
            pixelData.resize(pixelCount * sizeof(float));
            memcpy(pixelData.data(), pixels, pixelData.size());
        }
    };

    QString filePath;
    QString creatorApplication = "StellinaProcessor";
    QString lastError;
    
    CompressionType compressionType = CompressionType::ZLib;
    int compressionLevel = 6;
    ChecksumType checksumType = ChecksumType::SHA256;
    int verbosity = 1;
    
    QList<ImageData> images;
    QList<PropertyData> globalProperties;
    QList<PropertyData> pendingImageProperties;

    // Convert our enum to PCL enum
    pcl::XISFCompression::value_type getPCLCompression() const
    {
        switch (compressionType) {
	case CompressionType::None: return (pcl::XISFCompression::value_type) 0;
	case CompressionType::ZLib: return (pcl::XISFCompression::value_type) 1;
	case CompressionType::LZ4:  return (pcl::XISFCompression::value_type) 2;
	case CompressionType::ZSTD: return (pcl::XISFCompression::value_type) 4;
	default: return (pcl::XISFCompression::value_type) 1;
        }
    }

    pcl::XISFChecksum::value_type getPCLChecksum() const
    {
        switch (checksumType) {
            case ChecksumType::None:   return pcl::XISFChecksum::None;
            case ChecksumType::SHA1: return pcl::XISFChecksum::SHA1;
            case ChecksumType::SHA256: return pcl::XISFChecksum::SHA256;
            case ChecksumType::SHA512: return pcl::XISFChecksum::SHA512;
            default: return pcl::XISFChecksum::SHA256;
        }
    }

    pcl::Variant convertToPCLVariant(const QString& type, const QVariant& value)
    {
        if (type == "String") {
            return pcl::Variant(pcl::String(value.toString().toUtf8().constData()));
        }
        else if (type == "Float64") {
            return pcl::Variant(value.toDouble());
        }
        else if (type == "Float32") {
            return pcl::Variant(value.toFloat());
        }
        else if (type == "Int32") {
            return pcl::Variant(value.toInt());
        }
        else if (type == "UInt16") {
            return pcl::Variant(static_cast<uint16_t>(value.toUInt()));
        }
        else if (type == "Boolean") {
            return pcl::Variant(value.toBool());
        }
        else {
            // Default to string
            return pcl::Variant(pcl::String(value.toString().toUtf8().constData()));
        }
    }
};

// Implementation of public interface
SimplifiedXISFWriter::SimplifiedXISFWriter(const QString& filePath, CompressionType compression)
    : d(std::make_unique<SimplifiedXISFWriterPrivate>())
{
    d->filePath = filePath;
    d->compressionType = compression;
    
    qDebug() << "SimplifiedXISFWriter created for" << filePath 
             << "with compression" << static_cast<int>(compression);
}

SimplifiedXISFWriter::~SimplifiedXISFWriter() = default;

bool SimplifiedXISFWriter::addImage(const QString& id, const float* pixels, int width, int height, int channels)
{
    if (!pixels || width <= 0 || height <= 0 || channels <= 0) {
        d->lastError = QString("Invalid image parameters for %1").arg(id);
        qDebug() << d->lastError;
        return false;
    }

    try {
        // Create image data and copy the pending properties
        SimplifiedXISFWriterPrivate::ImageData imageData(id, pixels, width, height, channels);
        imageData.properties = d->pendingImageProperties;
        
        d->images.append(imageData);
        d->pendingImageProperties.clear();
        
        qDebug() << "Added image" << id << "(" << width << "x" << height << "x" << channels << ")";
        return true;
    }
    catch (const std::exception& e) {
        d->lastError = QString("Error adding image %1: %2").arg(id, e.what());
        qDebug() << d->lastError;
        return false;
    }
}

void SimplifiedXISFWriter::addProperty(const QString& id, const QString& type, const QVariant& value, const QString& comment)
{
    if (id.isEmpty() || type.isEmpty()) {
        qDebug() << "Invalid property: id or type is empty";
        return;
    }
    
    d->globalProperties.append(SimplifiedXISFWriterPrivate::PropertyData(id, type, value, comment, false));
    qDebug() << "Added global property" << id << "type" << type;
}

void SimplifiedXISFWriter::addImageProperty(const QString& id, const QString& type, const QVariant& value, const QString& comment)
{
    if (id.isEmpty() || type.isEmpty()) {
        qDebug() << "Invalid image property: id or type is empty";
        return;
    }
    
    d->pendingImageProperties.append(SimplifiedXISFWriterPrivate::PropertyData(id, type, value, comment, true));
    qDebug() << "Added image property" << id << "type" << type;
}

void SimplifiedXISFWriter::addWCSData(const FITSImage::Solution& solution)
{
    addImageProperty("WCS:CRVAL1", "Float64", solution.ra, "Reference RA (degrees)");
    addImageProperty("WCS:CRVAL2", "Float64", solution.dec, "Reference Dec (degrees)");
    addImageProperty("WCS:PIXSCALE", "Float64", solution.pixscale, "Pixel scale (arcsec/pixel)");
    addImageProperty("WCS:ORIENTAT", "Float64", solution.orientation, "Position angle (degrees)");
    addImageProperty("WCS:FIELDW", "Float64", solution.fieldWidth, "Field width (arcmin)");
    addImageProperty("WCS:FIELDH", "Float64", solution.fieldHeight, "Field height (arcmin)");
    
    QString parity = (solution.parity == FITSImage::NEGATIVE) ? "negative" : "positive";
    addImageProperty("WCS:PARITY", "String", parity, "Image parity");
    
    if (solution.raError > 0) {
        addImageProperty("WCS:RAERR", "Float64", solution.raError, "RA error (arcsec)");
        addImageProperty("WCS:DECERR", "Float64", solution.decError, "Dec error (arcsec)");
    }
    
    qDebug() << "Added WCS data: RA" << solution.ra << "Dec" << solution.dec;
}

void SimplifiedXISFWriter::addStellinaMetadata(const StellinaImageData& imageData)
{
    addImageProperty("Stellina:Altitude", "Float64", imageData.altitude, "Mount altitude (degrees)");
    addImageProperty("Stellina:Azimuth", "Float64", imageData.azimuth, "Mount azimuth (degrees)");
    addImageProperty("Stellina:Exposure", "Int32", imageData.exposureSeconds, "Exposure time (seconds)");
    addImageProperty("Stellina:Temperature", "Int32", imageData.temperatureKelvin, "Sensor temperature (K)");
    addImageProperty("Stellina:Binning", "String", imageData.binning, "Binning mode");
    addImageProperty("Stellina:BayerPattern", "String", imageData.bayerPattern, "CFA pattern");
    
    if (imageData.hasCalculatedCoords) {
        addImageProperty("Stellina:CalculatedRA", "Float64", imageData.calculatedRA, "Calculated RA (degrees)");
        addImageProperty("Stellina:CalculatedDec", "Float64", imageData.calculatedDec, "Calculated Dec (degrees)");
    }
    
    qDebug() << "Added Stellina metadata";
}

void SimplifiedXISFWriter::addProcessingHistory(const QString& step, const QJsonObject& parameters)
{
    QString timestamp = QDateTime::currentDateTimeUtc().toString(Qt::ISODate);
    addProperty(QString("Processing:%1:Timestamp").arg(step), "String", timestamp, "Processing timestamp");
    
    for (auto it = parameters.begin(); it != parameters.end(); ++it) {
        QString propId = QString("Processing:%1:%2").arg(step, it.key());
        addProperty(propId, "String", it.value().toString(), "Processing parameter");
    }
    
    qDebug() << "Added processing history for step:" << step;
}

void SimplifiedXISFWriter::setCompression(CompressionType compression, int level)
{
    d->compressionType = compression;
    d->compressionLevel = level;
    qDebug() << "Set compression to" << static_cast<int>(compression) << "level" << level;
}

void SimplifiedXISFWriter::setChecksum(ChecksumType checksum)
{
    d->checksumType = checksum;
    qDebug() << "Set checksum to" << static_cast<int>(checksum);
}

void SimplifiedXISFWriter::setCreatorApplication(const QString& appName)
{
    d->creatorApplication = appName;
}

void SimplifiedXISFWriter::setVerbosity(int level)
{
    d->verbosity = level;
}
bool SimplifiedXISFWriter::write()
{
    if (d->images.isEmpty()) {
        d->lastError = "No images to write";
        qDebug() << d->lastError;
        return false;
    }

    try {
        // Create PCL XISF writer
        pcl::XISFWriter writer;
        
        // Set up XISF options FIRST (before Create)
        pcl::XISFOptions xisfOptions;
        xisfOptions.compressionCodec = d->getPCLCompression();
        xisfOptions.compressionLevel = d->compressionLevel;
        xisfOptions.checksumAlgorithm = d->getPCLChecksum();
        xisfOptions.verbosity = d->verbosity;
        writer.SetOptions(xisfOptions);
        
        // Set up image options BEFORE Create
        pcl::ImageOptions imageOptions;
        imageOptions.bitsPerSample = 32; // 32-bit float
        imageOptions.ieeefpSampleFormat = true;
        imageOptions.embedProperties = true;
        imageOptions.embedRGBWS = true;
        
        // Set creator application BEFORE Create
        writer.SetCreatorApplication(pcl::String(d->creatorApplication.toUtf8().constData()));
        
        // NOW create the file (after all options are set)
        pcl::String pclPath(d->filePath.toUtf8().constData());
        writer.Create(pclPath, d->images.size());
        writer.SetImageOptions(imageOptions);
        
        // Add global properties
        for (const auto& prop : d->globalProperties) {
            pcl::Variant pclValue = d->convertToPCLVariant(prop.type, prop.value);
            writer.WriteProperty(pcl::IsoString(prop.id.toUtf8().constData()), pclValue);
        }
        
        // Write each image
        for (const auto& imageData : d->images) {
            // Set image ID
            writer.SetImageId(pcl::IsoString(imageData.id.toUtf8().constData()));
            
            // Add image-specific properties
            for (const auto& prop : imageData.properties) {
                pcl::Variant pclValue = d->convertToPCLVariant(prop.type, prop.value);
                writer.WriteImageProperty(pcl::IsoString(prop.id.toUtf8().constData()), pclValue);
            }
            
            // Create PCL Image from stored pixel data
            pcl::Image pclImage(imageData.width, imageData.height, 
                               imageData.channels == 1 ? pcl::ColorSpace::Gray : pcl::ColorSpace::RGB);
            
            const float* pixels = reinterpret_cast<const float*>(imageData.pixelData.constData());
            const float* src = pixels;
            
            for (int c = 0; c < imageData.channels; ++c) {
                for (int y = 0; y < imageData.height; ++y) {
                    for (int x = 0; x < imageData.width; ++x) {
                        pclImage.Pixel(x, y, c) = *src++;
                    }
                }
            }
            
            // Write the image
            writer.WriteImage(pclImage);
            
            qDebug() << "Wrote image" << imageData.id;
        }
        
        // Close the file
        writer.Close();
        
        qDebug() << "Successfully wrote PCL XISF file:" << d->filePath;
        return true;
    }
    catch (const pcl::Error& e) {
        d->lastError = QString("PCL Error: %1").arg(e.Message().c_str());
        qDebug() << d->lastError;
        return false;
    }
    catch (const std::exception& e) {
        d->lastError = QString("Standard error: %1").arg(e.what());
        qDebug() << d->lastError;
        return false;
    }
}

void SimplifiedXISFWriter::clear()
{
    d->images.clear();
    d->globalProperties.clear();
    d->pendingImageProperties.clear();
    d->lastError.clear();
    qDebug() << "Cleared all images and properties";
}

QString SimplifiedXISFWriter::lastError() const
{
    return d->lastError;
}

bool SimplifiedXISFWriter::hasImages() const
{
    return !d->images.isEmpty();
}

int SimplifiedXISFWriter::imageCount() const
{
    return d->images.size();
}
