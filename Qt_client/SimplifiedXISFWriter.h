// =============================================================================
// SIMPLIFIED XISF IMPLEMENTATION
// =============================================================================

#include <QXmlStreamWriter>
#include <QFile>
#include <QFileInfo>
#include <QDateTime>
#include <QCryptographicHash>
#include <QtEndian>
#include <QBuffer>
#include <QDebug>
#include <QJsonObject>
#include "struct.h"

// Simplified XISF writer that doesn't need PCL
class SimplifiedXISFWriter {
public:
    struct ImageBlock {
        QString id;
        int width;
        int height;
        int channels;
        QString sampleFormat;  // "Float32", "UInt16", etc.
        QString colorSpace;    // "Grayscale", "RGB"
        QByteArray data;       // Raw pixel data
        qint64 location;       // Offset in file
        qint64 size;           // Size in bytes
        QString compression;   // "none", "lz4", "zlib"
        
        ImageBlock() : width(0), height(0), channels(1), 
                      sampleFormat("Float32"), colorSpace("Grayscale"), 
                      location(0), size(0), compression("none") {}
    };
    
    struct Property {
        QString id;
        QString type;      // "String", "Float64", "Int32", "Boolean"
        QVariant value;
        QString comment;
        
        Property() = default;
        Property(const QString& _id, const QString& _type, const QVariant& _value, const QString& _comment = QString())
            : id(_id), type(_type), value(_value), comment(_comment) {}
    };

private:
    QString m_filePath;
    QList<ImageBlock> m_images;
    QList<Property> m_properties;
    QString m_creatorApplication;
    QString m_creationTime;
    
public:
    SimplifiedXISFWriter(const QString& filePath) 
        : m_filePath(filePath)
        , m_creatorApplication("StellinaProcessor")
        , m_creationTime(QDateTime::currentDateTimeUtc().toString(Qt::ISODate)) {
    }
    
    // Add image data (takes ownership of data)
    bool addImage(const QString& id, const float* pixels, int width, int height, int channels = 1) {
        ImageBlock image;
        image.id = id;
        image.width = width;
        image.height = height;
        image.channels = channels;
        image.sampleFormat = "Float32";
        image.colorSpace = (channels == 1) ? "Grayscale" : "RGB";
        
        // Convert to little-endian byte array
        qint64 pixelCount = width * height * channels;
        image.size = pixelCount * sizeof(float);
        image.data.resize(image.size);
        
        const float* src = pixels;
        float* dst = reinterpret_cast<float*>(image.data.data());
        
        for (qint64 i = 0; i < pixelCount; ++i) {
            dst[i] = qToLittleEndian(src[i]);
        }
        
        m_images.append(image);
        return true;
    }
    
    // Add property
    void addProperty(const QString& id, const QString& type, const QVariant& value, const QString& comment = QString()) {
        m_properties.append(Property(id, type, value, comment));
    }
    
    // Convenience methods for common metadata
    void addWCSData(const FITSImage::Solution& solution) {
        addProperty("WCS:CRVAL1", "Float64", solution.ra, "Reference RA (degrees)");
        addProperty("WCS:CRVAL2", "Float64", solution.dec, "Reference Dec (degrees)");
        addProperty("WCS:PIXSCALE", "Float64", solution.pixscale, "Pixel scale (arcsec/pixel)");
        addProperty("WCS:ORIENTAT", "Float64", solution.orientation, "Position angle (degrees)");
        addProperty("WCS:FIELDW", "Float64", solution.fieldWidth, "Field width (arcmin)");
        addProperty("WCS:FIELDH", "Float64", solution.fieldHeight, "Field height (arcmin)");
        
        QString parity = (solution.parity == FITSImage::NEGATIVE) ? "negative" : "positive";
        addProperty("WCS:PARITY", "String", parity, "Image parity");
        
        if (solution.raError > 0) {
            addProperty("WCS:RAERR", "Float64", solution.raError, "RA error (arcsec)");
            addProperty("WCS:DECERR", "Float64", solution.decError, "Dec error (arcsec)");
        }
    }
    
    void addStellinaMetadata(const StellinaImageData& imageData) {
        addProperty("Stellina:Altitude", "Float64", imageData.altitude, "Mount altitude (degrees)");
        addProperty("Stellina:Azimuth", "Float64", imageData.azimuth, "Mount azimuth (degrees)");
        addProperty("Stellina:Exposure", "Int32", imageData.exposureSeconds, "Exposure time (seconds)");
        addProperty("Stellina:Temperature", "Int32", imageData.temperatureKelvin, "Sensor temperature (K)");
        addProperty("Stellina:Binning", "String", imageData.binning, "Binning mode");
        addProperty("Stellina:BayerPattern", "String", imageData.bayerPattern, "CFA pattern");
        
        if (imageData.hasCalculatedCoords) {
            addProperty("Stellina:CalculatedRA", "Float64", imageData.calculatedRA, "Calculated RA (degrees)");
            addProperty("Stellina:CalculatedDec", "Float64", imageData.calculatedDec, "Calculated Dec (degrees)");
        }
    }
    
    void addProcessingHistory(const QString& step, const QJsonObject& parameters) {
        QString timestamp = QDateTime::currentDateTimeUtc().toString(Qt::ISODate);
        addProperty(QString("Processing:%1:Timestamp").arg(step), "String", timestamp, "Processing timestamp");
        
        for (auto it = parameters.begin(); it != parameters.end(); ++it) {
            QString propId = QString("Processing:%1:%2").arg(step, it.key());
            addProperty(propId, "String", it.value().toString(), "Processing parameter");
        }
    }

        // REPLACE your existing write() method with this version:
    bool write() {
        QFile file(m_filePath);
        if (!file.open(QIODevice::WriteOnly)) {
            qDebug() << "Failed to open file for writing:" << m_filePath;
            return false;
        }
        
        // Step 1: Generate initial XML to calculate size
        QByteArray xmlData;
        if (!writeXMLHeader(xmlData)) {
            return false;
        }
        
        // Step 2: Calculate block alignment (required for monolithic XISF)
        const qint64 blockAlignment = 4096;
        qint64 headerSize = xmlData.size();
        
        // Calculate padded header size (must be aligned to 4096 bytes)
        qint64 paddedHeaderSize = ((headerSize + blockAlignment - 1) / blockAlignment) * blockAlignment;
        
        // Step 3: Calculate image block locations after padded header
        qint64 currentOffset = paddedHeaderSize;
        for (auto& image : m_images) {
            image.location = currentOffset;
            currentOffset += image.size;
            // Align next block (though only one image for now)
            currentOffset = ((currentOffset + blockAlignment - 1) / blockAlignment) * blockAlignment;
        }
        
        // Step 4: Regenerate XML with correct locations
        xmlData.clear();
        if (!writeXMLHeader(xmlData)) {
            return false;
        }
        
        // Step 5: Pad header to alignment boundary
        headerSize = xmlData.size();
        paddedHeaderSize = ((headerSize + blockAlignment - 1) / blockAlignment) * blockAlignment;
        QByteArray padding(paddedHeaderSize - headerSize, 0);  // Zero-filled padding
        xmlData.append(padding);
        
        // Step 6: Write padded header
        if (file.write(xmlData) != xmlData.size()) {
            qDebug() << "Failed to write XML header";
            return false;
        }
        
        // Step 7: Write image data blocks at correct aligned positions
        for (const auto& image : m_images) {
            // Verify we're at the expected position
            if (file.pos() != image.location) {
                qDebug() << "Position mismatch: expected" << image.location << "but at" << file.pos();
                if (!file.seek(image.location)) {
                    qDebug() << "Failed to seek to image location" << image.location;
                    return false;
                }
            }
            
            if (file.write(image.data) != image.data.size()) {
                qDebug() << "Failed to write image data for" << image.id;
                return false;
            }
        }
        
        qDebug() << "Successfully wrote XISF file:" << m_filePath;
        qDebug() << "Header size:" << headerSize << "Padded to:" << paddedHeaderSize;
        qDebug() << "Image location:" << (m_images.isEmpty() ? 0 : m_images.first().location);
        
        return true;
    }
    
private:
    bool writeXMLHeader(QByteArray& xmlData) {
        QBuffer buffer(&xmlData);
        buffer.open(QIODevice::WriteOnly);
        
        QXmlStreamWriter xml(&buffer);
        xml.setAutoFormatting(true);
        xml.setAutoFormattingIndent(2);
        
        xml.writeStartDocument("1.0");
        
        // XISF root element - FIXED namespace handling
        xml.writeStartElement("xisf");
        xml.writeAttribute("version", "1.0");
        
        // Write namespace attributes explicitly (QXmlStreamWriter sometimes drops them)
        xml.writeAttribute("xmlns", "http://www.pixinsight.com/xisf");
        xml.writeAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
        xml.writeAttribute("xsi:schemaLocation", "http://www.pixinsight.com/xisf http://pixinsight.com/xisf/xisf-1.0.xsd");
        
        // Metadata
        xml.writeStartElement("Metadata");
        
        // Creator information
        xml.writeStartElement("Property");
        xml.writeAttribute("id", "XISF:CreatorApplication");
        xml.writeAttribute("type", "String");
        xml.writeAttribute("value", m_creatorApplication);
        xml.writeEndElement(); // Property
        
        xml.writeStartElement("Property");
        xml.writeAttribute("id", "XISF:CreationTime");
        xml.writeAttribute("type", "String");
        xml.writeAttribute("value", m_creationTime);
        xml.writeEndElement(); // Property
        
        // IMPORTANT: Add the required BlockAlignmentSize property
        xml.writeStartElement("Property");
        xml.writeAttribute("id", "XISF:BlockAlignmentSize");
        xml.writeAttribute("type", "UInt16");
        xml.writeAttribute("value", "4096");
        xml.writeAttribute("comment", "Block alignment size in bytes");
        xml.writeEndElement(); // Property
        
        // Write all custom properties with better formatting
        for (const auto& prop : m_properties) {
            xml.writeStartElement("Property");
            xml.writeAttribute("id", prop.id);
            xml.writeAttribute("type", prop.type);
            
            // Better numeric formatting
            QString valueStr;
            if (prop.type == "Float64") {
                valueStr = QString::number(prop.value.toDouble(), 'g', 15);
            } else if (prop.type == "Int32" || prop.type == "UInt16") {
                valueStr = QString::number(prop.value.toInt());
            } else if (prop.type == "Boolean") {
                valueStr = prop.value.toBool() ? "true" : "false";
            } else {
                valueStr = prop.value.toString();
            }
            
            xml.writeAttribute("value", valueStr);
            if (!prop.comment.isEmpty()) {
                xml.writeAttribute("comment", prop.comment);
            }
            xml.writeEndElement(); // Property
        }
        
        xml.writeEndElement(); // Metadata
        
        // Image elements
        for (const auto& image : m_images) {
            xml.writeStartElement("Image");
            xml.writeAttribute("id", image.id);
            xml.writeAttribute("geometry", QString("%1:%2:%3").arg(image.width).arg(image.height).arg(image.channels));
            xml.writeAttribute("sampleFormat", image.sampleFormat);
            xml.writeAttribute("colorSpace", image.colorSpace);
            xml.writeAttribute("location", QString("attachment:%1:%2").arg(image.location).arg(image.size));
            xml.writeEndElement(); // Image
        }
        
        xml.writeEndElement(); // xisf
        xml.writeEndDocument();
        
        // Debug: Let's see what XML we actually generated
        qDebug() << "Generated XML (first 500 chars):";
        qDebug() << xmlData.left(500);
        
        return !xml.hasError();
    }
};
