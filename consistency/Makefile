# Makefile for FITS WCS Consistency Checker
# 
# Dependencies:
# - cfitsio: brew install cfitsio (macOS) or apt-get install libcfitsio-dev (Ubuntu)
# - libnova: brew install libnova (macOS) or apt-get install libnova-dev (Ubuntu)

CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
TARGET = fits_wcs_checker
SOURCE = fits_wcs_checker.cpp

# Library paths and flags
LIBS = -lcfitsio -lnova -lm

# Platform-specific settings
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    # macOS with Homebrew
    CXXFLAGS += -I/opt/homebrew/include
    LDFLAGS += -L/opt/homebrew/lib
else ifeq ($(UNAME_S),Linux)
    # Linux
    CXXFLAGS += -I/usr/include
    LDFLAGS += -L/usr/lib
endif

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(TARGET) $(SOURCE) $(LIBS)

clean:
	rm -f $(TARGET)

install-deps-macos:
	brew install cfitsio libnova

install-deps-ubuntu:
	sudo apt-get update
	sudo apt-get install libcfitsio-dev libnova-dev

test: $(TARGET)
	@echo "Testing with sample FITS file..."
	@echo "Usage: ./$(TARGET) <fits_file_path>"

.PHONY: all clean install-deps-macos install-deps-ubuntu test
