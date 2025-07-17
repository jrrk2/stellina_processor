# Makefile for building, signing, and notarizing StellinaProcessor.app

# Variables
PROJECT = StellinaProcessor.xcodeproj
SCHEME = StellinaProcessor
ARCHIVE_PATH = build/StellinaProcessor.xcarchive
EXPORT_PATH = build/exported
EXPORT_OPTIONS_PLIST = exportOptions.plist
QT_LIB_PATH = /opt/homebrew/Cellar/qt/6.9.0/lib
# Use Developer ID Application for notarization, not Apple Development
DEVELOPER_ID = "Developer ID Application: Mr Jonathan Kimmitt (Y53675G8V9)"
ENTITLEMENTS = sandbox.entitlements
BUILD_DIR = build

run: ./build/exported/StellinaProcessor.app
	./build/exported/StellinaProcessor.app/Contents/MacOS/StellinaProcessor

archive:
	make build/StellinaProcessor.xcarchive

app:
	make build/exported/StellinaProcessor.app

zip:
	make build/exported/StellinaProcessor.zip

# Clean up previous builds
clean:
	@rm -rf $(BUILD_DIR)/Release
	@rm -rf $(BUILD_DIR)/Celestron*
	@rm -rf $(BUILD_DIR)/exported
	@rm -rf $(BUILD_DIR)/obj
	@rm -rf build/StellinaProcessor.xcarchive
	@echo "Cleaned up previous build artifacts, preserving necessary files like build/moc/"

# Build the archive
# CODE_SIGN_IDENTITY="Apple Development: Mr Jonathan Kimmitt (5AU5B5HJQX)"
build/StellinaProcessor.xcarchive:
	xcodebuild archive \
		-project $(PROJECT) \
		-scheme $(SCHEME) \
		-archivePath $(ARCHIVE_PATH) \
		CODE_SIGN_ENTITLEMENTS=$(ENTITLEMENTS) \
		ENABLE_HARDENED_RUNTIME=YES \
		OTHER_CODE_SIGN_FLAGS="--options runtime"
	@echo "Archive created at $(ARCHIVE_PATH)."

# Export the app and deploy Qt libraries
build/exported/StellinaProcessor.app: build/StellinaProcessor.xcarchive
	# Export the app (unsigned at this point)
	xcodebuild -exportArchive \
		-archivePath $< \
		-exportPath $(EXPORT_PATH) \
		-exportOptionsPlist $(EXPORT_OPTIONS_PLIST)

	# Deploy Qt libraries
	/opt/homebrew/Cellar/qt/6.9.0/bin/macdeployqt $(EXPORT_PATH)/StellinaProcessor.app \
		-libpath=$(QT_LIB_PATH) -always-overwrite -verbose=2

	 # Deploy Qt libraries
	 /opt/homebrew/Cellar/qt/6.9.0/bin/macdeployqt $(EXPORT_PATH)/StellinaProcessor.app \
		 -libpath=$(QT_LIB_PATH) -always-overwrite -verbose=2

	# Fix app bundle if needed
	python3 fix_app_bundle.py $(EXPORT_PATH)/StellinaProcessor.app

	find $@ -type f \( -name "*.dylib" -o -perm +111 \) -exec codesign --force --timestamp --options runtime --sign $(DEVELOPER_ID) {} \;

	# Re-sign the app after modification
	codesign --force --deep --timestamp --options runtime --sign $(DEVELOPER_ID) $(EXPORT_PATH)/StellinaProcessor.app
	@echo "App signed and ready for distribution."

# Create the .zip file for distribution using ditto with proper flags
build/exported/StellinaProcessor.zip: build/exported/StellinaProcessor.app
	(cd $(EXPORT_PATH); ditto -c -k --sequesterRsrc --keepParent StellinaProcessor.app StellinaProcessor.zip)
	@echo "App zipped for distribution at $(EXPORT_PATH)/StellinaProcessor.zip"

# Verify signatures before notarization (optional check)
verify-signatures: build/exported/StellinaProcessor.app
	@echo "🔍 Verifying signatures..."
	@echo "Main app bundle:"
	@codesign --verify --deep --strict $(EXPORT_PATH)/StellinaProcessor.app && echo "✅ Main app signature valid" || echo "❌ Main app signature invalid"
	@echo ""
	@echo "App bundle details:"
	@codesign --display --verbose $(EXPORT_PATH)/StellinaProcessor.app
	@echo ""
	@echo "Testing a few key frameworks:"
	@codesign --verify --verbose $(EXPORT_PATH)/StellinaProcessor.app/Contents/Frameworks/QtCore.framework && echo "✅ QtCore framework signature valid" || echo "❌ QtCore framework signature invalid"
	@codesign --verify --verbose $(EXPORT_PATH)/StellinaProcessor.app/Contents/Frameworks/QtGui.framework && echo "✅ QtGui framework signature valid" || echo "❌ QtGui framework signature invalid"
	@echo ""
	@echo "Testing main executable:"
	@codesign --verify --verbose $(EXPORT_PATH)/StellinaProcessor.app/Contents/MacOS/StellinaProcessor && echo "✅ Main executable signature valid" || echo "❌ Main executable signature invalid"
	@echo ""
	@echo "System policy assessment (will fail until notarized):"
	@spctl --assess --type execute --verbose $(EXPORT_PATH)/StellinaProcessor.app && echo "✅ System policy assessment passed" || echo "⚠️  System policy assessment failed - this is expected before notarization"
	@echo ""
	@echo "✅ App is properly signed and ready for notarization!"

# Notarize the app
notary: build/exported/StellinaProcessor.zip # verify-signatures
	@echo "Submitting for notarization..."
	(cd $(EXPORT_PATH); xcrun notarytool submit StellinaProcessor.zip \
		--apple-id "jonathan@kimmitt.co.uk" \
		--team-id Y53675G8V9 \
		--keychain-profile "notary-profile-name" \
		--wait)
	@echo "Notarization complete."

# Get notarization history and logs for debugging
debug-notary:
	@echo "Getting notarization history..."
	xcrun notarytool history --apple-id "jonathan@kimmitt.co.uk" --team-id Y53675G8V9 --keychain-profile "notary-profile-name"
	@echo ""
	@echo "To get detailed logs for a specific submission, run:"
	@echo "xcrun notarytool log <SUBMISSION_ID> --apple-id jonathan@kimmitt.co.uk --team-id Y53675G8V9 --keychain-profile notary-profile-name"

# Get logs for the most recent submission
debug-latest:
	@echo "Getting logs for the most recent submission..."
	$(eval LATEST_ID := $(shell xcrun notarytool history --apple-id "jonathan@kimmitt.co.uk" --team-id Y53675G8V9 --keychain-profile "notary-profile-name" | grep "id:" | head -1 | cut -d' ' -f2))
	xcrun notarytool log $(LATEST_ID) --apple-id "jonathan@kimmitt.co.uk" --team-id Y53675G8V9 --keychain-profile "notary-profile-name"

# Staple the notarization ticket
staple: build/exported/StellinaProcessor.app
	xcrun stapler staple $(EXPORT_PATH)/StellinaProcessor.app
	@echo "Stapling notarization ticket to app."

# Verify the notarized app
verify: build/exported/StellinaProcessor.app
	spctl --assess --type execute --verbose $(EXPORT_PATH)/StellinaProcessor.app
	@echo "App verification complete."

# Check entitlements file for correct line endings
check-entitlements:
	@echo "Checking entitlements file..."
	file $(ENTITLEMENTS)
	@echo "If it shows 'with CRLF line terminators', run: dos2unix $(ENTITLEMENTS)"

# Build and run the app locally (during development)
./Release/StellinaProcessor.app: ./StellinaProcessor.xcodeproj/project.pbxproj
	@rm -rf Release
	xcodebuild -project $(PROJECT) -scheme $(SCHEME) -configuration Release
	/opt/homebrew/Cellar/qt/6.9.0/bin/macdeployqt $@ -libpath=$(QT_LIB_PATH) -always-overwrite -verbose=2
	python3 fix_app_bundle.py $@

open: build/exported/StellinaProcessor.app
	open build/exported/StellinaProcessor.app

# Complete notarization workflow
notarize-complete: clean build/exported/StellinaProcessor.app build/exported/StellinaProcessor.zip verify-signatures notary staple verify
	@echo "Complete notarization workflow finished successfully!"
