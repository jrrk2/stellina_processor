#!/usr/bin/env python3
"""
FITS WCS Consistency Checker

This script reads a plate-solved FITS file and checks the consistency
of World Coordinate System (WCS) keywords, particularly focusing on
astrometry.net solved files.

Requirements:
- astropy
- numpy

Usage:
    python wcs_checker.py <fits_file_path>
"""

import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import warnings

def check_wcs_consistency(fits_file):
    """
    Check WCS consistency in a FITS file.
    
    Parameters:
    -----------
    fits_file : str
        Path to the FITS file
    """
    
    print(f"Analyzing FITS file: {fits_file}")
    print("=" * 60)
    
    try:
        # Open FITS file
        with fits.open(fits_file) as hdul:
            header = hdul[0].header
            data = hdul[0].data
            
            print(f"Image dimensions: {data.shape}")
            print(f"BITPIX: {header.get('BITPIX', 'N/A')}")
            print(f"Instrument: {header.get('INSTRUME', 'N/A')}")
            print(f"Date-Obs: {header.get('DATE-OBS', 'N/A')}")
            print()
            
            # Check if WCS exists
            wcs_keywords = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2']
            has_wcs = all(keyword in header for keyword in wcs_keywords)
            
            if not has_wcs:
                print("❌ No WCS information found in header")
                return False
            
            print("✅ WCS keywords found")
            
            # Create WCS object
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                wcs = WCS(header)
            
            # Basic WCS info
            print(f"WCS Type: {header['CTYPE1']}, {header['CTYPE2']}")
            print(f"Reference point (CRVAL): RA={header['CRVAL1']:.6f}°, Dec={header['CRVAL2']:.6f}°")
            print(f"Reference pixel (CRPIX): X={header['CRPIX1']:.2f}, Y={header['CRPIX2']:.2f}")
            print()
            
            # Check coordinate transformation matrix
            check_cd_matrix(header)
            
            # Check pixel scale consistency
            check_pixel_scale(header, wcs)
            
            # Check field center consistency
            check_field_center(header, wcs, data.shape)
            
            # Check SIP distortion if present
            check_sip_distortion(header)
            
            # Check coordinate system consistency (Alt/Az vs RA/Dec)
            check_coordinate_system_consistency(header)
            
            # Check astrometry.net specific fields
            check_astrometry_net_fields(header)
            
            # Validate coordinate transformation
            validate_coordinate_transform(wcs, data.shape)
            
            return True
            
    except Exception as e:
        print(f"❌ Error reading FITS file: {e}")
        return False

def check_cd_matrix(header):
    """Check CD matrix consistency"""
    print("CD Matrix Analysis:")
    print("-" * 20)
    
    cd_keys = ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
    cd_matrix = np.zeros((2, 2))
    
    missing_keys = []
    for i, key in enumerate(cd_keys):
        if key in header:
            cd_matrix[i//2, i%2] = header[key]
        else:
            missing_keys.append(key)
    
    if missing_keys:
        print(f"❌ Missing CD matrix elements: {missing_keys}")
        return
    
    print(f"CD1_1: {cd_matrix[0,0]:.2e}, CD1_2: {cd_matrix[0,1]:.2e}")
    print(f"CD2_1: {cd_matrix[1,0]:.2e}, CD2_2: {cd_matrix[1,1]:.2e}")
    
    # Calculate determinant
    det = np.linalg.det(cd_matrix)
    print(f"Determinant: {det:.2e}")
    
    if abs(det) < 1e-12:
        print("❌ CD matrix is singular (determinant ≈ 0)")
    else:
        print("✅ CD matrix is non-singular")
    
    # Calculate pixel scales
    scale_x = np.sqrt(cd_matrix[0,0]**2 + cd_matrix[1,0]**2) * 3600  # arcsec/pixel
    scale_y = np.sqrt(cd_matrix[0,1]**2 + cd_matrix[1,1]**2) * 3600  # arcsec/pixel
    
    print(f"Pixel scale X: {scale_x:.3f} arcsec/pixel")
    print(f"Pixel scale Y: {scale_y:.3f} arcsec/pixel")
    
    if abs(scale_x - scale_y) / max(scale_x, scale_y) > 0.01:
        print("⚠️  Significant difference in X/Y pixel scales (>1%)")
    else:
        print("✅ X/Y pixel scales are consistent")
    
    print()

def check_pixel_scale(header, wcs):
    """Check pixel scale consistency"""
    print("Pixel Scale Analysis:")
    print("-" * 20)
    
    # Get pixel scale from WCS
    try:
        pixel_scales = wcs.proj_plane_pixel_scales()
        avg_scale = np.mean(pixel_scales) * 3600  # Convert to arcsec/pixel
        print(f"WCS average pixel scale: {avg_scale:.3f} arcsec/pixel")
        
        # Compare with physical parameters if available
        if 'PIXSZ' in header and 'FOCAL' in header:
            pixsz = header['PIXSZ']  # microns
            focal = header['FOCAL']  # mm
            theoretical_scale = (pixsz / 1000) / focal * 206265  # arcsec/pixel
            
            print(f"Theoretical pixel scale: {theoretical_scale:.3f} arcsec/pixel")
            print(f"  (from PIXSZ={pixsz}μm, FOCAL={focal}mm)")
            
            diff_percent = abs(avg_scale - theoretical_scale) / theoretical_scale * 100
            print(f"Difference: {diff_percent:.1f}%")
            
            if diff_percent > 5:
                print("⚠️  Large difference between WCS and theoretical pixel scale")
            else:
                print("✅ WCS and theoretical pixel scales agree")
        
    except Exception as e:
        print(f"❌ Error calculating pixel scale: {e}")
    
    print()

def check_field_center(header, wcs, image_shape):
    """Check field center consistency"""
    print("Field Center Analysis:")
    print("-" * 20)
    
    # Get image center
    center_x = image_shape[1] / 2
    center_y = image_shape[0] / 2
    
    try:
        # Convert image center to world coordinates
        center_ra, center_dec = wcs.pixel_to_world_values(center_x, center_y)
        
        print(f"Image center pixel: ({center_x:.1f}, {center_y:.1f})")
        print(f"Image center RA/Dec: ({center_ra:.6f}°, {center_dec:.6f}°)")
        
        # Convert to HMS/DMS format
        coord = SkyCoord(center_ra*u.deg, center_dec*u.deg)
        print(f"Image center: {coord.to_string('hmsdms')}")
        
        # Compare with reference point
        crval1, crval2 = header['CRVAL1'], header['CRVAL2']
        crpix1, crpix2 = header['CRPIX1'], header['CRPIX2']
        
        print(f"Reference point (CRVAL): ({crval1:.6f}°, {crval2:.6f}°)")
        print(f"Reference pixel (CRPIX): ({crpix1:.2f}, {crpix2:.2f})")
        
        # Calculate separation
        ref_coord = SkyCoord(crval1*u.deg, crval2*u.deg)
        separation = coord.separation(ref_coord).arcsec
        
        print(f"Separation from reference: {separation:.1f} arcsec")
        
        # Check if reference pixel is reasonable
        if (0 <= crpix1 <= image_shape[1] and 0 <= crpix2 <= image_shape[0]):
            print("✅ Reference pixel is within image bounds")
        else:
            print("⚠️  Reference pixel is outside image bounds")
        
    except Exception as e:
        print(f"❌ Error calculating field center: {e}")
    
    print()

def check_sip_distortion(header):
    """Check SIP distortion polynomials"""
    print("SIP Distortion Analysis:")
    print("-" * 20)
    
    # Check for SIP keywords
    sip_keywords = ['A_ORDER', 'B_ORDER', 'AP_ORDER', 'BP_ORDER']
    has_sip = any(key in header for key in sip_keywords)
    
    if not has_sip:
        print("No SIP distortion coefficients found")
        print()
        return
    
    print("✅ SIP distortion coefficients found")
    
    for order_key in sip_keywords:
        if order_key in header:
            order = header[order_key]
            prefix = order_key.split('_')[0]
            print(f"{prefix} polynomial order: {order}")
            
            # Count coefficients
            coeff_count = 0
            max_coeff = 0
            for key in header:
                if key.startswith(f"{prefix}_") and key != order_key:
                    coeff_count += 1
                    try:
                        value = abs(header[key])
                        max_coeff = max(max_coeff, value)
                    except:
                        pass
            
            print(f"  Coefficients found: {coeff_count}")
            print(f"  Maximum coefficient: {max_coeff:.2e}")
    
    print()

def check_coordinate_system_consistency(header):
    """Check consistency between Alt/Az mount coordinates and RA/Dec WCS solution"""
    print("Coordinate System Consistency Analysis:")
    print("-" * 40)
    
    # Required keywords for comparison
    required_keys = ['STELLALT', 'STELLAZ', 'OBSLAT', 'OBSLONG', 'JULIAN', 
                     'LST', 'HOURANG', 'CRVAL1', 'CRVAL2']
    
    missing_keys = [key for key in required_keys if key not in header]
    if missing_keys:
        print(f"❌ Missing required keywords: {missing_keys}")
        print()
        return
    
    # Extract values
    stellalt = header['STELLALT']  # degrees
    stellaz = header['STELLAZ']    # degrees
    obslat = header['OBSLAT']      # degrees
    obslong = header['OBSLONG']    # degrees
    julian = header['JULIAN']      # Julian date
    lst = header['LST']            # Local sidereal time (hours)
    hourang = header['HOURANG']    # Hour angle (hours)
    crval1 = header['CRVAL1']      # RA (degrees)
    crval2 = header['CRVAL2']      # Dec (degrees)
    
    print("Mount and Observatory Data:")
    print(f"  Altitude: {stellalt:.6f}°")
    print(f"  Azimuth: {stellaz:.6f}°")
    print(f"  Observer Latitude: {obslat:.6f}°")
    print(f"  Observer Longitude: {obslong:.6f}°")
    print(f"  Julian Date: {julian:.8f}")
    print(f"  Local Sidereal Time: {lst:.8f} hours")
    print(f"  Hour Angle: {hourang:.8f} hours")
    print()
    
    print("WCS Solution:")
    print(f"  RA (CRVAL1): {crval1:.6f}°")
    print(f"  Dec (CRVAL2): {crval2:.6f}°")
    print()
    
    # Convert Alt/Az to RA/Dec using spherical astronomy
    try:
        calculated_ra, calculated_dec = altaz_to_radec(
            stellalt, stellaz, obslat, obslong, julian, lst
        )
        
        print("Calculated RA/Dec from Alt/Az:")
        print(f"  RA: {calculated_ra:.6f}°")
        print(f"  Dec: {calculated_dec:.6f}°")
        print()
        
        # Calculate differences
        ra_diff = abs(calculated_ra - crval1)
        dec_diff = abs(calculated_dec - crval2)
        
        # Handle RA wraparound (0°/360°)
        if ra_diff > 180:
            ra_diff = 360 - ra_diff
        
        print("Coordinate Differences:")
        print(f"  ΔRA: {ra_diff:.6f}° ({ra_diff * 3600:.1f} arcsec)")
        print(f"  ΔDec: {dec_diff:.6f}° ({dec_diff * 3600:.1f} arcsec)")
        
        # Convert to angular separation
        angular_sep = np.sqrt(
            (ra_diff * np.cos(np.radians(crval2)))**2 + dec_diff**2
        ) * 3600  # arcsec
        
        print(f"  Angular separation: {angular_sep:.1f} arcsec")
        print()
        
        # Evaluate consistency
        if angular_sep < 60:  # 1 arcminute
            print("✅ Excellent coordinate consistency (<1 arcmin)")
        elif angular_sep < 300:  # 5 arcminutes
            print("✅ Good coordinate consistency (<5 arcmin)")
        elif angular_sep < 1800:  # 30 arcminutes
            print("⚠️  Moderate coordinate consistency (<30 arcmin)")
        else:
            print("❌ Poor coordinate consistency (>30 arcmin)")
        
        # Check hour angle consistency
        ra_hours = crval1 / 15.0  # Convert RA to hours
        calculated_ha = lst - ra_hours
        
        # Normalize hour angle to [-12, 12] range
        while calculated_ha > 12:
            calculated_ha -= 24
        while calculated_ha < -12:
            calculated_ha += 24
        
        ha_diff = abs(calculated_ha - hourang)
        if ha_diff > 12:
            ha_diff = 24 - ha_diff
        
        print("Hour Angle Verification:")
        print(f"  Recorded hour angle: {hourang:.6f} hours")
        print(f"  Calculated hour angle: {calculated_ha:.6f} hours")
        print(f"  Difference: {ha_diff:.6f} hours ({ha_diff * 3600:.1f} arcsec)")
        
        if ha_diff < 0.01:  # ~36 arcsec
            print("✅ Hour angle is consistent")
        else:
            print("⚠️  Hour angle discrepancy detected")
        
    except Exception as e:
        print(f"❌ Error in coordinate conversion: {e}")
    
    print()

def altaz_to_radec(alt, az, lat, lon, julian, lst):
    """
    Convert Alt/Az coordinates to RA/Dec using spherical astronomy
    
    Parameters:
    -----------
    alt : float
        Altitude in degrees
    az : float
        Azimuth in degrees
    lat : float
        Observer latitude in degrees
    lon : float
        Observer longitude in degrees (not used in this calculation)
    julian : float
        Julian date (not used in this calculation)
    lst : float
        Local sidereal time in hours
    
    Returns:
    --------
    ra, dec : tuple
        Right ascension and declination in degrees
    """
    
    # Convert to radians
    alt_rad = np.radians(alt)
    az_rad = np.radians(az)
    lat_rad = np.radians(lat)
    
    # Calculate declination
    dec_rad = np.arcsin(
        np.sin(alt_rad) * np.sin(lat_rad) + 
        np.cos(alt_rad) * np.cos(lat_rad) * np.cos(az_rad)
    )
    
    # Calculate hour angle
    cos_ha = (np.sin(alt_rad) - np.sin(dec_rad) * np.sin(lat_rad)) / (
        np.cos(dec_rad) * np.cos(lat_rad)
    )
    
    # Clamp to valid range to avoid numerical errors
    cos_ha = np.clip(cos_ha, -1.0, 1.0)
    
    ha_rad = np.arccos(cos_ha)
    
    # Determine sign of hour angle based on azimuth
    if np.sin(az_rad) > 0:  # Azimuth > 180° (west)
        ha_rad = -ha_rad
    
    # Convert hour angle to hours
    ha_hours = ha_rad * 12.0 / np.pi
    
    # Calculate RA from LST and hour angle
    ra_hours = lst - ha_hours
    
    # Normalize RA to [0, 24) range
    while ra_hours < 0:
        ra_hours += 24
    while ra_hours >= 24:
        ra_hours -= 24
    
    # Convert to degrees
    ra_deg = ra_hours * 15.0
    dec_deg = np.degrees(dec_rad)
    
    return ra_deg, dec_deg

def check_astrometry_net_fields(header):
    """Check astrometry.net specific fields"""
    print("Astrometry.net Analysis:")
    print("-" * 20)
    
    # Check for astrometry.net signature
    anet_indicators = ['HISTORY Created by the Astrometry.net suite.',
                       'index id:', 'log odds:', 'nmatch:', 'scale:']
    
    is_anet = False
    for indicator in anet_indicators:
        for key in header:
            if isinstance(header[key], str) and indicator in header[key]:
                is_anet = True
                break
    
    if not is_anet:
        print("Not identified as astrometry.net solved")
        print()
        return
    
    print("✅ Identified as astrometry.net solved")
    
    # Extract solving statistics from comments
    for key in header:
        if key.startswith('COMMENT') and isinstance(header[key], str):
            comment = header[key]
            if 'log odds:' in comment:
                print(f"  {comment}")
            elif 'nmatch:' in comment:
                print(f"  {comment}")
            elif 'scale:' in comment:
                print(f"  {comment}")
            elif 'cpu time:' in comment:
                print(f"  {comment}")
    
    print()

def validate_coordinate_transform(wcs, image_shape):
    """Validate coordinate transformation at several points"""
    print("Coordinate Transform Validation:")
    print("-" * 35)
    
    # Test points: corners and center
    test_points = [
        (0, 0, "Bottom-left corner"),
        (image_shape[1]-1, 0, "Bottom-right corner"),
        (0, image_shape[0]-1, "Top-left corner"),
        (image_shape[1]-1, image_shape[0]-1, "Top-right corner"),
        (image_shape[1]/2, image_shape[0]/2, "Center")
    ]
    
    print("Testing coordinate transformations:")
    
    for x, y, description in test_points:
        try:
            # Forward transform: pixel to world
            ra, dec = wcs.pixel_to_world_values(x, y)
            
            # Reverse transform: world to pixel
            x_back, y_back = wcs.world_to_pixel_values(ra, dec)
            
            # Check round-trip accuracy
            dx = abs(x - x_back)
            dy = abs(y - y_back)
            
            print(f"  {description}:")
            print(f"    Pixel: ({x:.1f}, {y:.1f}) -> RA/Dec: ({ra:.6f}°, {dec:.6f}°)")
            print(f"    Round-trip error: Δx={dx:.3f}, Δy={dy:.3f} pixels")
            
            if max(dx, dy) > 0.1:
                print("    ⚠️  Large round-trip error")
            else:
                print("    ✅ Good round-trip accuracy")
            
        except Exception as e:
            print(f"    ❌ Transform failed: {e}")
    
    print()

def main():
    if len(sys.argv) != 2:
        print("Usage: python wcs_checker.py <fits_file_path>")
        sys.exit(1)
    
    fits_file = sys.argv[1]
    
    print("FITS WCS Consistency Checker")
    print("=" * 60)
    print()
    
    success = check_wcs_consistency(fits_file)
    
    if success:
        print("Analysis completed successfully!")
    else:
        print("Analysis failed - see errors above")
        sys.exit(1)

if __name__ == "__main__":
    main()
