#!/usr/bin/env python
##############################################################################
# MODULE:    i.hyper.specresamp
# AUTHOR(S): Created for hyperspectral spectral resampling
# PURPOSE:   Spectrally resample hyperspectral 3D raster to new wavelength sampling
# COPYRIGHT: (C) 2025 by the GRASS Development Team
# SPDX-License-Identifier: GPL-2.0-or-later
##############################################################################

# %module
# % description: Spectrally resample hyperspectral 3D raster to new wavelength sampling using Gaussian convolution
# % keyword: imagery
# % keyword: hyperspectral
# % keyword: resampling
# % keyword: spectral
# %end

# %option G_OPT_R3_INPUT
# % key: input
# % required: yes
# % description: Input hyperspectral 3D raster map (from i.hyper.import)
# % guisection: Input
# %end

# %option G_OPT_R3_OUTPUT
# % key: output
# % required: yes
# % description: Output resampled hyperspectral 3D raster map
# % guisection: Output
# %end

# %option
# % key: target_start
# % type: double
# % required: no
# % description: Target starting wavelength (nanometers)
# % guisection: Target Wavelengths
# %end

# %option
# % key: target_end
# % type: double
# % required: no
# % description: Target ending wavelength (nanometers)
# % guisection: Target Wavelengths
# %end

# %option
# % key: target_step
# % type: double
# % required: yes
# % description: Target wavelength step/sampling interval (nanometers)
# % guisection: Target Wavelengths
# %end

# %option
# % key: fwhm
# % type: double
# % required: no
# % description: Full Width at Half Maximum for Gaussian response (nm, default=target_step)
# % guisection: Processing
# %end

# %option
# % key: method
# % type: string
# % required: no
# % options: gaussian,nearest,linear
# % answer: gaussian
# % description: Resampling method
# % guisection: Processing
# %end

# %flag
# % key: n
# % description: Only use bands marked as valid (valid=1) in metadata
# % guisection: Processing
# %end

# %flag
# % key: i
# % description: Print information about input/output wavelengths without processing
# % guisection: Processing
# %end

import sys
import os
import numpy as np
import grass.script as gs
from grass.script import array as garray


def get_raster3d_info(raster3d):
    """Get information about 3D raster"""
    try:
        info = gs.raster3d_info(raster3d)
        return info
    except Exception as e:
        gs.fatal(f"Cannot get info for 3D raster {raster3d}: {e}")


def parse_wavelength_from_metadata(raster3d, band_num):
    """Parse wavelength and validity from band metadata"""
    band_name = f"{raster3d}#{band_num}"
    wavelength = None
    fwhm = None
    valid = True
    unit = "nm"
    
    try:
        result = gs.read_command('r.support', map=band_name, flags='n')
        
        for line in result.split('\n'):
            line = line.strip()
            if line.startswith('wavelength='):
                wavelength = float(line.split('=')[1])
            elif line.startswith('FWHM='):
                fwhm = float(line.split('=')[1])
            elif line.startswith('valid='):
                valid = int(line.split('=')[1]) == 1
            elif line.startswith('unit='):
                unit = line.split('=')[1].strip()
    except:
        pass
    
    return wavelength, fwhm, valid, unit


def convert_wavelength_to_nm(wavelength, unit):
    """Convert wavelength to nanometers"""
    unit = unit.lower().strip()
    
    if unit in ['nm', 'nanometer', 'nanometers']:
        return wavelength
    elif unit in ['um', 'µm', 'micrometer', 'micrometers', 'micron', 'microns']:
        return wavelength * 1000.0
    elif unit in ['m', 'meter', 'meters']:
        return wavelength * 1e9
    else:
        gs.warning(f"Unknown wavelength unit '{unit}', assuming nanometers")
        return wavelength


def get_all_band_wavelengths(raster3d, only_valid=False):
    """Extract all band wavelengths and metadata from 3D raster"""
    info = get_raster3d_info(raster3d)
    depths = int(info['depths'])
    
    bands = []
    
    gs.verbose(f"Scanning {depths} bands for wavelength metadata...")
    
    for i in range(1, depths + 1):
        wavelength, fwhm, valid, unit = parse_wavelength_from_metadata(raster3d, i)
        
        if wavelength is not None:
            wavelength_nm = convert_wavelength_to_nm(wavelength, unit)
            
            if only_valid and not valid:
                gs.verbose(f"Band {i}: {wavelength_nm} nm - SKIPPED (invalid)")
                continue
            
            bands.append({
                'band_num': i,
                'wavelength': wavelength_nm,
                'fwhm': fwhm if fwhm else 0,
                'valid': valid,
                'unit': unit
            })
            
            gs.verbose(f"Band {i}: {wavelength_nm} nm (FWHM: {fwhm}, valid: {valid})")
    
    if not bands:
        gs.fatal("No wavelength metadata found in 3D raster bands. "
                "Please use data imported with i.hyper.import or add wavelength metadata.")
    
    # Sort bands by wavelength
    bands.sort(key=lambda x: x['wavelength'])
    
    return bands


def gaussian_response(center_wl, wl_array, fwhm):
    """
    Generate Gaussian spectral response function
    
    Parameters:
    - center_wl: center wavelength of output band
    - wl_array: array of input wavelengths
    - fwhm: Full Width at Half Maximum
    """
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    response = np.exp(-0.5 * ((wl_array - center_wl) / sigma) ** 2)
    # Normalize so weights sum to 1
    response_sum = response.sum()
    if response_sum > 0:
        return response / response_sum
    else:
        return response


def nearest_neighbor_response(center_wl, wl_array):
    """Find nearest wavelength band"""
    idx = np.argmin(np.abs(wl_array - center_wl))
    response = np.zeros_like(wl_array)
    response[idx] = 1.0
    return response


def linear_interpolation_response(center_wl, wl_array):
    """Linear interpolation between two nearest bands"""
    if center_wl <= wl_array[0]:
        response = np.zeros_like(wl_array)
        response[0] = 1.0
        return response
    
    if center_wl >= wl_array[-1]:
        response = np.zeros_like(wl_array)
        response[-1] = 1.0
        return response
    
    # Find bracketing wavelengths
    idx_high = np.searchsorted(wl_array, center_wl)
    idx_low = idx_high - 1
    
    wl_low = wl_array[idx_low]
    wl_high = wl_array[idx_high]
    
    # Linear interpolation weights
    weight_high = (center_wl - wl_low) / (wl_high - wl_low)
    weight_low = 1.0 - weight_high
    
    response = np.zeros_like(wl_array)
    response[idx_low] = weight_low
    response[idx_high] = weight_high
    
    return response


def calculate_target_wavelengths(bands, target_start, target_end, target_step):
    """Calculate target output wavelengths"""
    input_wl = np.array([b['wavelength'] for b in bands])
    
    # If not specified, use input range
    if target_start is None:
        target_start = input_wl.min()
    if target_end is None:
        target_end = input_wl.max()
    
    # Generate target wavelengths
    target_wavelengths = np.arange(target_start, target_end + target_step, target_step)
    
    # Filter to only wavelengths within input range
    target_wavelengths = target_wavelengths[
        (target_wavelengths >= input_wl.min()) & 
        (target_wavelengths <= input_wl.max())
    ]
    
    if len(target_wavelengths) == 0:
        gs.fatal(f"No target wavelengths fall within input range "
                f"({input_wl.min():.1f} - {input_wl.max():.1f} nm)")
    
    return target_wavelengths


def print_resampling_info(input_bands, target_wavelengths, method, fwhm):
    """Print information about the resampling"""
    input_wl = np.array([b['wavelength'] for b in input_bands])
    
    gs.message("=" * 70)
    gs.message("Spectral Resampling Information:")
    gs.message("=" * 70)
    gs.message(f"Input bands: {len(input_bands)}")
    gs.message(f"Input wavelength range: {input_wl.min():.2f} - {input_wl.max():.2f} nm")
    gs.message(f"Input mean spacing: {np.mean(np.diff(input_wl)):.2f} nm")
    gs.message("")
    gs.message(f"Output bands: {len(target_wavelengths)}")
    gs.message(f"Output wavelength range: {target_wavelengths.min():.2f} - {target_wavelengths.max():.2f} nm")
    gs.message(f"Output spacing: {np.mean(np.diff(target_wavelengths)):.2f} nm")
    gs.message("")
    gs.message(f"Resampling method: {method}")
    if method == 'gaussian':
        gs.message(f"Gaussian FWHM: {fwhm:.2f} nm")
    gs.message("=" * 70)
    
    # Show first few and last few output wavelengths
    gs.message("Output wavelengths (nm):")
    if len(target_wavelengths) <= 10:
        gs.message("  " + ", ".join([f"{wl:.2f}" for wl in target_wavelengths]))
    else:
        gs.message("  " + ", ".join([f"{wl:.2f}" for wl in target_wavelengths[:5]]) + 
                  " ... " + 
                  ", ".join([f"{wl:.2f}" for wl in target_wavelengths[-5:]]))
    gs.message("=" * 70)


def resample_spectral(input_raster, input_bands, target_wavelengths, 
                     output_raster, method='gaussian', fwhm=None):
    """Perform spectral resampling"""
    
    # Get input wavelengths as array
    input_wl = np.array([b['wavelength'] for b in input_bands])
    n_input = len(input_bands)
    n_output = len(target_wavelengths)
    
    # Get raster info
    info = get_raster3d_info(input_raster)
    rows = int(info['rows'])
    cols = int(info['cols'])
    
    gs.message(f"Processing {rows} x {cols} pixels, {n_input} -> {n_output} bands...")
    gs.message(f"Total pixels to process: {rows * cols:,}")
    
    # Calculate response functions for each output band
    gs.verbose("Calculating spectral response functions...")
    responses = []
    for target_wl in target_wavelengths:
        if method == 'gaussian':
            response = gaussian_response(target_wl, input_wl, fwhm)
        elif method == 'nearest':
            response = nearest_neighbor_response(target_wl, input_wl)
        elif method == 'linear':
            response = linear_interpolation_response(target_wl, input_wl)
        responses.append(response)
    
    # Create temporary region for 3D raster operations
    gs.use_temp_region()
    gs.run_command('g.region', raster3d=input_raster)
    
    # Process each output band
    temp_maps = []
    
    for out_idx, (target_wl, response) in enumerate(zip(target_wavelengths, responses)):
        gs.percent(out_idx, n_output, 1)
        
        temp_map = f"tmp_specresamp_{os.getpid()}_{out_idx}"
        temp_maps.append(temp_map)
        
        # Build r.mapcalc expression for weighted sum
        terms = []
        for in_idx, weight in enumerate(response):
            if weight > 1e-6:  # Skip negligible weights
                band = input_bands[in_idx]
                band_name = f"{input_raster}#{band['band_num']}"
                terms.append(f"({weight:.10f} * {band_name})")
        
        if not terms:
            gs.fatal(f"No input bands contribute to output wavelength {target_wl:.2f} nm")
        
        expression = f"{temp_map} = " + " + ".join(terms)
        
        gs.run_command('r.mapcalc', expression=expression, quiet=True, overwrite=True)
        
        # Set metadata for this band
        gs.run_command('r.support', map=temp_map,
                      title=f"Band at {target_wl:.2f} nm",
                      units="reflectance", quiet=True)
        gs.write_command('r.support', map=temp_map, stdin=f"wavelength={target_wl}")
        gs.write_command('r.support', map=temp_map, stdin=f"unit=nm")
        gs.write_command('r.support', map=temp_map, stdin=f"valid=1")
        if method == 'gaussian':
            gs.write_command('r.support', map=temp_map, stdin=f"FWHM={fwhm}")
    
    gs.percent(1, 1, 1)
    
    # Stack temporary maps into 3D raster
    gs.message("Creating output 3D raster...")
    map_list = ",".join(temp_maps)
    
    gs.run_command('r.to.rast3', 
                  input=map_list,
                  output=output_raster,
                  overwrite=True)
    
    # Clean up temporary maps
    gs.message("Cleaning up temporary maps...")
    for temp_map in temp_maps:
        gs.run_command('g.remove', type='raster', name=temp_map, flags='f', quiet=True)
    
    # Set 3D raster metadata
    gs.run_command('r3.support', map=output_raster,
                  title="Spectrally resampled hyperspectral data",
                  history=f"Resampled from {input_raster} using {method} method")
    
    gs.message(f"Successfully created resampled 3D raster: {output_raster}")
    
    # Calculate and report statistics
    try:
        stats = gs.parse_command('r3.univar', map=output_raster, flags='g')
        gs.message("-" * 50)
        gs.message("Output Statistics:")
        gs.message(f"  Mean:   {float(stats['mean']):.4f}")
        gs.message(f"  Min:    {float(stats['min']):.4f}")
        gs.message(f"  Max:    {float(stats['max']):.4f}")
        gs.message(f"  StdDev: {float(stats['stddev']):.4f}")
        gs.message("-" * 50)
    except:
        pass


def main(options, flags):
    """Main function"""
    input_raster = options['input']
    output_raster = options['output']
    target_start = float(options['target_start']) if options['target_start'] else None
    target_end = float(options['target_end']) if options['target_end'] else None
    target_step = float(options['target_step'])
    fwhm = float(options['fwhm']) if options['fwhm'] else target_step
    method = options['method']
    only_valid = flags['n']
    info_only = flags['i']
    
    # Validate inputs
    if target_step <= 0:
        gs.fatal("Target wavelength step must be positive")
    
    if target_start and target_end and target_start >= target_end:
        gs.fatal("Target start wavelength must be less than target end wavelength")
    
    gs.message(f"Processing spectral resampling for: {input_raster}")
    
    # Get input band wavelengths
    input_bands = get_all_band_wavelengths(input_raster, only_valid=only_valid)
    
    # Calculate target wavelengths
    target_wavelengths = calculate_target_wavelengths(
        input_bands, target_start, target_end, target_step
    )
    
    # Print information
    print_resampling_info(input_bands, target_wavelengths, method, fwhm)
    
    if info_only:
        gs.message("Info mode: No output created.")
        return 0
    
    # Perform resampling
    resample_spectral(input_raster, input_bands, target_wavelengths,
                     output_raster, method=method, fwhm=fwhm)
    
    gs.message(f"Spectral resampling complete: {output_raster}")
    gs.message(f"Original bands: {len(input_bands)} -> Resampled bands: {len(target_wavelengths)}")
    
    return 0


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main(options, flags))
