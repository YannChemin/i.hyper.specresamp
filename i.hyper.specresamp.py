#!/usr/bin/env python
##############################################################################
# MODULE:    i.hyper.specresamp
# AUTHOR(S): Created for hyperspectral spectral resampling
# PURPOSE:   Spectrally resample hyperspectral 3D raster to new wavelength sampling
# COPYRIGHT: (C) 2025 by the GRASS Development Team
# SPDX-License-Identifier: GPL-2.0-or-later
##############################################################################

# ── ras3d standalone detection ────────────────────────────────────────────────
import os as _os
_RAS3D = False
if not _os.environ.get('GISBASE'):
    try:
        import importlib.util as _ilu
        if _ilu.find_spec('ras3d') and _ilu.find_spec('ras3d_grass_shim'):
            from ras3d_grass_shim import install as _r3_install
            _r3_install()
            _RAS3D = True
    except Exception:
        pass
# ─────────────────────────────────────────────────────────────────────────────

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
import ctypes
import numpy as np
import grass.script as gs
from grass.script import array as garray

# ---------------------------------------------------------------------------
# Fast Z-slice extraction via Rast3d_extract_z_slice() (ctypes)
# ---------------------------------------------------------------------------

_raster3d_lib = None


def _load_raster3d_lib():
    """Load libgrass_raster3d and its deps via ctypes (once per process)."""
    global _raster3d_lib
    if _raster3d_lib is not None:
        return _raster3d_lib

    gisbase = os.environ["GISBASE"]
    libdir = os.path.join(gisbase, "lib")

    # Load dependency chain with RTLD_GLOBAL so raster3d resolves their symbols.
    for name in ("libgrass_gis.so", "libgrass_raster.so"):
        ctypes.CDLL(os.path.join(libdir, name), ctypes.RTLD_GLOBAL)

    lib = ctypes.CDLL(
        os.path.join(libdir, "libgrass_raster3d.so"), ctypes.RTLD_GLOBAL
    )
    lib.Rast3d_extract_z_slice.restype = ctypes.c_int
    lib.Rast3d_extract_z_slice.argtypes = [
        ctypes.c_char_p,  # name3d
        ctypes.c_char_p,  # mapset3d ("" = search)
        ctypes.c_int,     # z  (0-based)
        ctypes.c_char_p,  # name2d
    ]

    # Initialise GRASS C environment for this process.
    libgis = ctypes.CDLL(os.path.join(libdir, "libgrass_gis.so"))
    libgis.G_gisinit(b"i.hyper.specresamp")

    _raster3d_lib = lib
    return lib


def extract_z_slice(name3d, band_num_1based, name2d):
    """Extract one band (1-based) from a 3D raster to a 2D raster.

    Uses Rast3d_extract_z_slice() which opens the map with RASTER3D_NO_CACHE
    and calls Rast3d_get_block() for tile-bulk reads — each tile is loaded
    exactly once instead of one function call per voxel.
    """
    if _RAS3D:
        import ras3d as _r3, ras3d_write as _r3w
        _h = _r3.open_cube(name3d)
        _arr = _r3.get_band(_h, band_num_1based - 1)
        from ras3d_grass_shim import get_band_cache
        get_band_cache()[name2d] = _arr
        _r3w.write_raster2d(_r3w.outpath(name2d), _arr, _h)
        _r3.close_cube(_h)
        return
    lib = _load_raster3d_lib()
    z = band_num_1based - 1  # convert 1-based band to 0-based z index
    ret = lib.Rast3d_extract_z_slice(
        name3d.encode(), b"", ctypes.c_int(z), name2d.encode()
    )
    if ret != 0:
        gs.fatal(
            f"Rast3d_extract_z_slice failed for band {band_num_1based} of {name3d}"
        )


def get_raster3d_info(raster3d):
    """Get information about 3D raster"""
    try:
        info = gs.raster3d_info(raster3d)
        return info
    except Exception as e:
        gs.fatal(f"Cannot get info for 3D raster {raster3d}: {e}")


_HYPER_JSON_BAND_CACHE = {}


def _load_hyper_json_bands(raster3d):
    """Read wavelength/fwhm/validity from i.hyper.import's JSON sidecar.

    i.hyper.import (HyperMetadata) stores band metadata at
    $MAPSET/grid3/<mapname>/hyper.json rather than in r3.support history,
    so that must be checked before falling back to per-band r.support probing.
    """
    if raster3d in _HYPER_JSON_BAND_CACHE:
        return _HYPER_JSON_BAND_CACHE[raster3d]

    import json as _json

    name, mapset = (raster3d.split('@', 1) if '@' in raster3d
                     else (raster3d, None))
    result = {}
    try:
        env = gs.gisenv()
        mapset = mapset or env['MAPSET']
        path = os.path.join(env['GISDBASE'], env['LOCATION_NAME'], mapset,
                            'grid3', name, 'hyper.json')
        if os.path.isfile(path):
            with open(path) as _fj:
                data = _json.load(_fj)
            b = data.get('bands') or {}
            wavelengths = b.get('wavelength')
            if wavelengths:
                fwhms = b.get('fwhm') or []
                valids = b.get('validity') or []
                for i, wl in enumerate(wavelengths):
                    result[i + 1] = (
                        float(wl),
                        float(fwhms[i]) if i < len(fwhms) else None,
                        bool(valids[i]) if i < len(valids) else True,
                        'nm',
                    )
    except Exception:
        result = {}

    _HYPER_JSON_BAND_CACHE[raster3d] = result
    return result


def parse_wavelength_from_metadata(raster3d, band_num):
    """Parse wavelength and validity from band metadata"""
    json_bands = _load_hyper_json_bands(raster3d)
    if band_num in json_bands:
        return json_bands[band_num]

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
    if _RAS3D:
        import json as _json
        for _sfx in ('', '.tif', '.tiff', '.h5', '.hdf5'):
            _base = raster3d.removesuffix(_sfx) if raster3d.endswith(_sfx) else raster3d
            _wlp = _base + '.wl.json'
            if _os.path.exists(_wlp):
                with open(_wlp) as _f:
                    _wl = _json.load(_f)
                _wl_nm = [w * 1000 if w < 10 else w for w in _wl]
                return [{'band_num': i+1, 'wavelength': wl, 'fwhm': None, 'valid': 1}
                        for i, wl in enumerate(_wl_nm)]
        import ras3d as _r3
        _h = _r3.open_cube(raster3d); _r = _r3.get_region(_h); _r3.close_cube(_h)
        return [{'band_num': i+1, 'wavelength': float(i+1), 'fwhm': None, 'valid': 1}
                for i in range(_r['depths'])]
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
    
    if _RAS3D:
        import ras3d as _r3
        import ras3d_write as _r3w
        from ras3d_grass_shim import get_band_cache
        # Pre-load all needed input slices into numpy arrays
        needed_in_indices = {
            in_idx
            for response in responses
            for in_idx, w in enumerate(response)
            if w > 1e-6
        }
        gs.message(f"[ras3d] Pre-loading {len(needed_in_indices)} input slices...")
        _h = _r3.open_cube(input_raster)
        slice_arrays = {}
        for in_idx in sorted(needed_in_indices):
            slice_arrays[in_idx] = _r3.get_band(_h, input_bands[in_idx]['band_num'] - 1)
        _r3.close_cube(_h)
        # Compute weighted sums in numpy and stack into cube
        gs.message(f"[ras3d] Computing {n_output} output bands...")
        out_bands = []
        _h2 = _r3.open_cube(input_raster)
        for out_idx, (target_wl, response) in enumerate(zip(target_wavelengths, responses)):
            gs.percent(out_idx, n_output, 1)
            accum = None
            for in_idx, weight in enumerate(response):
                if weight > 1e-6:
                    contrib = weight * slice_arrays[in_idx]
                    accum = contrib if accum is None else accum + contrib
            if accum is None:
                gs.fatal(f"No input bands contribute to output wavelength {target_wl:.2f} nm")
            out_bands.append(accum)
        gs.percent(1, 1, 1)
        gs.message("[ras3d] Writing output 3D cube...")
        _r3w.write_raster3d(
            _r3w.outpath(output_raster),
            out_bands,
            _h2,
            [float(wl) for wl in target_wavelengths],
        )
        _r3.close_cube(_h2)
        gs.message(f"Successfully created resampled 3D raster: {output_raster}")
    else:
        # Create temporary region for 3D raster operations
        gs.use_temp_region()
        gs.run_command('g.region', raster3d=input_raster)

        # ------------------------------------------------------------------
        # Pre-extract all needed input slices once via Rast3d_extract_z_slice.
        # This uses RASTER3D_NO_CACHE + Rast3d_get_block() so each tile at a
        # given Z level is read exactly once (tile-bulk path), instead of one
        # function call per voxel through r.mapcalc's #{band_num} accessor.
        # ------------------------------------------------------------------
        needed_in_indices = {
            in_idx
            for response in responses
            for in_idx, w in enumerate(response)
            if w > 1e-6
        }
        gs.message(f"Pre-extracting {len(needed_in_indices)} input slices...")
        slice_maps = {}  # in_idx -> temp 2D raster name
        for in_idx in sorted(needed_in_indices):
            band = input_bands[in_idx]
            slice_name = f"tmp_specresamp_in_{os.getpid()}_{in_idx}"
            extract_z_slice(input_raster, band['band_num'], slice_name)
            slice_maps[in_idx] = slice_name

        # Process each output band
        temp_maps = []

        for out_idx, (target_wl, response) in enumerate(zip(target_wavelengths, responses)):
            gs.percent(out_idx, n_output, 1)

            temp_map = f"tmp_specresamp_{os.getpid()}_{out_idx}"
            temp_maps.append(temp_map)

            # Build r.mapcalc expression for weighted sum using pre-extracted 2D slices
            terms = []
            for in_idx, weight in enumerate(response):
                if weight > 1e-6:  # Skip negligible weights
                    terms.append(f"({weight:.10f} * {slice_maps[in_idx]})")

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

        # Clean up input slice temp maps
        gs.message("Cleaning up temporary maps...")
        for name in slice_maps.values():
            gs.run_command('g.remove', type='raster', name=name, flags='f', quiet=True)
        # Clean up output band temp maps
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
