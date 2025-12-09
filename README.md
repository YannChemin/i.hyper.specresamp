## DESCRIPTION

*i.hyper.specresamp* performs spectral resampling of hyperspectral 
imagery stored as 3D raster maps (`raster_3d`) from the 
*i.hyper* module family. The module resamples hyperspectral data to a 
new wavelength sampling interval, useful for harmonizing datasets from different 
sensors or reducing spectral dimensionality.

The module reads wavelength metadata from the input 3D raster bands and 
generates a new 3D raster with bands at the specified target wavelength 
intervals. Three resampling methods are available to accommodate different 
use cases and accuracy requirements.

Spectral resampling is commonly needed when:

- Combining data from sensors with different spectral samplings
- Reducing spectral dimensionality for computational efficiency
- Matching hyperspectral data to multispectral sensor specifications
- Preparing data for algorithms that require specific wavelength bands
- Simulating sensor characteristics for validation studies

## RESAMPLING METHODS

### Gaussian Convolution (default)

The **gaussian** method uses Gaussian-weighted spectral response functions 
centered at each target wavelength. This is the most physically realistic method 
as it mimics how real sensors integrate energy across their spectral bands. The 
Full Width at Half Maximum (FWHM) parameter controls the bandwidth of the 
Gaussian response. By default, FWHM equals the target wavelength step, but it 
can be adjusted to simulate specific sensor characteristics.

This method is recommended for most applications, especially when simulating 
multispectral sensors or when physical accuracy is important.

### Nearest Neighbor

The **nearest** method selects the input band closest to each target 
wavelength. This is the fastest method but produces discontinuous results. Use 
this for quick testing or when you simply want to extract specific wavelength 
bands without interpolation.

### Linear Interpolation

The **linear** method linearly interpolates between the two nearest input 
bands. This provides smooth results and is faster than Gaussian convolution, but 
doesn't account for realistic sensor response characteristics. Useful when 
computational speed is important and physical accuracy is less critical.

## NOTES

The module requires wavelength metadata in the input 3D raster. This metadata 
is automatically added by *i.hyper.import*. The required metadata fields 
are: **wavelength**, **unit**, and optionally **FWHM** and **valid**.

Output wavelengths are automatically constrained to the range of input 
wavelengths. If you specify target_start or target_end outside the input range, 
they will be adjusted to the input range boundaries.

The **-n** flag restricts processing to only bands marked as valid 
(`valid=1`) in the metadata. This is useful for excluding problematic 
bands affected by atmospheric absorption or sensor artifacts.

The **-i** flag provides information mode that displays the resampling 
plan without creating output. Use this to verify wavelength mappings and the 
number of output bands before processing.

### FWHM Parameter

The FWHM (Full Width at Half Maximum) parameter controls the bandwidth of the 
Gaussian response function. Smaller FWHM values produce narrower spectral 
responses (more selective), while larger values produce broader responses (more 
smoothing). The default FWHM equals the target wavelength step, which is 
appropriate for most applications.

For simulating specific sensors, set FWHM to match the sensor's spectral 
bandwidth. For example, Landsat 8 OLI bands have FWHM values ranging from 
60-180 nm.

### Memory Considerations

The module processes one output band at a time using *r.mapcalc*, which 
is memory-efficient for large datasets. Temporary 2D rasters are created for each 
output band and then stacked into the final 3D raster using *r.to.rast3*.

## EXAMPLES

```bash
# Resample PRISMA data (400-2400nm @ ~3nm) to 10nm spacing
i.hyper.specresamp input=prisma \
                   output=prisma_10nm \
                   target_step=10

# This will create output from 400-2400nm at 10nm intervals
# Using Gaussian convolution with FWHM=10nm
```

```bash
# Resample to coarser 30nm spacing for faster processing
i.hyper.specresamp input=enmap \
                   output=enmap_30nm \
                   target_step=30 \
                   method=gaussian \
                   fwhm=30
```

```bash
# Extract only VNIR range (400-1000nm) at 5nm spacing
i.hyper.specresamp input=prisma \
                   output=prisma_vnir_5nm \
                   target_start=400 \
                   target_end=1000 \
                   target_step=5 \
                   method=gaussian
```

```bash
# Quick band selection using nearest neighbor (fast)
i.hyper.specresamp input=tanager \
                   output=tanager_20nm \
                   target_step=20 \
                   method=nearest
```

```bash
# Simulate Landsat 8 OLI band 5 (NIR: 851-879nm, center=865nm)
# Using narrow range and appropriate FWHM
i.hyper.specresamp input=prisma \
                   output=prisma_oli_b5 \
                   target_start=865 \
                   target_end=865 \
                   target_step=1 \
                   method=gaussian \
                   fwhm=28
```

```bash
# Preview resampling plan before processing
i.hyper.specresamp input=enmap \
                   output=test \
                   target_step=15 \
                   -i

# Output shows:
# - Input wavelength range and number of bands
# - Output wavelength range and number of bands
# - Resampling method parameters
# - List of output wavelengths
```

```bash
# Use only valid bands and linear interpolation
i.hyper.specresamp input=prisma \
                   output=prisma_12nm_valid \
                   target_step=12 \
                   method=linear \
                   -n
```

```bash
# Complete workflow: import, resample, calculate albedo
# Step 1: Import hyperspectral data
i.hyper.import input=/data/PRISMA.he5 \
               product=prisma \
               output=prisma_full

# Step 2: Resample to 10nm for faster processing
i.hyper.specresamp input=prisma_full \
                   output=prisma_10nm \
                   target_step=10

# Step 3: Calculate albedo from resampled data
i.hyper.albedo input=prisma_10nm \
               output=prisma_albedo \
               weighting=solar

# Step 4: Compare processing time vs original
# (resampled data processes much faster)
```

```bash
# Resample SWIR region to match VNIR sampling
# If VNIR is at 5nm and SWIR at 10nm, standardize to 10nm
i.hyper.specresamp input=enmap_swir \
                   output=enmap_swir_10nm \
                   target_step=10 \
                   method=gaussian \
                   fwhm=10
```

## SEE ALSO

[i.hyper.import](i.hyper.import.html),
[i.hyper.albedo](i.hyper.albedo.html),
[i.hyper.rgb](i.hyper.rgb.html),
[i.hyper.composite](i.hyper.composite.html),
[i.hyper.preproc](i.hyper.preproc.html),
[r.mapcalc](https://grass.osgeo.org/grass-stable/manuals/r.mapcalc.html),
[r.to.rast3](https://grass.osgeo.org/grass-stable/manuals/r.to.rast3.html),
[r3.univar](https://grass.osgeo.org/grass-stable/manuals/r3.univar.html),
[r3.support](https://grass.osgeo.org/grass-stable/manuals/r3.support.html)

## REFERENCES

- Guanter, L., et al. (2015). The EnMAP spaceborne imaging spectroscopy 
  mission for earth observation. *Remote Sensing*, 7(7), 8830-8857.
- Loizzo, R., et al. (2018). PRISMA: The Italian hyperspectral mission. 
  *IGARSS 2018-2018 IEEE International Geoscience and Remote Sensing 
  Symposium*, 175-178.
- Chander, G., Markham, B. L., & Helder, D. L. (2009). Summary of current 
  radiometric calibration coefficients for Landsat MSS, TM, ETM+, and EO-1 
  ALI sensors. *Remote Sensing of Environment*, 113(5), 893-903.
- Huete, A., et al. (2002). Overview of the radiometric and biophysical 
  performance of the MODIS vegetation indices. *Remote Sensing of 
  Environment*, 83(1-2), 195-213.

## AUTHORS

Created for the i.hyper module family

Based on spectral resampling concepts from hyperspectral remote sensing community
