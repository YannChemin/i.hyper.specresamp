"""Tests that libras3d is functional for i.hyper.specresamp standalone mode."""
import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_ras3d_common import (
    WYVERN_PATH, TANAGER_PATH,
    skip_without_ras3d, skip_without_wyvern, skip_without_tanager,
    open_cube_checked, assert_band_valid, install_ras3d_shim, make_wl_sidecar,
)

@skip_without_ras3d
def test_shim_installs():
    install_ras3d_shim()
    import grass.script as gs
    assert hasattr(gs, 'raster3d_info')

@skip_without_ras3d
@skip_without_wyvern
def test_extract_z_slice_ras3d(tmp_path):
    install_ras3d_shim()
    os.environ['RAS3D_OUTDIR'] = str(tmp_path)
    sys.path.insert(0, '/home/yann/dev/i.hyper.specresamp')
    from i_hyper_specresamp import extract_z_slice
    extract_z_slice(WYVERN_PATH, 5, 'specr_slice_5')
    from ras3d_grass_shim import get_band_cache
    assert 'specr_slice_5' in get_band_cache()
    assert_band_valid(get_band_cache()['specr_slice_5'], 'specresamp slice 5')

@skip_without_ras3d
@skip_without_wyvern
def test_wavelength_sidecar():
    install_ras3d_shim()
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    sidecar, _ = make_wl_sidecar(WYVERN_PATH, r['depths'])
    ras3d.close_cube(h)
    sys.path.insert(0, '/home/yann/dev/i.hyper.specresamp')
    from i_hyper_specresamp import get_all_band_wavelengths
    bands = get_all_band_wavelengths(WYVERN_PATH)
    assert len(bands) == r['depths']
    os.unlink(sidecar)

@skip_without_ras3d
@skip_without_tanager
def test_open_hdf5_and_bands():
    import ras3d
    h, r = open_cube_checked(TANAGER_PATH)
    arr = ras3d.get_band(h, 0)
    assert_band_valid(arr, 'Tanager band 0')
    arr = ras3d.get_band(h, 212)
    assert_band_valid(arr, 'Tanager band 212')
    ras3d.close_cube(h)
