import io

from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from gwcs import wcs
from gwcs import coordinate_frames as cf

import numpy as np
from astropy.io import fits
import asdf
from asdf import fits_embed, AsdfFile
from asdf.asdf import open_asdf


hdul = fits.open("wepkrgS20170927S0015.fits", mode='update')

shift = models.Shift(-1.*u.pix) & models.Shift(-2.*u.pix)

#shift = models.Chebyshev2D(2, 2)

#shift = models.Polynomial2D(2, 2)

detector_frame = cf.Frame2D(name="detector", axes_names=("x", "y"),
                            unit=(u.pix, u.pix))
sky_frame = cf.CelestialFrame(reference_frame=coord.ICRS(), name='icrs',
                              unit=(u.deg, u.deg))

pipeline = [(detector_frame, shift),
            (sky_frame, None)]
wcsobj = wcs.WCS(pipeline)

tree = {"wcs": wcsobj}

#fa = fits_embed.AsdfInFits(hdul, tree)
#fa = fits_embed.AsdfInFits(fits.HDUList(), tree)
#fa.write_to('test.fits')
af = AsdfFile(tree)

fd = io.BytesIO()

#af._write_tree(tree, fd, pad_blocks=False)
af.write_to(fd)
fd.seek(0)
wcsbytes = fd.read()
fd.close()
# af.write_to('testfull.asdf')
af.close()

fmt = 'A{0}'.format(len(wcsbytes))

col = fits.Column(name='gWCS', format=fmt, array=np.array([wcsbytes]),
                  ascii=True)
#print(col.array)
hdu = fits.TableHDU.from_columns([col], name='WCS')

hdul.append(hdu)

hdul.flush()

# Read back in:

hdul = fits.open("wepkrgS20170927S0015.fits")

hdu = hdul['WCS']

col = hdu.data['gWCS']

wcsbytes = col.tobytes()

wcstext = wcsbytes.decode('utf-32')
wcstext = wcstext[:-1] + '\n'
#print('[{0}]'.format(wcstext[-1]))

#fd.write('\n'.encode('utf-8'))

fd = io.BytesIO()
# #fd.write('#ASDF 1.0.0\n'.encode('utf-8'))
# #print(wcsbytes.decode('utf-8'))
fd.write(wcstext.encode('utf-8'))
fd.seek(0)
#print(fd.read().decode())

#fd2 = open('testenc.txt', 'wb')
#fd2.write(fd.read())
#fd2.close()

# # Show that ASDF file is actually valid:
# # - This works properly when the ASDF version is correct.
# # af = open_asdf('testfull.asdf')
# #fd = asdf.generic_io.RealFile('testfull.asdf', 'r')
# #fd = asdf.generic_io.get_file('testfull.asdf', 'r')
# fd = asdf.generic_io.get_file(fd, 'r')
# print(fd)
# print(fd.read())
# #af = open_asdf(fd)
# #print(af.tree)

# # Open valid ASDF file and convert it to a file-like Bytes object, the same as
# # from the problematic FITS extension. This works.
# fd = io.BytesIO()
# ogg = open('testfull.asdf', 'rb')
# # fd.write(ogg.read().encode('utf-8'))
# fd.write(ogg.read())
# fd.seek(0)
# ogg.close()
# #print(fd.read().decode('utf-8'))
# #print(fd.read())


# # Try parsing the file-like object as ASDF:
# # print(type(fd))
af = open_asdf(fd)
# # #af = open_asdf('testfull.asdf')
print(af.tree)


#print(hdu.name)

# t = fa.tree.to_tree(fa.tree['wcs'], )
# print(type(t))
# print(dir(t))
# print(t)
# print('FAT', fa.tree)


#/home/jturner/anaconda3/envs/wavecal/lib/python3.6/site-packages/asdf/fits_embed.py
#fa.write_to("testit.fits", overwrite=True)

#fits.info("imaging_with_wcs_in_asdf.fits")


# %YAML 1.1
# %TAG ! tag:stsci.edu:asdf/
# --- !core/asdf-1.1.0
# asdf_library: !core/software-1.0.0 {author: Space Telescope Science Institute, homepage: 'http://github.com/spacetelescope/asdf',
#   name: asdf, version: 2.4.3a1.dev10+g1b41c6d.d20191029}
# history:
#   extensions:
#   - !core/extension_metadata-1.0.0
#     extension_class: asdf.extension.BuiltinExtension
#     software: {name: asdf, version: 2.4.3a1.dev10+g1b41c6d.d20191029}
#   - !core/extension_metadata-1.0.0
#     extension_class: astropy.io.misc.asdf.extension.AstropyAsdfExtension
#     software: {name: astropy, version: 4.0.dev25566}
#   - !core/extension_metadata-1.0.0
#     extension_class: gwcs.extension.GWCSExtension
#     software: {name: gwcs, version: 0.11.1a1.dev50+g2fc6944}
#   - !core/extension_metadata-1.0.0
#     extension_class: astropy.io.misc.asdf.extension.AstropyExtension
#     software: {name: astropy, version: 4.0.dev25566}
# wcs: !<tag:stsci.edu:gwcs/wcs-1.0.0>
#   name: ''
#   steps:
#   - !<tag:stsci.edu:gwcs/step-1.0.0>
#     frame: !<tag:stsci.edu:gwcs/frame2d-1.0.0>
#       axes_names: [x, y]
#       axes_order: [0, 1]
#       axis_physical_types: ['custom:x', 'custom:y']
#       name: detector
#       unit: [!unit/unit-1.0.0 pixel, !unit/unit-1.0.0 pixel]
#     transform: !transform/concatenate-1.1.0
#       forward:
#       - !transform/shift-1.2.0
#         offset: !unit/quantity-1.1.0 {unit: !unit/unit-1.0.0 pixel, value: -1.0}
#       - !transform/shift-1.2.0
#         offset: !unit/quantity-1.1.0 {unit: !unit/unit-1.0.0 pixel, value: -2.0}
#   - !<tag:stsci.edu:gwcs/step-1.0.0>
#     frame: !<tag:stsci.edu:gwcs/celestial_frame-1.0.0>
#       axes_names: [lon, lat]
#       axes_order: [0, 1]
#       axis_physical_types: [pos.eq.ra, pos.eq.dec]
#       name: icrs
#       reference_frame: !<tag:astropy.org:astropy/coordinates/frames/icrs-1.1.0>
#         frame_attributes: {}
#       unit: [!unit/unit-1.0.0 deg, !unit/unit-1.0.0 deg]

# #ASDF 1.0.0
# #ASDF_STANDARD 1.2.0
# %YAML 1.1
# ---
# wcs: !<tag:stsci.edu:gwcs/wcs-1.0.0>
#   name: ''
#   steps:
#   - !<tag:stsci.edu:gwcs/step-1.0.0>
#     frame: !<tag:stsci.edu:gwcs/frame2d-1.0.0>
#       axes_names: [x, y]
#       name: detector
#       unit: [!<tag:stsci.edu:asdf/unit/unit-1.0.0> pixel, !<tag:stsci.edu:asdf/unit/unit-1.0.0> pixel]
#     transform: !<tag:stsci.edu:asdf/transform/concatenate-1.1.0>
#       forward:
#       - !<tag:stsci.edu:asdf/transform/shift-1.2.0>
#         offset: !<tag:stsci.edu:asdf/unit/quantity-1.1.0> {unit: !<tag:stsci.edu:asdf/unit/unit-1.0.0> pixel,
#           value: -1.0}
#       - !<tag:stsci.edu:asdf/transform/shift-1.2.0>
#         offset: !<tag:stsci.edu:asdf/unit/quantity-1.1.0> {unit: !<tag:stsci.edu:asdf/unit/unit-1.0.0> pixel,
#           value: -2.0}
#   - !<tag:stsci.edu:gwcs/step-1.0.0>
#     frame: !<tag:stsci.edu:gwcs/celestial_frame-1.0.0>
#       axes_names: [lon, lat]
#       name: icrs
#       reference_frame: {type: ICRS}
#       unit: [!<tag:stsci.edu:asdf/unit/unit-1.0.0> deg, !<tag:stsci.edu:asdf/unit/unit-1.0.0> deg]
