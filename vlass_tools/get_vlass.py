import os
import matplotlib.pyplot as plt
import requests
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import fits
import astropy.time
from astropy import table
from astropy.utils.data import download_file
import pandas as pd

archive_url = 'https://archive-new.nrao.edu/vlass'
tilelist_url = archive_url + '/VLASS_dyn_summary.php'

# TODO: use CADC cutout service
# http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync?ID=ad:VLASS/VLASS1.1.ql.T29t05.J112120%2B763000.10.2048.v1.I.iter1.image.pbcor.tt0.subim.fits&CIRCLE=168.3471998536797+76.18699791158396+0.01

def get_coverage(co, epoch=None, tilesfile='/Users/claw/data/vlass/VLASS_dyn_summary_21feb20.txt'):
    """ Given a ra, dec, return the name of a VLASS tile.
    Sexagesimal (hh:mm:ss, dd:mm:ss) coords expected.
    Can input coordinate as (RA, Dec) tuple or astropy SkyCoord.
    """

    if not isinstance(co, coord.SkyCoord):
        ra, dec = co
        if ra > 24:
            print('Warning: RA should be in hours')
        co = coord.SkyCoord(ra, dec, unit=(u.hour, u.deg))

    if os.path.exists(tilesfile):
        print(f'Using tiles file from disk {tilesfile}')
        with open(tilesfile) as fp:
            lines_all = fp.readlines()
        rows_imaged = filter(lambda x: 'imaged' in x, lines_all)
    else:
        print('Using requests to get tiles file')
        res = requests.get(tilelist_url)
        rows_imaged = filter(lambda x: 'imaged' in x,
                             res.content.decode().split('\n'))

    rows = []
    for row in rows_imaged:
        name, decmin, decmax, ramin, ramax, epoch0, date, *_ = row.split()
        if epoch is not None:
            if epoch != epoch0:
                continue
        if (co.ra.hour >= float(ramin)) and (co.ra.hour < float(ramax)) and \
           (co.dec.deg >= float(decmin)) and (co.dec.deg < float(decmax)):
            rows.append([name, decmin, decmax, ramin, ramax, epoch0, date])

    if len(rows):
        tab = table.Table(names=('name', 'decmin', 'decmax', 'ramin', 'ramax', 'epoch', 'date'), rows=rows)
        return tab
    else:
        return


def get_tilename(co):
    tab = get_coverage(co)
    assert len(tab) == 1
    return tab['name'].tolist()[0]


def get_filename(co, tilename=None):
    """ Given coord, find tile and fits directory for VLASS data.
    Calls get_tilename as intermediate (getting tile, e.g., T20t18).
    """

    if tilename is None:
        tilename = get_tilename(co)
    res = requests.get('{0}/quicklook/VLASS1.1/{1}/'
                       .format(archive_url, tilename))

    lines = filter(lambda x: 'href' in x and tilename in x,
                   res.content.decode().split('\n'))

    sep0 = 360*u.deg
    file0 = ''
    for line in lines:
        filename = line.split("\"")[1][:-1]
        center = filename.split(tilename)[-1].split('.')[1]
        co2 = coord.SkyCoord('{0}:{1}:{2}'
                             .format(center[1:3], center[3:5], center[5:7]),
                             '{0}:{1}:{2}'
                             .format(center[7:10], center[10:12], center[12:14]),
                             unit=(u.hour, u.deg))
        if co.separation(co2) < sep0:
            sep0 = co.separation(co2)
            file0 = filename

    assert (file0 != '') and (sep0 != 360*u.deg), 'No closest file found'

    print('Found a file {0} offset by {1}'.format(file0, sep0))

    return file0


def get_fitsname(co):
    """ Given coord, find fits file in VLASS
    """

    tilename = get_tilename(co)
    file0 = get_filename(co, tilename=tilename)
    base_url = '{0}/quicklook/VLASS1.1/{1}/{2}'.format(archive_url, tilename, file0)
    res = requests.get(base_url)

    fitsfile = list(filter(lambda x: file0 in x and 'tt0.subim.fits' in x,
                           res.content.decode().split('\n')))[0].split('\"')[1]

    return '{0}/{1}'.format(base_url, fitsfile)


def get_fits(co, outname=None):
    """ Given coord, download VLASS fits file.
    """

    fitsfile = get_fitsname(co)
    if outname is None:
        outname = fitsfile.split('/')[-1]

    if not os.path.exists(outname):
        image_file = download_file(fitsfile, cache=False)
        os.rename(image_file, outname)
    else:
        print('File {0} already downloaded'.format(outname))

    return outname


def parse_tab():
    ofektab = pd.read_table('ofek_table1.txt', delim_whitespace=True,
                            na_filter=False, comment='#', index_col=False,
                            names=['ra0', 'ra1', 'ra2', 'dec0', 'dec1', 'dec2',
                                   'Sf', 'eSf', 'L','z', 'offset'])
    for i in range(len(ofektab)):
        row = ofektab.iloc[i]
        co = coord.SkyCoord((row.ra0, row.ra1, row.ra2),
                            (row.dec0, row.dec1, row.dec2),
                            unit=(u.hour, u.deg))
        try:
            fitsname = get_fits(co)
            print(fitsname)
        except:
            pass
#        except FileNotFoundError:
#            print('Skipping coord at {0}'.format(co))


def make_reg():
    ofektab = pd.read_table('ofek_table1.txt', delim_whitespace=True,
                            na_filter=False, comment='#', index_col=False,
                            names=['ra0', 'ra1', 'ra2', 'dec0', 'dec1', 'dec2',
                                   'Sf', 'eSf', 'L', 'z', 'offset'])

    with open('ds9.reg', 'a') as fp:
        for i in range(len(ofektab)):
            row = ofektab.iloc[i]
            fp.write("circle({0}:{1}:{2},{3}:{4}:{5},30\")\n"
                     .format(int(row.ra0), int(row.ra1), row.ra2,
                             int(row.dec0), int(row.dec1), row.dec2))


def lc():
    t_first = astropy.time.Time(1994.7, format='decimalyear')
    t_nvss = astropy.time.Time(1995.4, format='decimalyear')
    t_vlass = astropy.time.Time(2018.0, format='decimalyear')
    s_first = 21.1*u.millijansky
    s_nvss = 18.5*u.millijansky
    s_vlass = 0.1*u.millijansky
    f_first = 1.4*u.gigahertz
    f_nvss = 1.4*u.gigahertz
    f_vlass = 3.0*u.gigahertz
    ss = lambda t, alpha: s_first*(t/t_first.value)**alpha

    plt.plot([t_first.value, t_nvss.value, t_vlass.value],
             [s_first.value, s_nvss.value, s_vlass.value], '.')
    plt.loglog()
    plt.show()
