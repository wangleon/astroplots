#!/usr/bin/env python
from distutils.core import setup

setup(
    name         = 'keplermap',
    version      = '0.1',
    description  = 'Plot the skymap of Kepler fields',
    author       = 'Liang Wang',
    author_email = 'lwang@mpe.mpg.de',
    license      = 'BSD',
    packages     = [
                    'keplermap',
                   ],
    package_data = {'keplermap':['data/*']},
    )
