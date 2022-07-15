import os
from setuptools import setup, find_packages


PRODUCT_NAME = 'mindyourps'
with open('requirements.txt') as reqs:
    DEPENDENCIES = reqs.read().strip().split('\n')

VERSION = os.environ.get('GITVERSION_MAJORMINORPATCH', '1.0.0')

setup(
    name=PRODUCT_NAME,
    version=VERSION,
    author='Penedos et al 2022',
    description=(
        'Functions for excluding relatedness between samples based on sequence'
        ' and date. Manuscript at https://www.thelancet.com/journals/ebiom/'
        'article/PIIS2352-3964(22)00173-6/fulltext?rss=yes'),
    install_requires=DEPENDENCIES,
    packages=find_packages()
)
