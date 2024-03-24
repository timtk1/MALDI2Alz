from setuptools import setup, find_packages

from pyImagingMSpec import __version__

setup(name='pyImagingMSpec',
      version=__version__,
      description='Python library for processing imaging mass spectrometry data',
      url='https://github.com/alexandrovteam/pyImagingMSpec',
      author='Alexandrov Team, EMBL',
      packages=find_packages())
