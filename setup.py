from setuptools import setup, find_packages

from pyLip import __version__

setup(name='pyLip',
      version=__version__,
      description='Python library for processing individual mass spectra',
      url='https://github.com/alexandrovteam/pyLip',
      author='Andrew Palmer, EMBL',
      packages=find_packages())
