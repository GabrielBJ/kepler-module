# Import the libraries needed to run the setup 

from setuptools import setup, find_packages

setup(name = "odekepler",
      description = "This produces a two-body simulation for the earth",
      author = "J. Gabriel Balarezo",
      author_email = "balarezog961@gmail.com",
      license = "BSD",
      version = "1.1",
      packages = find_packages(),
      install_requires = ['numpy', 'matplotlib', 'pillow', 'requests', 'scienceplots'])

