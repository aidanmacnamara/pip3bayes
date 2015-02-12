from setuptools import setup, find_packages
import os
import glob

datadir = os.path.join('share','data')
datafiles = [(datadir, [f for f in glob.glob(os.path.join(datadir, '*'))])]

setup(
    name='pip3bayes',
    version='0.1.3',
    description='Using BayesSB with chemical dimerization',
    url='',
    author='Aidan MacNamara',
    author_email='aidan.macnamara@ebi.ac.uk',
    license='',
    zip_safe=False,
    package_dir = {'':'src'},
    packages = ["pip3bayes"],
    data_files = datafiles,
)