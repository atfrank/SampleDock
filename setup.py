import os
from setuptools import setup, find_packages

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

with open('LICENSE') as f:
    license = f.read()

setup(
    name='SampleDock',
    version='0.1.0',
    description='Molecular design framework the merges generative AI and molecular docking',
    author='Aaron Frank and Ziqiao Xu',
    author_email='afrankz@umich.edu',
    url='https://github.com/atfrank/SampleDock',
    license=license,
    packages=find_packages()
)