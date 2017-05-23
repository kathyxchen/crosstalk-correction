import os
from setuptools import setup


with open(os.path.join(os.path.dirname(__file__), "README.rst")) as readme:
    long_description = readme.read()

# Allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setup(
    name="crosstalk-correction",
    version="1.0.2",
    packages=["crosstalk_correction"],
    include_package_data=True,
    license="LICENSE",
    description="Python 3 implementation of maximum impact estimation "
                "(Donato et al., 2013)",
    long_description=long_description,
    url="https://github.com/kathyxchen/crosstalk-correction",
    author="Kathy Chen",
    author_email="chen.kathleenm@gmail.com",
    classifiers=[
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=["numpy"],
)
