#! /usr/bin/env python
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name="neighbors",
    version='0.1',
    author="Semyeong Oh",
    author_email="semyeong.oh@gmail.com",
    packages=["neighbors"],
    url="",
    description="",
    # long_description=rd("README.md"),
    # install_requires=[],
    include_package_data=True
)
