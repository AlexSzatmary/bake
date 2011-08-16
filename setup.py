#!/usr/bin/python

from setuptools import setup, find_packages
setup(
    name = "bake",
    version = "0.0.0",
    packages = find_packages(),
    entry_points="""
    # -*- Entry points: -*-
    [console_scripts]
    bake=bake.cmdline:main
    """
)
