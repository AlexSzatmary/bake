#!/usr/bin/python

from bake import __version__

from setuptools import setup, find_packages
setup(
    name = "bake",
    version = __version__,
    packages = find_packages(),
    entry_points="""
    # -*- Entry points: -*-
    [console_scripts]
    bake=bake.cmdline:main
    """,
    author = "Alex C. Szatmary",
    author_email = "alex.szatmary@gmail.com",
    description = "bake runs the same code a number of times, varying "\
	"parameters.",
    license = "MIT",
)
