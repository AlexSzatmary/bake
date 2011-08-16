#!/usr/bin/python

from setuptools import setup, find_packages
setup(
    name = "bake",
    version = "0.0.1",
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
