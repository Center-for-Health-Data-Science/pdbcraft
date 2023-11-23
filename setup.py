#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    setup.py
#
#    pdbcraft setup.
#
#    Copyright (C) 2023 Valentina Sora 
#                       <sora.valentina1@gmail.com>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program. 
#    If not, see <http://www.gnu.org/licenses/>.


# Standard library
from setuptools import setup


# The name of the project
name = "pdbcraft"

# The URL where to find the project
url = \
    f"https://github.com/Center-for-Health-Data-Science/" \
    f"{name}"

# The project's author(s)
author = "Valentina Sora"

# The project's version
version = "0.0.1"

# A brief description of the project
description = \
    "A lightweight Python package for parsing and manipulating " \
    "PDB files."

# Which packages are included
packages = ["pdbcraft"]

# Run the setup
setup(name = name,
      url = url,
      author = author,
      version = version,
      description = description,
      packages = packages)