#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    _util.py
#
#    Private miscellanea utilities.
#
#    Copyright (C) 2025 Valentina Sora 
#                       <sora.valentina1@gmail.com>
#
#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later
#    version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#    See the GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General
#    Public License along with this program; if not, write to the
#    Free Software Foundation, Inc., 51 Franklin Street, Fifth
#    Floor, Boston, MA  02110-1301, USA.


#######################################################################


# Set the module's description.
__doc__ = "Private miscellanea utilities."


#######################################################################


def sort_keys_at_depths(d,
                        depths):
    """Recursively sort the keys of a dictionary at specific depths.
    
    Parameters
    ----------
    d : :class:`dict`
        The dictionary to be processed.

    depths : :class:`list`
        A list of depths at which to sort the keys.

    Returns
    -------
    d : :class:`dict`
        The dictionary with sorted keys at the specified depths.
    """

    # Define the recursive step.
    def recurse(d,
                depths,
                current_depth):

        # If the current object is not a dictionary
        if not isinstance(d, dict):

            # Simply return the object.
            return d

        # If the current depth is one of the target depths
        if current_depth in depths:

            # Sort the dictionary keys.
            return dict(\
                        sorted(\
                            (k, \
                             recurse(\
                                d = v,
                                depths = depths,
                                current_depth = current_depth + 1)) \
                            for k, v in d.items()))
        # Otherwise       
        else:

            # Keep the current keys' order and recurse into the values.
            return {k: recurse(d = v, 
                               depths = depths,
                               current_depth = current_depth + 1) \
                    for k, v in d.items()}

    #-----------------------------------------------------------------#

    # Return the result of the recursive step.
    return recurse(d = d,
                   depths = depths,
                   current_depth = 1)
