#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    mmcifparsers.py
#
#    Parsers for mmCIF files.
#
#    Copyright (C) 2024 Valentina Sora 
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
__doc__ = "Parsers for mmCIF files."


#######################################################################


# Import from the standard library.
from collections import defaultdict
import copy
import logging as log
import os
import re
import shlex
# Import from 'pdbcraft'.
from .structure import Structure
from . import _defaults


#######################################################################


# Get the module's logger.
logger = log.getLogger(__name__)


#######################################################################


class MMCIFParser:

    """
    A parser for mmCIF files.
    """

    # Initialize the number to be used for unnamed models in the
    # structure built from the mmCIF file.
    UNKNOWN_MODEL = 1

    # Initialize the character to be used for unnamed chains in the
    # structure built from the mmCIF file.
    UNKNOWN_CHAIN_ID = ""

    # Initialize the character to be used for unnamed segments in the
    # structure built from the mmCIF file.
    UNKNOWN_SEG_ID = ""

    # Initialize the number to be used for undefined residue numbers
    # in the structure built from the mmCIF file.
    UNKNOWN_RES_SEQ = 1

    # Initialize the character to be used for undefined atom types in
    # the structure build from the mmCIF file.
    UNKNOWN_TYPE_SYMBOL = ""


    ###################################################################


    def _format_token(self,
                      token):
        """Format string tokens parsed mmCIF file.

        Parameters
        ----------
        token : ``str``
            The token as a string.

        Returns
        -------
        token : ``str``, ``int``, ``float``
            The token as a string, integer, or floating-point
            number.
        """

        # If there is an underscore in the token (as it happens for the
        # '_struct.conn.ptnr1_symmetry' and _struct.conn.ptnr2_symmetry'
        # values)
        if "_" in token:

            # Just return it as it is to avoid it being converted to
            # an integer.
            return token

        #-------------------------------------------------------------#

        # Try to convert the token to an integer.
        try:

            # If successful, return the token as an integer.
            return int(token)

        #-------------------------------------------------------------#

        # If the token could not be converted to an integer
        except ValueError:

            # Try to convert the token to a floating-point number.
            try:

                # If successful, return the token as a floating-point
                # number.
                return float(token)

            #---------------------------------------------------------#

            # If the token could not be converted to a floating-point
            # number.
            except ValueError:

                # If the character is one of those identifying an
                # unassigned value.
                if token in _defaults.MMCIF_UNASSIGNED_CHARACTERS:

                    # Convert the token to an empty string and
                    # return it.
                    return ""
            
                # Just return the token.
                return token


    def _yield_tokens(self,
                      file):
        """Parse all lines of a mmCIF file and yield their constituent
        tokens as strings, one by one.

        Parameters
        ----------
        file : ``str``
            The mmCIF file to be parsed.

        Yields
        ------
        ``str``
            A string representing a token.
        """

        # Open the mmCIF file.
        with open(file, "r") as f:

            # Set a flag for whether we encountered a multi-line quoted
            # string and initialize it to False.
            is_multiline = False

            # Set a variable to store the index of the last multi-line
            # multi-line encountered (for logging purposes).
            multiline_start = None

            # Set an empty list that will store the content of the
            # current multi-line quoted string while parsing it.
            multiline_content = []

            #---------------------------------------------------------#

            # For each line and corresponding index.
            for line_ix, line in enumerate(f):

                #-----------------------------------------------------#

                # If the line is a comment line
                if line.startswith(_defaults.MMCIF_COMMENT_CHARACTER):

                    # Go to the next line.
                    continue

                #-----------------------------------------------------#

                # If the line is the beginning of a data block
                elif line.startswith(_defaults.MMCIF_DATA_DIRECTIVE):

                    # Get the data block's name.
                    block_name = line.split("_")[1].rstrip()

                    # Yield first the directive, then the block's
                    # name (as separate tokens).
                    yield _defaults.MMCIF_DATA_DIRECTIVE
                    yield block_name

                    # Go to the next line ('yield' does not
                    # work exactly as 'return', so we would keep
                    # going if we did not put 'continue' here).
                    continue

                #-----------------------------------------------------#

                # If the line is the beginning of a multi-line quoted
                # string
                elif line.startswith(\
                    _defaults.MMCIF_MULTILINE_CHARACTER):

                    # Strip the line of any trailing white spaces or
                    # newline characters.
                    line = \
                        line.rstrip().lstrip(\
                            _defaults.MMCIF_MULTILINE_CHARACTER)

                    #-------------------------------------------------#

                    # If the line is otherwise empty (= end of a
                    # multi-line quoted string)
                    if not line:

                        # A multi-line has ended. Therefore, yield its
                        # content, separated by a white space.
                        yield " ".join(multiline_content)

                        # Reset the content of the multi-line section
                        # to be used for the next multi-line section
                        # encountered.
                        multiline_content = []

                        # Set the index storing the start of a multi-
                        # line to None.
                        multiline_start = None

                        # Set the flag specifying whether we are inside
                        # a multi-line section to False.
                        is_multiline = False

                        # Go to the next line.
                        continue

                    #-------------------------------------------------#

                    # If the line is not empty
                    else:

                        # If the multi-line starting index is not
                        # defined (= we are at the beginning of a
                        # multi-line quoted string)
                        if multiline_start is None:

                            # Turn on the flag indicating that we are
                            # parsing a multi-line section.
                            is_multiline = True

                            # Save the current index as the start of
                            # the multi-line section.
                            multiline_start = line_ix

                #-----------------------------------------------------#

                # If we are inside a multi-line section
                if is_multiline:

                    # Add the content of the current line to the
                    # content of the multi-line section.
                    multiline_content.append(line.rstrip())

                #-----------------------------------------------------#

                # Otherwise
                else:

                    # Yield the line's content, one token at a time.
                    yield from shlex.split(line)

            #---------------------------------------------------------#

            # If a multi-line was never properly closed
            if is_multiline:

                # Raise an error.
                errstr = \
                    "Missing semicolon in multi-line quoted string " \
                    f"starting in line {multiline_start}."
                raise ValueError(errstr)


    def _parse_as_dict(self,
                       file):
        """Parse a mmCIF file and return a nested dictionary with the
        information in it. The keys of the outer dictionary will
        represent the names of the categories encountered in the
        file (for instance, '_atom_site'), while the keys of the
        inner dictionaries will represent the data items collected
        for each category (for instance, 'label_atom_id' for the
        '_atom_site' category).

        Parameters
        ----------
        file : ``str``
            The mmCIF file.

        Returns
        -------
        data : ``dict``
            The dictionary containing information parsed from the file.
        """

        # Initialize an empty dictionary to store the data parsed from
        # the file.
        data = {}

        #-------------------------------------------------------------#

        # Initialize the current category we are parsing to None.
        current_category = None

        # Initialize the current data item we are parsing to None.
        current_item = None

        #-------------------------------------------------------------#

        # Set a flag for whehter we are entering a 'loop_' section,
        # and initialize it to False.
        is_loop = False

        # Initialize a variable storing the index of the current line
        # containing values parsed from a 'loop_' section to 0.
        current_loop_line_index = 0

        # Initialize a variable storing the index of the latest item
        # found in a 'loop_' section to 0.
        current_loop_item_index = 0

        # Initialize an empty list that will store the data items
        # found in a 'loop_' section.
        current_loop_items = []

        #-------------------------------------------------------------#

        # For each token parsed from the file
        for token in self._yield_tokens(file = file):

            #---------------------------------------------------------#

            # If the token marks the start of a 'loop_' section
            # (the directives are case-insensitive, so we make sure
            # that we capture the token even if it contains capital
            # letters)
            if token.lower() == _defaults.MMCIF_LOOP_DIRECTIVE:

                # Set the 'loop_' flag to True.
                is_loop = True

                # Reset the index keeping track of the index of the
                # last line parsed from the current 'loop_' section.
                current_loop_line_index = 0

                # Reset the index keeping track of the items found
                # in the current 'loop_' section.
                current_loop_item_index = 0

                # Reset the list of items parsed from a 'loop_'
                # section.
                current_loop_items = []

                # Go to the next token.
                continue

            #---------------------------------------------------------#

            # If we have entered a 'loop_' section
            if is_loop:

                # Set an expression to evaluate whether the token
                # is a data item. The conditions that must be met are:
                # - The token starts with an underscore.
                # - AND we are at the start of the portion of the 
                #       'loop_' section containig data items
                #       OR the current line of the 'loop_' section is
                #         a line containing data items.
                is_data_item = \
                    token.startswith("_") and \
                    (current_loop_item_index == 0 or \
                        current_loop_line_index % \
                        current_loop_item_index == 0)

                #-----------------------------------------------------#

                # If the token is a data item
                if is_data_item:

                    # Get the name of the category the 'loop_' section
                    # belongs to and the data item.
                    category, item = token.split(".")

                    # If there is no category with such a name in the
                    # data dictionary
                    if category not in data:

                        # Add the category to the data dictionary.
                        data[category] = {}

                    # If we have already parsed lines containing values
                    # for the data items in this 'loop_' section
                    # (remember that we have already ascertained that
                    # the current token is a data item)
                    if current_loop_line_index > 0:

                        # The 'loop_' section has ended, and the
                        # current data item does not belong to it.
                        is_loop = False

                    # Otherwise, we have encountered a new data item
                    # belonging to the current 'loop_' section.
                    else:

                        # Add the item to the final dictionary.
                        data[category][item] = []

                        # Add the current item to the list of items
                        # found for the current 'loop_' section.
                        current_loop_items.append(item)

                        # Update the index of the last item found in
                        # the 'loop_' section.
                        current_loop_item_index += 1

                        # Go to the next line.
                        continue

                #-----------------------------------------------------#

                # Otherwise (= we are parsing the 'loop_' section's
                # values)
                else:

                    # Find the data item the value belongs to.
                    item = \
                        current_loop_items[\
                            current_loop_line_index % \
                            current_loop_item_index]

                    # Append the current token to the values associated
                    # with that data item item.
                    data[category][item].append(\
                        self._format_token(token = token))

                    # Update the value of the index keeping track of
                    # the number of lines containing values that we
                    # have parsed.
                    current_loop_line_index += 1

                    # Go to the next line.
                    continue

            #---------------------------------------------------------#

            # The following code gets executed when we are not in a
            # 'loop_' section.

            # If the current item is unspecified (= the current token
            # is a data item)
            if current_item is None:

                # If the token starts with an underscore
                if token.startswith("_"):

                    # Get the name of the category and of the data
                    # item.
                    category, item = token.split(".")

                    # If the category is not in the final dictionary
                    # yet.
                    if category not in data:

                        # Add it to the final dictionary.
                        data[category] = {}

                    # Update the variable keeping track of the
                    # category we are currently parsing.
                    current_category = category

                    # Update the variable keeping track of the data
                    # item we are currently parsing.
                    current_item = item

            #---------------------------------------------------------#

            # Otherwise (= the current token is a value)
            else:

                # Store the current token as the value associated
                # with the current data item.
                data[current_category][current_item] = \
                    [self._format_token(token = token)]

                # Reset the variable keeping track of the category we
                # are currently parsing.
                current_section = None

                # Reset the variable keeping track of the data item we
                # are currently parsing.
                current_key = None

        #-------------------------------------------------------------#

        # Return the dictionary.
        return data


    def _check_category(self,
                        data,
                        category_name,
                        strict):
        """Check whether the values reported for each data item in a
        given data category comply with the mmCIF format
        specifications.

        Parameters
        ----------
        data : ``dict``
            A dictionary containing the data parsed from the mmCIF
            file.

        category_name : ``str``
            The name of the data category to be checked.

        strict : ``bool``
            Whether we are running in 'strict' mode or not.
            
            Running in 'strict' mode means that any violation of the
            format specifications raises an error, while, in non-strict
            mode, most violations only trigger a warning.
        """

        # Inform the user we are checking the category.
        infostr = \
            f"Now checking the values in the '{category_name}' " \
            "category..."
        logger.info(infostr)

        #-------------------------------------------------------------#

        # Get the string representing the supported "unassigned"
        # characters (to be used in log messages).
        unassigned_str = \
            ", ".join([f"'{c}'" for c \
                       in _defaults.MMCIF_UNASSIGNED_CHARACTERS])

        #-------------------------------------------------------------#

        # Get the category's data.
        category = data[category_name]

        #-------------------------------------------------------------#

        # For each data item (of those to be checked)
        for item_name in _defaults.MMCIF_CATEGORIES[category_name]:

            #---------------------------------------------------------#

            # Get:
            # - The item's data type.
            # - The controlled vocabulary associated with the item,
            #   if any.
            # - Whether the item must be present in the category.
            # - The default value for the item, if it is undefined.
            item_type, _, item_voc, item_mandatory, item_default = \
                _defaults.MMCIF_CATEGORIES[category_name][item_name]

            #---------------------------------------------------------#

            # If the data item is not present in the category
            if item_name not in category:

                # If the item is mandatory
                if item_mandatory:

                    # Raise an error.
                    errstr = \
                        f"Missing item '{item_name}' in category " \
                        f"'{category_name}'."
                    raise ValueError(errstr)

                # Otherwise
                else:

                    # Go to the next field.
                    continue

            #---------------------------------------------------------#

            # Get the item's values .
            item_vals = category[item_name]

            # Evaluate whether any of the values is of an invalid type
            # and it is not an empty string.
            #
            # We check whether the values are empty strings and not
            # '?' or '.' because those have already been converted to
            # empty strings while parsing the mmCIF file.
            is_any_invalid = \
                any(not isinstance(v, item_type) and v != "" \
                    for v in item_vals)

            # If any invalid values were found
            if is_any_invalid:

                # Raise an error (we raise regardless of whether we
                # are in 'strict' mode or not because, otherwise,
                # we would encounter an error when trying to write
                # out the structure with invalid data types).
                errstr = \
                    "Invalid data type found in item " \
                    f"'{category_name}.{item_name}'. " \
                    f"All values must be either of type " \
                    f"'{item_type}' or any of the characters " \
                    "representing unassigned values " \
                    f"({unassigned_str})."
                raise TypeError(errstr)

            #---------------------------------------------------------#

            # If there is a default value associated with the item
            # (we already made sure the item is not mandatory).
            if item_default is not None:

                # Substitute all empty strings found in the values
                # with the default one before proceeding with the
                # checks.
                item_vals = \
                    [val if val != "" else item_default \
                     for val in item_vals]

                # Warn the user that the substitution has happened
                warnstr = \
                    f"All undefined values in the item " \
                    f"'{category_name}.{item_name}' were replaced " \
                    f"with the default value '{item_default}'."
                logger.warning(warnstr)

            #---------------------------------------------------------#

            # If the item does not have an associated controlled
            # vocabulary.
            if item_voc is None:

                # Go to the next field.
                continue

            #---------------------------------------------------------#

            # If the item has an associated controlled vocabulary
            elif isinstance(item_voc, set):

                # Evaluate whether any of the values does not respect
                # the controlled vocabulary and it is not an empty
                # string.
                #
                # We check whether the values are empty strings and not
                # '?' or '.' because those have already been
                # converted to empty strings while parsing the mmCIF
                # file.
                is_any_not_voc = \
                    any(v not in item_voc and v != "" \
                        for v in item_vals)

                # If any of the values does not appear in the
                # controlled vocabulary
                if is_any_not_voc:

                    # Set a string containing the items not belonging
                    # to the controlled vocabulary.
                    not_voc = \
                        ", ".join(set(\
                            [f"'{v}'" for v in item_vals \
                             if v not in item_voc and v != ""]))

                    # Set a string representing the controlled
                    # vocabulary.
                    item_voc_str = \
                        ", ".join([f"'{i}'" for i in item_voc])

                    # Set the warning/error string.
                    errstr = \
                        "Invalid value found in the item " \
                        f"'{category_name}.{item_name}': {not_voc}. " \
                        "Values should either belong to the " \
                        "controlled vocabulary associated with the " \
                        f"item ({item_voc_str}) or be one of the " \
                        "characters representing unassigned " \
                        f"values ({unassigned_str})."

                    # If we are in 'strict' mode
                    if strict:

                        # Raise an error.
                        raise ValueError(errstr)

                    # Otherwise
                    else:

                        # Warn the user.
                        logger.warning(errstr)

                        # Continue.
                        continue

            #---------------------------------------------------------#

            # If the item is a pointer to another item
            elif isinstance(item_voc, tuple):

                # Get the category and the item the current item
                # points to.
                pointed_category, pointed_item = item_voc
                
                # If the category it points to is not in the data
                if pointed_category not in data:

                    # Set the warning/error string
                    errstr = \
                        f"No category '{pointed_category}' was " \
                        "found. This category should be present " \
                        "when the item " \
                        f"'{category_name}.{item_name}' is " \
                        f"present since the item " \
                        f"'{category_name}.{item_name}' is a " \
                        f"pointer to the item " \
                        f"'{pointed_category}.{pointed_item}'."

                    # If we are in 'strict' mode
                    if strict:

                        # Raise an error.
                        raise ValueError(errstr)

                    # Otherwise
                    else:

                        # Warn the user.
                        logger.warning(errstr)

                        # Continue
                        continue

                # If the item it points to is not in the data
                if pointed_item not in data[pointed_category]:

                    # Set the warning/error string.
                    errstr = \
                        f"No item '{pointed_item}' was found " \
                        f"in category '{pointed_category}'. " \
                        "This item should be present when " \
                        f"the item '{category_name}.{item_name}' " \
                        "is present since the item " \
                        f"'{category_name}.{item_name}' is a " \
                        f"pointer to the item " \
                        f"'{pointed_category}.{pointed_item}'."

                    # If we are in 'strict' mode
                    if strict:

                        # Raise an error.
                        raise ValueError(errstr)

                    # Otherwise
                    else:

                        # Warn the user.
                        logger.warning(errstr)

                        # Continue.
                        continue

                # Get the controlled vocabulary.
                item_voc = data[pointed_category][pointed_item]

                # Evaluate whether any of the values does not
                # respect the controlled vocabulary and it is not
                # an empty string.
                #
                # We check whether the values are empty strings and
                # not '?' or '.' because those have already been
                # converted to empty strings while parsing the mmCIF
                # file.
                is_any_not_voc = \
                    any(v not in item_voc and v != "" \
                        for v in item_vals)

                # If any of the values does not belong to the
                # vocabulary
                if is_any_not_voc:

                    # Set a string containing the items not belonging
                    # to the controlled vocabulary.
                    not_voc = \
                        ", ".join(set(\
                            [f"'{v}'" for v in item_vals \
                             if v not in item_voc and v != ""]))

                    # Set a string containing the controlled
                    # vocabulary.
                    item_voc_str = \
                        ", ".join([f"'{i}'" for i in item_voc])

                    # Set the warning/error string
                    errstr = \
                        "Invalid value found in item " \
                        f"'{category_name}.{item_name}': {not_voc}. " \
                        "Values should either belong to the " \
                        "controlled vocabulary associated with the " \
                        f"item ({item_voc_str}) or be one of the " \
                        "characters representing unassigned " \
                        f"values ({unassigned_str})."

                    # If we are in 'strict' mode
                    if strict:

                        # Raise an error.
                        raise ValueError(errstr)

                    # Otherwise
                    else:

                        # Warn the user.
                        logger.warning(errstr)

            #---------------------------------------------------------#

            # If the controlled vocabulary is a list
            elif isinstance(item_voc, list):

                # The list contains the boundaries of the numerical
                # interval of values the field can contain.
                min_val, max_val = item_voc

                # Evaluate whether any of the values falls outside
                # of the boundaries.
                is_any_outide = \
                    any((isinstance(val, str) \
                            and val != "") or \
                        (not isinstance(val, str) \
                            and not min_val <= val <= max_val) \
                        for val in item_vals)

                # If there are values outside of the boundaries
                if is_any_outide:

                    # Set the warning/error string.
                    errstr = \
                        "Invalid value found in item " \
                        f"'{category_name}.{item_name}'. " \
                        f"Values should be between {min_val} and " \
                        f"{max_val}."

                    # If we are in 'strict' mode
                    if strict:

                        # Raise an error.
                        raise ValueError(errstr)

                    # Otherwise
                    else:

                        # Warn the user.
                        logger.warning(errstr)

            #---------------------------------------------------------#

            # Update the items' values.
            category[item_name] = item_vals

        #-------------------------------------------------------------#

        # Inform the user we have finished checking the category.
        infostr = \
            "Done checking the values in the category " \
            f"'{category_name}'."
        logger.info(infostr)

        # Return the category
        return category


    def _get_atom_data(self,
                       data,
                       use_auth_chain_ids,
                       use_auth_res_seqs,
                       use_auth_res_names,
                       use_auth_atom_names,
                       only_alt_loc,
                       strict):
        """Get a the atomic coordinates and associated data from
        the '_atom_site' data category in a mmCIF file.

        Parameters
        ----------
        data : ``dict``
            A dictionary containing the data parsed from the file.

        use_auth_chain_ids : ``bool``
            Use the authors-defined chain IDs instead of the
            automatically assigned ones, if available.

        use_auth_res_seqs : ``bool``
            Use the authors-defined residue numbers instead of the
            automatically assigned ones, if available.

        use_auth_res_names : ``bool``
            Use the authors-defined residue names instead of the
            automatically assigned ones, if available.

        use_auth_atom_names : ``bool``
            Use the authors-defined atom names instead of the
            automatically assigned ones, if available.

        only_alt_loc : ``str``, ``"A"``
            For atoms having multiple alternate locations, keep only
            those with ``only_alt_loc`` alternate location.

        strict : ``bool``
            Whether to run in strict mode.

            Running in strict mode means that:

            - An error will be raised if any of the categories or
              data items to which some '_atom_site' items point to are
              missing from the file.

            - An error will be raised if any of the values of a data
              item in the '_atom_site' category that must respect a
              controlled vocabulary (according to the mmCIF format
              specifications) does not respect the vocabulary.

            - An error will be raised if any of the values of a 
              numerical data item in the '_atom_site' category that
              must be within certain boundaries (according to the
              mmCIF format specifications) falls outside the
              boundaries.

            If ``False``, only a warning will be issued in the
            previous scenarios.

            An error will be raised if any of the values of a data
            item in the '_atom_site' category is of the wrong type
            regardless of whether we are in 'strict' mode or not.

        Returns
        -------
        atom_data : ``dict``
            A dictionary containing the atomic coordinates and
            associated data.
        """

        # Set the name of the category that contains the atomic
        # coordinates.
        category_name = "_atom_site"

        #-------------------------------------------------------------#

        # If the category is not in the data
        if category_name not in data:

            # Raise an error
            errstr = \
                f"No atomic coordinates found (missing " \
                f"category '{category_name}')."
            raise ValueError(errstr)

        #-------------------------------------------------------------#

        # Check the data in the category.
        category = \
            self._check_category(data = data,
                                 category_name = "_atom_site",
                                 strict = strict)

        #-------------------------------------------------------------#

        # Initialize an empty dictionary to store the non-atom-level
        # attributes.
        non_atom_attrs = {}

        # Create a copy of the category that will contain all atom-
        # level attributes. This will be modified in place.
        atom_attrs = copy.deepcopy(category)

        # Initialize an empty dictionary to store the data that will
        # be used to build the 'Structure'.
        atom_data = {}

        # Initialize an empty dictionary to store, for each model,
        # the identifiers corresponding to a specific chain, residue,
        # and atom (to correctly identify bonds later on).
        models2ids = {}

        #-------------------------------------------------------------#

        # Set a variable to store the current model.
        current_model = None

        # Set a variable to store the current chain ID.
        current_chain_id = None

        # Set a variable to store the current residue.
        current_res = None

        #-------------------------------------------------------------#

        # Get the atoms' serials (in mmCIF files, they do not restart
        # for each model).
        non_atom_attrs["id"] = atom_attrs.pop("id")

        #-------------------------------------------------------------#

        # Get the headers of the PDB-like sections.
        non_atom_attrs["group_PDB"] = atom_attrs.pop("group_PDB")

        #-------------------------------------------------------------#

        # There are no segment IDs in mmCIF files, so just use the
        # "unknown" character to identify the only segment present
        # in chains.
        non_atom_attrs["seg_id"] = \
            [self.UNKNOWN_SEG_ID] * len(non_atom_attrs["id"])

        #-------------------------------------------------------------#

        # Get the atoms' insertion codes.
        non_atom_attrs["pdbx_PDB_ins_code"] = \
            atom_attrs.pop("pdbx_PDB_ins_code")

        #-------------------------------------------------------------#

        # Set the name of the item containing the models.
        item_models = "pdbx_PDB_model_num"

        # If there is a category containing the models
        if item_models in category:

            # Get the models and save them the data dictionary.
            non_atom_attrs[item_models] = atom_attrs.pop(item_models)

        # If no category containing the models was found
        else:

            # Warn the user.
            warnstr = \
                f"{category_name}: no models' numbers found " \
                f"(missing item '{category_name}.{item_models}'). " \
                "All atoms will be considered as belonging to one " \
                f"model, numbered {self.UNKNOWN_MODEL}."
            logger.warning(warnstr)

            # Set it to a list of "unknown" characters and save
            # it in the data dictionary.
            non_atom_attrs[item_models] = \
                [self.UNKNOWN_MODEL] * len(data["id"])

        #-------------------------------------------------------------#

        # Get the automatically assigned chain IDs.
        non_atom_attrs["label_asym_id"] = \
            atom_attrs.pop("label_asym_id")

        # If the user requested the authors-defined IDs
        if use_auth_chain_ids:

            # If there are authors-defined IDs
            if "auth_asym_id" in category:

                # Use them.
                non_atom_attrs["label_asym_id"] = \
                    atom_attrs.pop("auth_asym_id")

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested " \
                    "to use the chain IDs defined by the " \
                    "authors, but none were found (missing item " \
                    f"'{category_name}.auth_asym_id'). " \
                    "The values associated with the item " \
                    f"'{category_name}.label_asym_id' will be " \
                    "used instead."
                logger.warning(warnstr)

        #-------------------------------------------------------------#

        # Get the automatically assigned residues' numbers.
        non_atom_attrs["label_seq_id"] = \
            atom_attrs.pop("label_seq_id")

        # If the user requested the authors-defined numbers
        if use_auth_res_seqs:

            # If there are authors-defined residues' numbers
            if "auth_seq_id" in category:

                # Use them.
                non_atom_attrs["label_seq_id"] = \
                    atom_attrs.pop("auth_seq_id")

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested " \
                    "to use the residues' numbers defined by " \
                    "the authors, but none were found (missing item " \
                    f"'{category_name}.auth_seq_id'). The values " \
                    "associated with the item " \
                    f"'{category_name}.label_seq_id' will be " \
                    "used instead."
                logger.warning(warnstr)

        #-------------------------------------------------------------#

        # Get the names of the residues.
        non_atom_attrs["label_comp_id"] = \
            atom_attrs.pop("label_comp_id")

        # If the user requested the authors-defined names
        if use_auth_res_names:
            
            # If there are authors-defined names
            if "auth_comp_id" in category:

                # Use them.
                non_atom_attrs["label_comp_id"] = \
                    atom_attrs.pop("auth_comp_id")

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested to " \
                    "use the residues' names defined by the " \
                    "authors, but none were found (missing item " \
                    f"'{category_name}.auth_comp_id'). The values " \
                    "associated with the item " \
                    f"'{category_name}.label_comp_id' will be " \
                    "used instead."
                logger.warning(warnstr)

        #-------------------------------------------------------------#

        # If the user requested the authors-defined names
        if use_auth_atom_names:

            # If there are authors-defined names
            if "auth_atom_id" in category:

                # Use them.
                atom_attrs["label_atom_id"] = \
                    atom_attrs.pop("auth_atom_id")

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested to use " \
                    "the atoms' names defined by the authors, but " \
                    "none were found (missing item " \
                    f"'{category_name}.auth_atom_id'). The values " \
                    "associated with the item  " \
                    f"'{category_name}.label_atom_id' will be used " \
                    "instead."
                logger.warning(warnstr)

        #-------------------------------------------------------------#

        # Set the name of the item containing the symbols for the
        # atoms' types.
        item_types_symbols = "type_symbol"

        # If there are symbols available
        if not item_types_symbols in category:

            # Warn the user.
            warnstr = \
                f"{category_name}: no symbols representing the " \
                "atoms' types were found (missing item " \
                f"'{category_name}.{item_types_symbols}'). " \
                "All atom types will be represented with the symbol " \
                f"'{self.UNKNOWN_TYPE_SYMBOL}'."
            logger.warning(warnstr)

            # Set it to a list of "unknown" characters.
            atom_attrs["type_symbol"] = \
                [self.UNKNOWN_TYPE_SYMBOL] * len(non_atom_attrs["id"])

        #-------------------------------------------------------------#

        # For each atom
        for atom_ix in range(len(non_atom_attrs["id"])):

            #---------------------------------------------------------#

            # Get the model the atom belongs to.
            model = non_atom_attrs["pdbx_PDB_model_num"][atom_ix]

            #---------------------------------------------------------#

            # Get the header.
            header = non_atom_attrs["group_PDB"][atom_ix]

            # If the header is HETATM
            if header == "HETATM":

                # Set the heteroresidue flag to True.
                is_het = True

            # If the header is ATOM
            elif header == "ATOM":

                # Set the heteroresidue flag to False.
                is_het = False

            #---------------------------------------------------------#

            # Get the ID of the chain the atom belongs to.
            chain_id = non_atom_attrs["label_asym_id"][atom_ix]

            #---------------------------------------------------------#

            # Get the ID of the segment the atom belongs to.
            seg_id = non_atom_attrs["seg_id"][atom_ix]

            #---------------------------------------------------------#

            # Get the original sequence number of the residue the
            # atom belongs to.
            res_seq_orig = non_atom_attrs["label_seq_id"][atom_ix]

            #---------------------------------------------------------#

            # Get the insertion code of the residue the atom belongs
            # to.
            res_i_code = non_atom_attrs["pdbx_PDB_ins_code"][atom_ix]

            #---------------------------------------------------------#

            # Get the name of the residue the atom belongs to.
            res_name = non_atom_attrs["label_comp_id"][atom_ix]

            #---------------------------------------------------------#

            # If the model is in the dictionary
            if model not in atom_data:

                # Add it.
                atom_data[model] = \
                    {"_items" : {},
                     "_attributes" : {}}

            #---------------------------------------------------------#

            # If the chain ID is not in the dictionary
            if chain_id not in atom_data[model]["_items"]:

                # Add it.
                atom_data[model]["_items"][chain_id] = \
                    {"_items" : {},
                     "_attributes" : {}}

            #---------------------------------------------------------#

            # If the segment ID is not in the dictionary
            if seg_id not in atom_data[model]["_items"][\
                chain_id]["_items"]:

                # Add it.
                atom_data[model]["_items"][chain_id]["_items"][\
                    seg_id] = \
                        {"_items" : {},
                         "_attributes" : {}}

            #---------------------------------------------------------#

            # If the residue's number is undefined
            if not res_seq_orig and res_seq_orig != 0:

                # If we are in a different residue compared to
                # the one the previous atom belongs to, and that
                # residue originally had an undefined sequence number,
                # too
                if not current_res[3] and current_res[2] != res_name:

                    # Set the current residue's number to
                    # the previous residue's number, plus 1.
                    res_seq = current_res[0] + 1

                # Otherwise
                else:

                    # Set it the "unknown" number.
                    res_seq = self.UNKNOWN_RES_SEQ

            # Otherwise
            else:

                # Set it to the residue's original number.
                res_seq = res_seq_orig

            #---------------------------------------------------------#
            
            # If the residue's number is not in the dictionary
            if (res_seq, res_i_code, res_name) not in \
                atom_data[model]["_items"][chain_id]["_items"][\
                    seg_id]["_items"]:

                # If the residue's original number was undefined
                if not res_seq_orig and res_seq_orig != 0:

                    # Warn the user.
                    warnstr = \
                        f"{category_name}: no sequence number " \
                        f"was found for residue '{res_name}' in " \
                        f"chain '{chain_id}', model {model}. " \
                        f"The number {res_seq} was assigned to it."
                    logger.warning(warnstr)

                # Add the residue to the dictionary.
                atom_data[model]["_items"][chain_id]["_items"][\
                    seg_id]["_items"][\
                        (res_seq, res_i_code, res_name)] = \
                            {"_items" : {},
                             "_attributes" : {"is_het" : is_het}}

            #---------------------------------------------------------#

            # Initialize an empty dictionary to store the current
            # atom's attributes.
            atom_attrs_dict = {}

            #---------------------------------------------------------#

            # Set a flag for whether we should skip the atom.
            skip_atom = False

            #---------------------------------------------------------#

            # For each item and associated values
            for item, values in atom_attrs.items():

                # Get the value for the current atom.
                value = values[atom_ix]

                # If the current item is the atom's alternate
                # location
                if item == "label_alt_id":

                    # If the alternate location is not the one we want
                    # to keep
                    if value not in ("", only_alt_loc):

                        # We will skip the atom.
                        skip_atom = True

                        # Exit the loop.
                        break

                    # Otherwise
                    else:

                        # Set the alternate location for the
                        # atom to an empty string.
                        atom_attrs_dict[item] = ""

                # Otherwise
                else:

                    # Add the item and the corresponding value to
                    # the dictionary.
                    atom_attrs_dict[item] = value

            #---------------------------------------------------------#

            # If we should skip the atom
            if skip_atom:

                # Skip it and go to the next atom.
                continue

            #---------------------------------------------------------#

            # If we are now in another model with respect to the one
            # the previous atom belongs to
            if model != current_model:

                # Reset the count for the atoms' serial numbers.
                serial = 1

            #---------------------------------------------------------#

            # If we are now in another chain with respect to the one
            # he previous atom belongs to, but we are still in the
            # same model
            if chain_id != current_chain_id and model == current_model:

                # Add 2 to the current atom's serial (to accommodate
                # for TER records when we write out the structure
                # in PDB format).
                serial += 2
                
            #---------------------------------------------------------#

            # Add the data for the current atom to the dictionary.
            atom_data[model]["_items"][chain_id]["_items"][seg_id][\
                "_items"][(res_seq, res_i_code, res_name)]["_items"][\
                    serial] = \
                        {"_attributes" : atom_attrs_dict}

            #---------------------------------------------------------#

            # Get the identifier for the current atom to be used in the
            # dictionary mapping each model to its atoms' unique
            # identifiers.
            atom_id_table = \
                (chain_id,
                 (res_seq_orig, res_i_code, res_name),
                 atom_attrs["label_atom_id"][atom_ix])

            # If we do not have the model in the dictionary mapping
            # each model to the atoms' identifiers
            if model not in models2ids:

                # Create an entry for the model in the dictionary
                # and add the current atom's identifier.
                models2ids[model] = [atom_id_table]

            # Otherwise
            else:

                # Just add the atom's identifier to the dictionary.
                models2ids[model].append(atom_id_table)

            #---------------------------------------------------------#

            # Update the current model.
            current_model = model

            # Update the current chain ID.
            current_chain_id = chain_id

            # Update the current residue.
            current_res = (res_seq, res_i_code, res_name, res_seq_orig)

            # Update the count for the atoms' serial numbers.
            serial += 1

        #-------------------------------------------------------------#

        # Return:
        # - The dictionary containing the atomic coodinates (and
        #   associated data).
        # - The dictionary that will be used to correctly assign
        #   connectivity records to the atoms.
        return atom_data, models2ids
    

    def _get_conect_data(self,
                         data,
                         models2ids,
                         use_auth_chain_ids,
                         use_auth_res_seqs,
                         use_auth_res_names,
                         use_auth_atom_names,
                         strict):
        """Get the connectivity records from the '_struct_conn'
        data category in a mmCIF file.

        Parameters
        ----------
        data : ``dict``
            A dictionary containing the data parsed from the file.

        models2ids : ``dict``
            A dictionary mapping unique atoms' identifiers to the
            serial numbers the atoms have in the different models.

        use_auth_chain_ids : ``bool``
            Use the authors-defined chain IDs instead of the
            automatically assigned ones, if available.

        use_auth_res_seqs : ``bool``
            Use the authors-defined residue numbers instead of the
            automatically assigned ones, if available.

        use_auth_res_names : ``bool``
            Use the authors-defined residue names instead of the
            automatically assigned ones, if available.

        use_auth_atom_names : ``bool``
            Use the authors-defined atom names instead of the
            automatically assigned ones, if available.

        strict : ``bool``
            Whether to run in strict mode, meaning that:

            - An error will be raised if any of the categories or
              items to which some '_struct_conn' items point to are
              missing from the file.

            - An error will be raised if any of the values of a data
              item in the '_struct_conn' category that must respect a
              controlled vocabulary (according to the mmCIF format's
              specifications) does not respect the vocabulary.

            - An error will be raised if any of the values of a 
              numerical data item in the '_struct_conn' category that
              must be within certain boundaries (according to the
              mmCIF format's specifications) falls outside the
              boundaries.

            If ``False``, only a warning will be issued in the
            previous scenarios.

            An error will be raised if any of the values of
            a data item in the '_struct_conn' category is of
            the wrong type, regardless of whether we are
            in 'strict' mode or not.

        Returns
        -------
        conect_data : ``dict``
            A dictionary containing the connectivity data.
        """

        # Initialize an empty dictionary to store the data parsed from
        # the section.
        conect_data = defaultdict(lambda: defaultdict(dict))

        #-------------------------------------------------------------#

        # Get the name of the category that contains the connectivity
        # data.
        category_name = "_struct_conn"

        #-------------------------------------------------------------#

        # If there is no category containing connectivity data
        if category_name not in data:

            # Return an empty dictionary
            return conect_data

        #-------------------------------------------------------------#

        # Check the category.
        category = \
            self._check_category(data = data,
                                 category_name = "_struct_conn",
                                 strict = strict)

        #-------------------------------------------------------------#

        # Get the bonds' types.
        bonds_types = category["conn_type_id"]

        #-------------------------------------------------------------#

        # Get the chain IDs of the atoms involved in the bonds.
        atoms1_chain_ids = category["ptnr1_label_asym_id"]
        atoms2_chain_ids = category["ptnr2_label_asym_id"]

        # If the user requested the authors-defined IDs
        if use_auth_chain_ids:

            # If there are authors-defined IDs for the first atoms
            if "ptnr1_auth_asym_id" in category:

                # Use them.
                atoms1_chain_ids = category["ptnr1_auth_asym_id"]

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested to " \
                    "use the chain IDs defined by the authors, but " \
                    "none were found (missing item " \
                    f"'{category_name}.ptnr1_auth_asym_id'). The " \
                    "values associated with the item " \
                    f"'{category_name}.ptnr1_label_asym_id' will " \
                    "be used instead."
                logger.warning(warnstr)

            # If there are authors-defined IDs for the second atoms
            if "ptnr2_auth_asym_id" in category:

                # Use them.
                atoms2_chain_ids = category["ptnr2_auth_asym_id"]

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested to use " \
                    "the chain IDs defined by the authors, but " \
                    "none were found (missing item " \
                    f"'{category_name}.ptnr2_auth_asym_id'). The " \
                    "values associated with the item " \
                    f"'{category_name}.ptnr2_label_asym_id' will " \
                    "be used instead."
                logger.warning(warnstr)

        #-------------------------------------------------------------#

        # Get the residues' numbers of the atoms involved in the bonds.
        atoms1_res_seqs = category["ptnr1_label_seq_id"]
        atoms2_res_seqs = category["ptnr2_label_seq_id"]

        # If the user requested the authors-defined numbers
        if use_auth_res_seqs:

            # If there are numbers for the first atoms
            if "ptnr1_auth_seq_id" in category:

                # Use them.
                atoms1_res_seqs = category["ptnr1_auth_seq_id"]

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested to use " \
                    "the residues' numbers defined by the authors, " \
                    "but none were found (missing " \
                    f"item '{category_name}.ptnr1_auth_seq_id'). " \
                    "The values associated with the item " \
                    f"'{category_name}.ptnr1_label_seq_id' will " \
                    "be used instead."
                logger.warning(warnstr)

            # If there are numbers for the second atoms
            if "ptnr2_auth_seq_id" in category:

                # Use them.
                atoms2_res_seqs = category["ptnr2_auth_seq_id"]

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested to use " \
                    "the residues' numbers defined by the authors, " \
                    "but none were found (missing " \
                    f"item '{category_name}.ptnr2_auth_seq_id'). " \
                    "The values associated with the item " \
                    f"'{category_name}.ptnr2_label_seq_id' will " \
                    "be used instead."
                logger.warning(warnstr)

        #-------------------------------------------------------------#

        # Get the residues' names of the atoms involved in the bonds.
        atoms1_res_names = category["ptnr1_label_comp_id"]
        atoms2_res_names = category["ptnr2_label_comp_id"]

        # If the user requested the authors-defined names
        if use_auth_res_names:

            # If there are authors-defined names for the first atoms
            if "ptnr1_auth_comp_id" in category:

                # Use them.
                atoms1_res_names = category["ptnr1_auth_comp_id"]

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested to use " \
                    "the residues' names defined by the authors, " \
                    "but none were found (missing " \
                    f"item '{category_name}.ptnr1_auth_comp_id'). " \
                    "The values associated with the item " \
                    f"'{category_name}.ptnr1_label_comp_id' will " \
                    "be used instead."
                logger.warning(warnstr)

            # If there are authors-defined names for the second atoms
            if "ptnr2_auth_comp_id" in category:

                # Use them.
                atoms2_res_names = category["ptnr2_auth_comp_id"]

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested to use " \
                    "the residues' names defined by the authors, " \
                    "but none were found (missing " \
                    f"item '{category_name}.ptnr2_auth_comp_id'). " \
                    "The values associated with the item " \
                    f"'{category_name}.ptnr2_label_comp_id' will " \
                    "be used instead."
                logger.warning(warnstr)

        #-------------------------------------------------------------#

        # Get the names of the atoms involved in the bonds.
        atoms1_names = category["ptnr1_label_atom_id"]
        atoms2_names = category["ptnr2_label_atom_id"]

        # If the user requested the authors-defined names
        if use_auth_atom_names:

            # If there are authors-defined names for the first atoms
            if "ptnr1_auth_atom_id" in category:

                # Use them.
                atoms1_names = category["ptnr1_auth_atom_id"]

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested to use " \
                    "the atoms' names defined by the authors, but " \
                    "none were found (missing item " \
                    f"'{category_name}.ptnr1_auth_atom_id'). " \
                    "The values associated with the item  " \
                    f"'{category_name}.ptnr1_label_atom_id' will " \
                    "be used instead."
                logger.warning(warnstr)

            # If there are authors-defined names for the second atoms
            if "ptnr1_auth_atom_id" in category:

                # Use them.
                atoms2_names = category["ptnr2_auth_atom_id"]

            # Otherwise
            else:

                # Warn the user.
                warnstr = \
                    f"{category_name}: the user requested to use " \
                    "the atoms' names defined by the authors, but " \
                    "none were found (missing item " \
                    f"'{category_name}.ptnr2_auth_atom_id'). " \
                    "The values associated with the item  " \
                    f"'{category_name}.ptnr2_label_atom_id' will " \
                    "be used instead."
                logger.warning(warnstr)

        #-------------------------------------------------------------#

        # If there are insertion codes specified for the first atoms
        # in the bonds, use them. Otherwise, set the insertion codes
        # to empty strings.
        atoms1_icodes = \
            category["pdbx_ptnr1_PDB_ins_code"] \
            if "pdbx_ptnr1_PDB_ins_code" in category \
            else [""] * len(atoms1_names)

        # If there are insertion codes specified for the second atoms
        # in the bonds, use them. Otherwise, set the insertion codes
        # to empty strings.
        atoms2_icodes = \
            category["pdbx_ptnr2_PDB_ins_code"] \
            if "pdbx_ptnr2_PDB_ins_code" in category \
            else [""] * len(atoms1_names)

        #-------------------------------------------------------------#

        # Zip all data that will be used to identify the bonds.
        zipped_data = \
            zip(atoms1_chain_ids, atoms2_chain_ids,
                atoms1_res_seqs, atoms2_res_seqs,
                atoms1_icodes, atoms2_icodes,
                atoms1_res_names, atoms2_res_names,
                atoms1_names, atoms2_names,
                bonds_types)

        #-------------------------------------------------------------#

        # For each piece of data
        for atom_ix, \
            (atom1_chain_id, atom2_chain_id, \
             atom1_res_seq, atom2_res_seq, \
             atom1_icode, atom2_icode, \
             atom1_res_name, atoms2_res_name, \
             atom1_name, atom2_name, \
             bond_type) in enumerate(zipped_data):

            #---------------------------------------------------------#

            # Build the unique identifier for the first atom involved
            # in the bond.
            atom1_id = \
                (atom1_chain_id,
                 (atom1_res_seq, atom1_icode, atom1_res_name),
                 atom1_name)

            # Build the unique identifier for the second atom involved
            # in the bond.
            atom2_id = \
                (atom2_chain_id,
                 (atom2_res_seq, atom2_icode, atoms2_res_name),
                 atom2_name)

            #---------------------------------------------------------#

            # For each model found in the mmCIF file
            for model in models2ids:

                # If both atoms are present in the current model
                if atom1_id in models2ids[model] \
                and atom2_id in models2ids[model]:

                    # Build the ID that will be used in the
                    # connectivity data for the first atom.
                    atom1 = \
                        (atom1_id[0], "", atom1_id[1], atom1_id[2])

                    # Build the ID that will be used in the
                    # connectivity data for the second atom.
                    atom2 = \
                        (atom2_id[0], "", atom2_id[1], atom2_id[2])

                    # Initialize an empty dictionary to store the
                    # current bond's data.
                    bond_data = {}

                    # For each data item in the category
                    for item in category:

                        # If the item is one of those containing
                        # 'ptnr'
                        if "ptnr" in item:

                            # If it is not one of those regarding
                            # the symmetry
                            if item not in \
                                ("ptnr1_symmetry", "ptnr2_symmetry"):

                                # Ignore the item
                                continue

                        # Add the item and its value to the bond's
                        # data.
                        bond_data[item] = category[item][atom_ix]

                    # Add the first atom to the dictionary.
                    conect_data[model][atom1][atom2] = bond_data

                    # Add the second atom to the dictionary.
                    conect_data[model][atom2][atom1] = bond_data

        #-------------------------------------------------------------#

        # If no connectivity data were found
        if not conect_data:
            
            # Return None.
            return None

        # Otherwise
        else:

            # Convert the nested default dictionary to a standard
            # nested dictionary.
            conect_data = \
                {model : \
                    {atom1 : \
                        {atom2 : data for atom2, data \
                         in conect_data[model][atom1].items()} \
                     for atom1 in conect_data[model]} \
                 for model in conect_data}

            # Return the connectivity data.
            return conect_data


    ###################################################################


    def parse(self,
              file,
              only_alt_loc = "A",
              use_auth_chain_ids = True,
              use_auth_res_seqs = True,
              use_auth_res_names = True,
              use_auth_atom_names = True,
              strict = False):
        """Parse a mmCIF file.

        Parameters
        ----------
        file : ``str``
            The mmCIF file to be parsed.

        only_alt_loc : ``str``, ``"A"``
            For atoms having multiple alternate locations, keep only
            those with ``only_alt_loc`` alternate location.

        use_auth_chain_ids : ``bool``, ``True``
            Use the authors-defined chain IDs instead of the
            automatically assigned ones, if available.

        use_auth_res_seqs : ``bool``, ``True``
            Use the authors-defined residue numbers instead of the
            automatically assigned ones, if available.

        use_auth_res_names : ``bool``, ``True``
            Use the authors-defined residue names instead of the
            automatically assigned ones, if available.

        use_auth_atom_names : ``bool``, ``True``
            Use the authors-defined atom names instead of the
            automatically assigned ones, if available.

        strict : ``bool``, ``False``
            Whether to run in strict mode, meaning that:

            - An error will be raised if any of the categories or
              data items to which some '_atom_site' items point to are
              missing from the file.

            - An error will be raised if any of the values of a data
              item in the '_atom_site' category that must respect a
              controlled vocabulary (according to the mmCIF format's
              specifications) does not respect the vocabulary.

            - An error will be raised if any of the values of a 
              numerical data item in the '_atom_site' category that
              must be within certain boundaries (according to the
              mmCIF format's specifications) falls outside the
              boundaries.

            If ``False``, only a warning will be issued in the
            previous scenarios.

            An error will be raised if any of the values of
            a data item in the '_atom_site' category is of
            the wrong type regardless of whether we are
            in 'strict' mode or not.

        Returns
        -------
        struct : :class:`pdbcraft.structure.Structure`
            The parsed structure.
        """

        #-------------------------------------------------------------#

        # Inform the user that we are parsing the file.
        infostr = f"Now parsing mmCIF file '{file}'..."
        logger.info(infostr)

        # Parse the raw data from the file as a dictionary.
        data = self._parse_as_dict(file = file)

        #-------------------------------------------------------------#

        # Get the 'atom_data' dictionary and the table that will be
        # used to correctly identify the bonds, if recorded in the
        # data.
        atom_data, models2ids = \
            self._get_atom_data(\
                data = data,
                use_auth_chain_ids = use_auth_chain_ids,
                use_auth_res_seqs = use_auth_res_seqs,
                use_auth_res_names = use_auth_res_names,
                use_auth_atom_names = use_auth_atom_names,
                only_alt_loc = only_alt_loc,
                strict = strict)

        #-------------------------------------------------------------#

        # If no atomic coordinates were found
        if not atom_data:

            # Raise an error.
            errstr = \
                "Empty mmCIF file: no '_atom_site' category found " \
                f"in '{file}'."
            raise ValueError(errstr)

        #-------------------------------------------------------------#

        # Get the dictionary containing the connectivity data.
        conect_data = \
            self._get_conect_data(\
                data = data,
                models2ids = models2ids,
                use_auth_chain_ids = use_auth_chain_ids,
                use_auth_res_seqs = use_auth_res_seqs,
                use_auth_res_names = use_auth_res_names,
                use_auth_atom_names = use_auth_atom_names,
                strict = strict)

        #-------------------------------------------------------------#

        # Get the name of the file (without the path leading to it,
        # if any is present).
        file_name, _ = os.path.splitext(os.path.basename(file))

        #-------------------------------------------------------------#

        # Create the structure.
        struct = \
            Structure(atom_data = atom_data,
                      conect_data = conect_data,
                      name = file_name)

        # Renumber the atoms since parts of different chains may have
        # been separated in the mmCIF file.
        struct._update_atom_and_conect()

        # Inform the user about the structure's creation.
        infostr = \
            f"The structure '{file_name}' was created from '{file}.'"
        logger.info(infostr)

        #-------------------------------------------------------------#

        # Return the structure.
        return struct
