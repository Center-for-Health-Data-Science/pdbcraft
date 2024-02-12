#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    mmcifwriters.py
#
#    Writers for mmCIF files.
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


# Standard libaray
import itertools
import logging as log
# pdbcraft
from . import _defaults


# Get the module's logger
logger = log.getLogger(__name__)


class MMCIFWriter:

    """
    A writer for mmCIF files.
    """

    #-------------------------- Attributes ---------------------------#


    # The character to be assigned to a value when the value is
    # unknown
    UNKNOWN_CHARACTER = "?"


    #----------------------- Private methods -------------------------#


    def _get_format_strings_atom_site(self,
                                      struct):
        """Get the format strings to be used for data items
        in the '_atom_site' category.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written out.

        Returns
        -------
        fmt_strings : ``dict``
            A dictionary mapping the data items to the format
            strings used for them.
        """

        # Initialize a dictionary containing, for the '_atom_site'
        # data category, the names of the non-atom-level items
        # (plus the 'id' item) mapped to the default width they
        # will have in the output mmCIF file
        i2w = \
            {"group_PDB" : 6,
             "id" : 1,
             "pdbx_PDB_model_num" : 1,
             "label_asym_id" : 1,
             "label_seq_id" : 1,
             "pdbx_PDB_ins_code" : 1,
             "label_comp_id" : 1,
             "auth_seq_id" : 1,
             "auth_comp_id" : 1,
             "auth_asym_id" : 1,
             "auth_atom_id" : 1}

        #-------------------------------------------------------------#

        # Initialize a variable to store the current atom's
        # unique id
        atom_id = 1

        #-------------------------------------------------------------#

        # For each model in the structure and associated data
        for mod, mod_d in struct.atom_data.items():

            # If the string representing the current model
            # is longer than the longest one found so far
            if len(str(mod)) > i2w["pdbx_PDB_model_num"]:

                # Update the width of the corresponding data item
                i2w["pdbx_PDB_model_num"] = len(str(mod))

            #---------------------------------------------------------#

            # For each chain and associated data
            for ch, ch_d in mod_d["_items"].items():

                # If the string representing the current chain
                # is longer than the longest one found so far
                if len(str(ch)) > i2w["label_asym_id"]:

                    # Update the width of the corresponding data
                    # item
                    i2w["label_asym_id"] = len(str(ch))

                    # Update the width of the corresponding 'auth'
                    # data item
                    i2w["auth_asym_id"] = len(str(ch))

                #-----------------------------------------------------#
                
                # For each segment and associated data
                for seg, seg_d in ch_d["_items"].items():

                    # We do nothing here since there is no concept
                    # of 'segment' is mmCIF files

                    #-------------------------------------------------#

                    # For each residue and associated data
                    for res, res_d in seg_d["_items"].items():

                        # Get the residue's number, insertion
                        # code, and name
                        res_seq, i_code, res_name = res

                        #---------------------------------------------#

                        # If the string representing the current
                        # residue's number is longer than the
                        # longest one found so far
                        if len(str(res_seq)) > i2w["label_seq_id"]:

                            # Update the width of the corresponding
                            # data item
                            i2w["label_seq_id"] = len(str(res_seq))

                            # Update the width of the corresponding
                            # 'auth' data item
                            i2w["auth_seq_id"] = len(str(res_seq))

                        #---------------------------------------------#

                        # If the string representing the current
                        # residue's insertion code is longer than
                        # the longest one found so far
                        if len(str(i_code)) > i2w["pdbx_PDB_ins_code"]:

                            # Update the width of the corresponding
                            # data item
                            i2w["pdbx_PDB_ins_code"] = len(str(i_code))
                        
                        #---------------------------------------------#

                        # If the string representing the current
                        # residue's name is longer than the longest
                        # one found so far
                        if len(str(res_name)) > i2w["label_comp_id"]:
                            
                            # Update the width of the corresponding
                            # data item
                            i2w["label_comp_id"] = len(str(res_name))

                            # Update the width of the corresponding
                            # 'auth' data item
                            i2w["auth_comp_id"] = len(str(res_name))

                        #---------------------------------------------#

                        # For each atom and associated data
                        for atom, atom_d in res_d["_items"].items():

                            # Update the variable storing the
                            # atom's unique ID
                            atom_id += 1

                            # For each attribute of the atom
                            for attr, val \
                            in atom_d["_attributes"].items():

                                # If the attribute is not
                                # in the dictionary yet
                                if attr not in i2w:

                                    # Add it, together with the
                                    # length of the string
                                    # representation of the
                                    # associated value. If the
                                    # string is empty, add 1
                                    # instead since in the output
                                    # mmCIF file all empty strings
                                    # will be replaced by the
                                    # 'unknown' character (which
                                    # has a length of 1)
                                    i2w[attr] = \
                                        len(str(val)) \
                                        if len(str(val)) > 0 \
                                        else 1

                                # Otherwise
                                else:

                                    # If the string representing
                                    # the current attribute's value
                                    # is longer than the longest
                                    # one found so far
                                    if len(str(val)) > i2w[attr]:

                                        # Update the width of the
                                        # corresponding data item
                                        i2w[attr] = len(str(val))

                                        # If the attribute is the
                                        # 'label_atom_id' attribute
                                        if attr == "label_atom_id":

                                            # Update the corresponding
                                            # 'auth' item as well
                                            i2w["auth_atom_id"] = \
                                                len(str(val))

        #-------------------------------------------------------------#

        # Update the width of the item containing the atoms'
        # unique IDs
        i2w["id"] = len(str(atom_id))

        #-------------------------------------------------------------#

        # Sort the data items in the dictionary
        i2w = \
            {k: i2w[k]+1 for k \
                in _defaults.MMCIF_CATEGORIES["_atom_site"] \
             if k in i2w}

        #-------------------------------------------------------------#

        # Return the dictionary
        return i2w


    def _get_format_strings_struct_conn(self,
                                        struct):
        """Get the format strings to be used for data items
        in the '_struct_conn' category.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written out.

        Returns
        -------
        fmt_strings : ``dict``
            A dictionary mapping the data items to the format
            strings used for them.
        """

        # Initialize a dictionary containing, for the '_struct_conn'
        # data category, the names of the non-atom-level items
        # (plus the 'id' item) mapped to the default width they
        # will have in the output mmCIF file
        i2w = \
            {"conn_id" : 1,
             "conn_type_id" : 1,
             "ptnr1_label_asym_id" : 1,
             "ptnr1_label_comp_id" : 1,
             "ptnr1_label_seq_id" : 1,
             "ptnr1_label_atom_id" : 1,
             "ptnr2_label_asym_id" : 1,
             "ptnr2_label_comp_id" : 1,
             "ptnr2_label_seq_id" : 1,
             "ptnr2_label_atom_id" : 1,
             "ptnr1_auth_asym_id" : 1,
             "ptnr1_auth_comp_id" : 1,
             "ptnr1_auth_seq_id" : 1,
             "ptnr2_auth_asym_id" : 1,
             "ptnr2_auth_comp_id" : 1,
             "ptnr2_auth_seq_id" : 1,
             "pdbx_ptnr1_PDB_ins_code" : 1,
             "pdbx_ptnr2_PDB_ins_code" : 1,
             "pdbx_ptnr1_label_alt_id" : 1,
             "pdbx_ptnr2_label_alt_id" : 1}

        #-------------------------------------------------------------#

        # For each model in the connectivity data
        for mod in struct.conect_data:

            # For each atom and bonded atoms
            for atom1, atoms_bonded in struct.conect_data[mod].items():

                # Get the chain the atom is on
                len_atom1_ch = len(str(atom1[0]))

                # If the length of the string representation of
                # the current chain is longer than the longest
                # one found so far
                if len_atom1_ch > i2w["ptnr1_label_asym_id"]:

                    # Update the corresponding 'label' and 'auth'
                    # data items
                    i2w["ptnr1_label_asym_id"] = len_atom1_ch
                    i2w["ptnr1_auth_asym_id"] = len_atom1_ch

                #-----------------------------------------------------#

                # Get the sequence number of the residue the atom
                # is on
                len_atom1_res_seq = len(str(atom1[2][0]))

                # If the length of the string representation of
                # the current residue's sequence number is longer
                # than the longest one found so far      
                if len_atom1_res_seq > i2w["ptnr1_label_seq_id"]:

                    # Update the corresponding 'label' and 'auth'
                    # data items
                    i2w["ptnr1_label_seq_id"] = len_atom1_res_seq
                    i2w["ptnr1_auth_seq_id"] = len_atom1_res_seq

                #-----------------------------------------------------#

                # Get the insertion code of the residue the atom
                # is on
                len_atom1_i_code = len(str(atom1[2][1]))

                # If the length of the string representation of
                # the current residue's insertion code is longer
                # than the longest one found so far  
                if len_atom1_i_code > i2w["pdbx_ptnr1_PDB_ins_code"]:

                    # Update the corresponding 'label' data item
                    i2w["pdbx_ptnr1_PDB_ins_code"] = len_atom1_i_code

                #-----------------------------------------------------#

                # Get the name of the residue the atom is on
                len_atom1_res_name = len(str(atom1[2][2]))

                # If the length of the string representation of
                # the current residue's name is longer
                # than the longest one found so far  
                if len_atom1_res_name > i2w["ptnr1_label_comp_id"]:

                    # Update the corresponding 'label' and 'auth'
                    # data items
                    i2w["ptnr1_label_comp_id"] = len_atom1_res_name
                    i2w["ptnr1_auth_comp_id"] = len_atom1_res_name

                #-----------------------------------------------------#

                # Get the atom's name
                len_atom1_label_id = len(str(atom1[3]))

                # If the length of the string representation of
                # the current atom's name is longer
                # than the longest one found so far  
                if len_atom1_label_id > i2w["ptnr1_label_atom_id"]:

                    # Update the corresponding 'label' data item
                    i2w["ptnr1_label_atom_id"] = len_atom1_label_id 

                #-----------------------------------------------------#

                # For each bonded atom and associated information
                # about the bond
                for atom2, bond_attrs in atoms_bonded.items():

                    # Get the chain the atom is on
                    len_atom2_ch = len(str(atom2[0]))

                    # If the length of the string representation of
                    # the current chain is longer than the longest
                    # one found so far
                    if len_atom2_ch > i2w["ptnr2_label_asym_id"]:

                        # Update the corresponding 'label' and 'auth'
                        # data items
                        i2w["ptnr2_label_asym_id"] = len_atom2_ch
                        i2w["ptnr2_auth_asym_id"] = len_atom2_ch

                    #-------------------------------------------------#

                    # Get the sequence number of the residue the atom
                    # is on
                    len_atom2_res_seq = len(str(atom2[2][0]))

                    # If the length of the string representation of
                    # the current residue's sequence number is longer
                    # than the longest one found so far      
                    if len_atom2_res_seq > i2w["ptnr2_label_seq_id"]:

                        # Update the corresponding 'label' and 'auth'
                        # data items
                        i2w["ptnr2_label_seq_id"] = len_atom2_res_seq
                        i2w["ptnr2_auth_seq_id"] = len_atom2_res_seq

                    #-------------------------------------------------#

                    # Get the insertion code of the residue the atom
                    # is on
                    len_atom2_i_code = len(str(atom2[2][1]))

                    # If the length of the string representation of
                    # the current residue's insertion code is longer
                    # than the longest one found so far  
                    if len_atom2_i_code > \
                        i2w["pdbx_ptnr2_PDB_ins_code"]:

                        # Update the corresponding 'label' data item
                        i2w["pdbx_ptnr2_PDB_ins_code"] = \
                            len_atom2_i_code

                    #-------------------------------------------------#

                    # Get the name of the residue the atom is on
                    len_atom2_res_name = len(str(atom2[2][2]))

                    # If the length of the string representation of
                    # the current residue's name is longer
                    # than the longest one found so far  
                    if len_atom2_res_name > i2w["ptnr2_label_comp_id"]:

                        # Update the corresponding 'label' and 'auth'
                        # data items
                        i2w["ptnr2_label_comp_id"] = len_atom2_res_name
                        i2w["ptnr2_auth_comp_id"] = len_atom2_res_name

                    #-------------------------------------------------#

                    # Get the atom's name
                    len_atom2_label_id = len(str(atom2[3]))

                    # If the length of the string representation of
                    # the current atom's name is longer
                    # than the longest one found so far  
                    if len_atom2_label_id > i2w["ptnr2_label_atom_id"]:

                        # Update the corresponding 'label' data item
                        i2w["ptnr2_label_atom_id"] = len_atom2_label_id 

                    #-------------------------------------------------#

                    # For each bond attribute and its value
                    for bond_attr, val in bond_attrs.items():

                        # If the attribute is not in the dictionary
                        # yet
                        if bond_attr not in i2w:

                            # Add it
                            i2w[bond_attr] = \
                                len(str(val)) if len(str(val)) > 0 \
                                else 1

                        # Otherwise
                        else:

                            # If the string representation of the
                            # attribute's value is longer than the
                            # longest one found so far
                            if len(str(val)) > i2w[bond_attr]:

                                # Update it
                                i2w[bond_attr] = len(str(val))

        #-------------------------------------------------------------#

        # Sort the data items in the dictionary
        i2w = \
            {k: i2w[k] for k \
                in _defaults.MMCIF_CATEGORIES["_struct_conn"] \
             if k in i2w}

        #-------------------------------------------------------------#

        # Return the dictionary
        return i2w


    def _get_format_strings(self,
                            struct):
        """Get the format strings to be used for each data item
        in the '_atom_site' and '_struct_conn' data categories
        in the output mmCIF file.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written out.

        Returns
        -------
        fmt_strings : ``dict``
            A dictionary mapping the data items to the format
            strings used for them.
        """

        # Get the format strings for the '_atom_site' category
        i2w_atom = \
            self._get_format_strings_atom_site(struct = struct)

        # Get the format strings for the '_struct_conn' category
        i2w_con = \
            self._get_format_strings_struct_conn(struct = struct)

        #-------------------------------------------------------------#

        # Merge them in a dictionary for easy access
        item2width = \
            {"_atom_site" : i2w_atom,
             "_struct_conn" : i2w_con}

        #-------------------------------------------------------------#

        # Initialize an empty dictionary to store the format strings
        fmt_strings = {category : {} for category in item2width}

        # Set the template format strings for floating point numbers,
        # integers, and strings
        fmt_string_float = "{{:<{item_width}.{item_float_digits}f}}"
        fmt_string_int = "{{:<{item_width}d}}"
        fmt_string_str = "{{:<{item_width}s}}"

        #-------------------------------------------------------------#

        # For each category and associated dictionary of data items/
        # item's widths
        for category, i2w in item2width.items():

            # Get the informaton about each item (such as the value's
            # data type and, if it a floating point number, how many
            # decimal digits should be written out, etc.)
            items_info = _defaults.MMCIF_CATEGORIES[category]

            # For each item'a name and width
            for item_name, item_width in i2w.items():

                # Get the value's data type and the number of decimal
                # digits to be written out, if it is a floating-point
                # number
                item_type, item_float_digits = \
                    items_info[item_name][:2]

                #-----------------------------------------------------#

                # If the value is a floating point number
                if item_type == float:

                    # Set the corresponding format string
                    fmt_string = \
                        fmt_string_float.format(\
                            item_width = item_width,
                            item_float_digits = item_float_digits)

                #-----------------------------------------------------#

                # If the value is an integer
                elif item_type == int:

                    # Set the corresponding format string
                    fmt_string = \
                        fmt_string_int.format(\
                            item_width = item_width)

                #-----------------------------------------------------#

                # If the value is a string
                elif item_type == str:

                    # Set the corresponding format string
                    fmt_string = \
                        fmt_string_str.format(\
                            item_width = item_width)

                #-----------------------------------------------------#

                # Add the string to the dictionary of format strings
                fmt_strings[category][item_name] = fmt_string

        #-------------------------------------------------------------#      

        # Return the format strings
        return fmt_strings


    def _write_items(self,
                     file_handle,
                     items,
                     fmt_strings):
        """Write out the data items.

        Parameters
        ----------
        file_handle : file_handle
            The file handle.

        items : ``dict``
            A dictionary mapping the data items to their values.

        fmt_strings : ``dict``
            A dictionary mapping the data items to the format
            strings used for them.
        """

        # For each item and associated format string
        for item_ix, (item, fmt_str) in enumerate(fmt_strings.items()):

            #---------------------------------------------------------#

            # Get the value of the item for the current atom
            value = \
                items[item] \
                if items[item] is not None \
                and items[item] not in ("", " ") \
                else self.UNKNOWN_CHARACTER

            #---------------------------------------------------------#

            # Try to write out the value using the format string
            try:
                
                file_handle.write(fmt_str.format(value))

            # If an exception arises - the value is unknown,
            # so format strings for floating-point numbers
            # and integer fail since the character used
            # for unknown values is a string
            except ValueError:

                # Modify the format string so that it accepts
                # strings
                fmt_str = fmt_str[:-2] + "s}"

                # Write out the value
                file_handle.write(fmt_str.format(value))

            #---------------------------------------------------------#

            # If we are at the last item
            if item_ix == (len(fmt_strings)-1):

                # Terminate the line
                file_handle.write("\n")

                # Break out of the loop
                break

            #---------------------------------------------------------#

            # Write a whitespace between each item
            file_handle.write(" ")


    def _write_atom_type(self,
                         struct,
                         file_handle):
        """Write the ``Structure``'s unique atim types to the
        '_atom_type' category.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written.

        file_handle : file handle
            The file handle.
        """

        # Get the unique type symbols per model
        unique_type_symbols_per_model = \
            struct._get(action = "get_unique",
                        level = "atoms",
                        selected = \
                            {"atoms" : \
                                {"_attributes" : ["type_symbol"]}},
                        squeeze = "models",
                        elements_type = "attributes")

        #-------------------------------------------------------------#

        # Get the unique type symbols
        unique_type_symbols = \
            sorted(list(set(\
                itertools.chain.from_iterable(\
                    [v for v \
                     in unique_type_symbols_per_model.values()]))))

        #-------------------------------------------------------------#

        # If there are any type symbols
        if unique_type_symbols:

            #---------------------------------------------------------#

            # Write out the 'loop_' directive
            file_handle.write(\
                f"{_defaults.MMCIF_LOOP_DIRECTIVE}\n")

            #---------------------------------------------------------#

            # Write out the name of the data category
            file_handle.write("_atom_type.symbol\n")

            #---------------------------------------------------------#

            # For each type symbol
            for type_symbol in unique_type_symbols:

                # Write it out
                file_handle.write(f"{type_symbol}\n")

            #---------------------------------------------------------#

            # Write out a comment character
            file_handle.write("#\n")


    def _write_atom_site(self,
                         struct,
                         file_handle,
                         fmt_strings):
        """Write the ``Structure``'s atomic coordinates to the
        '_atom_site' category.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written.

        file_handle : file handle
            The file handle.

        fmt_strings : ``dict``
            A dictionary mapping the data items to the format
            strings used for them.
        """

        # Set a counter for the atoms' unique IDs (they are
        # continuous in mmCIF files)
        atom_id = 1

        #-------------------------------------------------------------#

        # Write out the 'loop_' directive
        file_handle.write(\
            f"{_defaults.MMCIF_LOOP_DIRECTIVE}\n")

        #-------------------------------------------------------------#

        # For each item in the '_atom_site' section
        for item in fmt_strings:

            # Write out its name
            file_handle.write(f"_atom_site.{item}\n")

        #-------------------------------------------------------------#

        # For each model and associated data
        for i, (mod, mod_d) in enumerate(struct.atom_data.items()):

            # For each chain and associated data
            for ch, ch_d in mod_d["_items"].items():

                # For each segment and associated data
                for seg, seg_d in ch_d["_items"].items():

                    # For each residue and associated data
                    for (res_seq, res_i_code, res_name), res_d \
                        in seg_d["_items"].items():

                        # Initialize an empty dictionary to store
                        # the values for each item
                        items = {item : None for item in fmt_strings}

                        # Add the model to the dictionary of
                        # values to be used in the current
                        # record
                        items["pdbx_PDB_model_num"] = mod

                        # Add the chain ID to the dictionary of
                        # values to be used in the current record
                        items["label_asym_id"] = ch

                        # Add the chain ID to the dictionary of values
                        # to be used in the current record also for
                        # the 'auth' record
                        items["auth_asym_id"] = ch

                        # Add the residue's number to the
                        # dictionary of values to be used
                        # in the current record
                        items["label_seq_id"] = res_seq

                        # Add the residue's number to the
                        # dictionary of values to be used
                        # in the current record also for
                        # the 'auth' records
                        items["auth_seq_id"] = res_seq

                        # Add the residue's insertion code
                        # to the dictionary of values to
                        # be used in the current record
                        items["pdbx_PDB_ins_code"] = res_i_code

                        # Add the residue's name to the
                        # dictionary of values to be used
                        # in the current record
                        items["label_comp_id"] = res_name

                        # Add the residue's name to the
                        # dictionary of values to be used
                        # in the current record also for
                        # the 'auth' records
                        items["auth_comp_id"] = res_name

                        # If we need to write a HETATM record
                        if res_d["_attributes"]["is_het"]:
                            
                            # Set the corresponding 'group_PDB'
                            # for the current record
                            items["group_PDB"] = "HETATM"

                        # If we need to write an ATOM record
                        else:

                            # Set the corresponding 'group_PDB'
                            # for the current record
                            items["group_PDB"] = "ATOM"

                        # For each atom and associated data
                        for atom, atom_d in res_d["_items"].items():

                            # Set the atom's unique ID in the
                            # dictionary of values to be used
                            # for the current record
                            items["id"] = atom_id

                            # Update the dictionary of values to
                            # be used for the current record with
                            # the atom's attributes (we are
                            # guaranteed that all atoms'
                            # attributes are among the items,
                            # because of how the items'
                            # dictionary is initialized)
                            items.update(\
                                {item : val for item, val \
                                 in atom_d["_attributes"].items() \
                                 if item in items})

                            # Try to update the 'auth' item
                            # corresponding to 'label_atom_id'
                            try:
                                
                                items.update(\
                                    {"auth_atom_id" : \
                                        atom_d["_attributes"][\
                                            "label_atom_id"]})
                            
                            # If the atom has no 'label_atom_id'
                            # attribute
                            except KeyError:
                                
                                # Raise an error
                                errstr = \
                                    "Atoms must have a " \
                                    "'label_atom_id' attribute " \
                                    "to write a valid mmCIF " \
                                    "file."
                                raise ValueError(errstr)

                            # Write out the data items
                            self._write_items(\
                                file_handle = file_handle,
                                items = items,
                                fmt_strings = fmt_strings)

                            # Update the atoms' IDs' counter
                            atom_id += 1

        #-------------------------------------------------------------#

        # Write out a comment line to terminate the section
        file_handle.write("#\n")


    def _write_struct_conn(self,
                           struct,
                           file_handle,
                           fmt_strings):
        """Write the ``Structure``'s connectivity data to the
        '_struct_conn' category.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written.

        file_handle : file handle
            The file handle.

        fmt_strings : ``dict``
            A dictionary mapping the data items to the format
            strings used for them.
        """

        # If the structure does not have associated connectivity
        # data
        if not any(struct.conect_data.values()):

            # Simply return
            return

        #-------------------------------------------------------------#

        # Initialize an empty set to keep track of which
        # bonds have been written out, so that we do not
        # write them out twice
        written = set()

        #-------------------------------------------------------------#

        # Initialize an empty list to store the sorted
        # connectivity records
        sorted_conect = []

        #-------------------------------------------------------------#

        # Write out the 'loop_' directive
        file_handle.write(\
            f"{_defaults.MMCIF_LOOP_DIRECTIVE}\n")

        #-------------------------------------------------------------#

        # For each item
        for item in fmt_strings:

            # Write out its name
            file_handle.write(f"_struct_conn.{item}\n")

        #-------------------------------------------------------------#

        # For each model in the existent connectivity data
        for i, (mod, mod_data) in \
            enumerate(struct.conect_data.items()):

            # For each connectivity record for the current model
            # (= the atom the record refers to and the record
            # itself, containing the atoms bonded to the
            # first one and information about the bond)
            for atom1, atoms_bonded in mod_data.items():

                # For each bonded atom and attributes associated
                # with the bond
                for atom2, bond_attrs in atoms_bonded.items():

                    # Initialize an empty dictionary to store the
                    # values for each data item
                    items = {item : None for item in fmt_strings}

                    # Set the value of the 'label' and 'auth'
                    # items storing the chain the atom is on
                    items["ptnr1_label_asym_id"] = atom1[0]
                    items["ptnr1_auth_asym_id"] = atom1[0]

                    # Set the value of the 'label' and 'auth'
                    # items storing the sequence number of the
                    # residue the atom is on
                    items["ptnr1_label_seq_id"] = atom1[2][0]
                    items["ptnr1_auth_seq_id"] = atom1[2][0]

                    # Set the value of the item storing the insertion
                    # code of the residue the atom is on
                    items["pdbx_ptnr1_PDB_ins_code"] = atom1[2][1]

                    # Set the value of the 'label' and 'auth'
                    # items storing the name of the residue
                    # residue the atom is on
                    items["ptnr1_label_comp_id"] = atom1[2][2]
                    items["ptnr1_auth_comp_id"] = atom1[2][2]

                    # Set the value of the item storing the atom's
                    # name
                    items["ptnr1_label_atom_id"] = atom1[3]

                    # Set the value of the item storing the atom's
                    # alternative location
                    items["pdbx_ptnr1_label_alt_id"] = ""
                    
                    # Set the value of the 'label' and 'auth'
                    # items storing the chain the atom is on
                    items["ptnr2_label_asym_id"] = atom2[0]
                    items["ptnr2_auth_asym_id"] = atom2[0]

                    # Set the value of the 'label' and 'auth'
                    # items storing the sequence number of the
                    # residue the atom is on
                    items["ptnr2_label_seq_id"] = atom2[2][0]
                    items["ptnr2_auth_seq_id"] = atom2[2][0]

                    # Set the value of the item storing the insertion
                    # code of the residue the atom is on
                    items["pdbx_ptnr2_PDB_ins_code"] = atom2[2][1]

                    # Set the value of the 'label' and 'auth'
                    # items storing the name of the residue
                    # residue the atom is on
                    items["ptnr2_label_comp_id"] = atom2[2][2]
                    items["ptnr2_auth_comp_id"] = atom2[2][2]

                    # Set the value of the item storing the atom's
                    # name
                    items["ptnr2_label_atom_id"] = atom2[3]

                    # Set the value of the item storing the atom's
                    # alternative location
                    items["pdbx_ptnr2_label_alt_id"] = ""

                    #-------------------------------------------------#

                    # For each bond attribute and associated value
                    for bond_attr, val in bond_attrs.items():

                        # Save the value
                        items[bond_attr] = val

                    #-------------------------------------------------#

                    # Update the list storing the connectivity
                    # records to be sorted
                    sorted_conect.append(\
                        [items["conn_type_id"],
                         int(items["id"][6:]),
                         (atom1, atom2),
                         items])

            #---------------------------------------------------------#

            # Sort the list with the connectivity records by the
            # type of bond and the bond's index
            sorted_conect = \
                [elem[2:] for elem in \
                    sorted(sorted_conect, 
                            key = lambda x: (x[0], x[1]))]
            
            #---------------------------------------------------------#

            # For each pair of bonded atom and the bond's details
            for (atom1, atom2), items in sorted_conect:

                # If the bond has not been written yet
                if (atom1, atom2) not in written:

                    # Write it
                    self._write_items(\
                        file_handle = file_handle,
                        items = items,
                        fmt_strings = fmt_strings)

                    # Update the set of bonds that have
                    # already been written
                    written.add((atom1, atom2))
                    written.add((atom2, atom1))

            #---------------------------------------------------------#

            # If there are multiple models in the structure
            if i > 0:
                
                # Warn the user
                warnstr = \
                    "The structure contains several models. " \
                    "However, only connectivity data from " \
                    "the first model were written, since " \
                    "there is no way to distinguish between " \
                    "models in the '_struct_conn' data " \
                    "category. Please ensure that the " \
                    "connectivity data are valid for all models."
                logger.warning(warning)

                # Exit the loop
                break


    #------------------------ Public methods -------------------------#


    def write(self,
              struct,
              file,
              write_conect_data = True):
        """Write a ``Structure`` to a mmCIF file.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written.

        file : ``str``
            The file where the structure will be written.

        write_conect_data : ``bool``, ``True``
            Write the connectivity data associated with the
            structure, if any are present.
        """
        
        # Open the output mmCIF file
        with open(file, "w") as out:

            # Get the format strings to be used
            fmt_strings = \
                self._get_format_strings(struct = struct)

            #-----------------------------------------------------#

            # Write out the 'data_' directive and the structure's
            # name
            out.write(\
                f"{_defaults.MMCIF_DATA_DIRECTIVE}{struct.name}\n")

            # Write out a comment line
            out.write("#\n")

            #-----------------------------------------------------#

            # Write out the '_atom_type' category
            self._write_atom_type(\
                struct = struct,
                file_handle = out)

            #-----------------------------------------------------#

            # Write out the '_atom_site' category
            self._write_atom_site(\
                struct = struct,
                file_handle = out,
                fmt_strings = fmt_strings["_atom_site"])

            #-----------------------------------------------------#

            # If we need to write out connectivity data
            if write_conect_data:

                # Write those out
                self._write_struct_conn(\
                    struct = struct,
                    file_handle = out,
                    fmt_strings = fmt_strings["_struct_conn"])
