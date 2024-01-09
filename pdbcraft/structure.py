#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    structure.py
#
#    Utilities to manipulate atomic structures.
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


# Standard libaray
import collections
import copy
import logging as log
import os
import re


# Get the module's logger
logger = log.getLogger(__name__)


class Structure:

    """
    Class implementing a structure parsed from a PDB file.
    """

    # At what depth in the 'atom_data' dictionary different
    # items in the structure's hierarchy are
    _ITEM2DEPTH = \
        {# Models' level
         "model" : 1,
         # Chains' level
         "chain" : 2,
         # Segments' level
         "segment" : 3, 
         # Residues' level
         "residue" : 4,
         # Residues' attributes' level (for instance, 'res_name')
         "residue_attribute" : 5,
         # Atoms' level
         "atom" : 6,
         # Atoms' attributes' level (for instance, 'atom_name')
         "atom_attribute" : 7}

    # A mapping between residues' three-letter names and residues'
    # one-letter names (canonical residues + histidine
    # protonation states)
    _RESNAMES_3TO1 = \
        {"ALA" : "A", "CYS" : "C", "ASP" : "D", "GLU" : "E",
         "PHE" : "F", "GLY" : "G", "HIS" : "H", "ILE" : "I",
         "LYS" : "K", "LEU" : "L", "MET" : "M", "ASN" : "N",
         "PRO" : "P", "GLN" : "Q", "ARG" : "R", "SER" : "S",
         "THR" : "T", "VAL" : "V", "TRP" : "W", "TYR" : "Y",
         "HIE" : "H", "HID" : "H", "HIP" : "H"}

    # Which residues are considered protein residues
    _RESNAMES_PROTEIN = \
        ["ALA", "CYS", "ASP", "GLU",
         "PHE", "GLY", "HIS", "ILE",
         "LYS", "LEU", "MET", "ASN",
         "PRO", "GLN", "ARG", "SER",
         "THR", "VAL", "TRP", "TYR",
         "HIE", "HID", "HIP"]

    # Which residues are considered DNA residues
    _RESNAMES_DNA = ["DA", "DC", "DG", "DT"]

    # Which residues are considered RNA residues
    _RESNAMES_RNA = ["A", "C", "G", "U"]


    #------------------------ Initialization -------------------------#


    def __init__(self,
                 atom_data = None,
                 conect_data = None,
                 name = None):
        """Initialize a ``Structure`` instance.

        Parameters
        ----------
        atom_data : ``dict``, optional
            A nested dictionary containing the structure

        conect_data : ``dict``, optional
            A nested dictionary containing the CONECT data
            for the structure.

        name : ``str``, optional
            The name of the structure.

        Notes
        -----
        To create an empty ``Structure``, leave both ``atom_data``
        and ``conect_data`` blank.
        """

        #---------------- Set the atomic coordinates -----------------#


        # If no atomic coordinates were provided
        if atom_data is None:

            # Set the atomic coordinates to an empty dictionary
            self._atom_data = {}

        # Otherwise
        else:
            
            # Set the atomic coordinates to the dictionary
            # provided
            self._atom_data = atom_data


        #-------------------- Set the CONECT data --------------------#


        # If no CONECT data were provided
        if conect_data is None:

            # If no atoms' coordinates were provided
            if atom_data is None:

                # Set the CONECT data to an empty dictionary
                self._conect_data = {}

            # Otherwise
            else:

                # Set the CONECT data to a dictionary containing
                # an empty dictionary for each of the models
                # found in the atomic coordinates
                self._conect_data = \
                    {mod : {} for mod in self.atom_data}

        # Otherwise
        else:
            
            # Set the CONECT data to the dictionary provided
            self._conect_data = conect_data


        #----------------------- Set the name ------------------------#


        # Set the name of the structure
        self._name = name


    #-------------------------- Properties ---------------------------#


    @property
    def atom_data(self):
        """A dictionary containing the atomic coordinates.
        """
        
        return self._atom_data


    @atom_data.setter
    def atom_data(self,
                  value):
        """Raise an error if the user tries to modify the
        'atom_data' attribute.
        """

        errstr = \
            "The 'atom_data' attribute cannot be directly modified."
        raise ValueError(errstr)


    @property
    def conect_data(self):
        """A dictionary containing the CONECT data associated
        with the structure.
        """
        
        return self._conect_data


    @conect_data.setter
    def conect_data(self,
                    value):
        """Raise an error if the user tries to modify the
        'conect_data' attribute.
        """

        errstr = \
            "The 'conect_data' attribute cannot be directly modified."
        raise ValueError(errstr)


    @property
    def name(self):
        """The name of the structure.
        """
        
        return self._name


    #----------------------- "Dunder" methods ------------------------#


    def __eq__(self,
               other):
        """Return ``True`` if all attributes of ``self`` are equal
        to the attributes or ``other`` (the structures' names can be
        different).

        Parameters
        ----------
        other: ``pdbcraft.structure.Structure``
            Another structure.
        """

        # If 'other' is not a 'Structure'
        if not isinstance(other, self.__class__):

            # Raise an error
            errstr = \
                "Equality and inequality operations are " \
                "supported only between " \
                f"'{self.__class__.__name__}' instances."
            raise NotImplementedError(errstr)

        # Check the equality
        return (self.atom_data == other.atom_data \
                and self.conect_data == other.conect_data)

    def __ne__(self,
               other):
        """Return ``True`` if any of the data attributes of ``self``
        is different from the corresponding attribute of ``other``
        (the structures' names can be different).

        Parameters
        ----------
        other: ``pdbcraft.structure.Structure``
            Another structure.
        """

        # Check the inequality
        return not self.__eq__(other = other)

    

    #------------------------ Private methods ------------------------#


    def _update_conect_data(self,
                            mapping):
        """Update the CONECT data for the structure.

        Parameters
        ----------
        mapping : ``dict``
            A dictionary mapping the old atom numbering to
            the new atom numbering.
        """

        # If the structure has any associated CONECT data
        if self.conect_data:

            # Initialize a dictionary to store the updated CONECT
            # data, creating an empty dictionary for each model
            # in the structure
            conect_data_updated = \
                {mod : {} for mod in self.conect_data}

            # For each model in the CONECT data
            for mod in self.conect_data:

                # For each CONECT record for the current model
                # (= the atom the record refers to and the record
                # itself, containing the atoms bonded to the
                # first one)
                for atom_record, record \
                    in self.conect_data[mod].items():

                    # Try to convert the atom's serial number
                    # to the new numbering
                    try:

                        atom_record_updated = \
                            mapping[mod][atom_record]

                    # If the atom was not found in the mapping
                    # (because, for instance, it was removed
                    # from the structure)
                    except KeyError:

                        # Go to the next atom
                        continue

                    # Initialize an empty list to store the new
                    # record
                    record_updated = []

                    # For each atom in the record
                    for atom_bonded in record:

                        # Try to convert the atom's serial number
                        # to the new numbering
                        try:

                            record_updated.append(\
                                mapping[mod][atom_bonded])

                        # If the atom was not found (because, for
                        # instance, it was removed from the structure)
                        except KeyError:

                            # Go to the next atom
                            continue

                    # If the length of the new record is equal to the
                    # length of the old record (= all atoms were added,
                    # meaning that no atoms were missing from the
                    # current record)
                    if len(record_updated) == len(record):

                        # Add the record to the dictionary of
                        # updated CONECT data
                        conect_data_updated[mod][\
                            atom_record_updated] = record_updated

            # Update the CONECT data for the structure
            self._conect_data = conect_data_updated


    def _add_bonds(self,
                   bonds):
        """Add bonds between atoms (it modifies the CONECT data,
        while the atomic coordinates remain unchanged).

        Parameters
        ----------
        bonds : ``list``
            A list of tuples, each representing a bond between two
            atoms and whether the bond is a ``"single"``, 
            ``"double"``, or ``"triple"`` bond.

            An atom can be represented either by its unique
            serial number or by its "path", namely a tuple
            containing the model, chain, segment, and residue
            the atom belongs to, plus the atom's name. This is
            based on the assumption that, within each residue,
            atoms must have unique names.

            For instance, a single bond between atom 1 and atom 2
            could be represented as ``(1, 2, "single")`` or as
            ``((1, "A", "", (1, " "), "N"), 
            (1, "A", "", (1, " "), "CA"),
            "single")``, assuming that both atoms belong to
            model 1, chain A, an unnamed segment, residue 1 with
            insertion code " " (= no insertion code), and one is
            named N, while the other is named CA. The model's
            specification is necessary because atoms in the same
            chain, segment, and residue may have different serial
            numbers in different models.

            Each bond should be specified only once (for instance,
            we should not specify a bond between atom 2 and atom 1
            if we have already specified a bond between atom 1 and
            atom 2) since, for each bond specified, the CONECT records
            for both atoms involved in the bond will be updated.
        """

        # Create a copy of the current CONECT data that will
        # be updated with the new bonds
        conect_data_updated = copy.deepcopy(self.conect_data)

        # For each bond (first atom, second atom, and bond type)
        for atom_1, atom_2, bond_type in bonds:

            # If the first atom is defined as a path to an atom
            if isinstance(atom_1, tuple):

                # Get the atom's serial number from its path
                atom_1 = self._get_atom_serial(path = atom_1)

            # If the second atom is defined as a path to an atom
            if isinstance(atom_2, tuple):

                # Get the atom's serial number from its path
                atom_2 = self._get_atom_serial(path = atom_2)

            # If the bond is a single bond
            if bond_type == "single":

                # Create the CONECT records representing
                # the new bond for both atoms
                new_record_1, new_record_2 = \
                    [atom_2], [atom_1]

            # If the bond is a double bond
            elif bond_type == "double":

                # Create the CONECT records representing
                # the new bond for both atoms
                new_record_1, new_record_2 = \
                    [atom_2, atom_2], [atom_1, atom_1]

            # If the bond is a triple bond
            elif bond_type == "triple":
                
                # Create the CONECT records representing
                # the new bond for both atoms
                new_record_1, new_record_2 = \
                    [atom_2, atom_2, atom_2], [atom_1, atom_1, atom_1]

            # For each model's CONECT records
            for mod in self.conect_data:

                # If we already have a CONECT record for the
                # first atom
                if atom_1 in conect_data_updated[mod]:

                    # Update it with the new record
                    record_1 = \
                        [*conect_data_updated[mod][atom_1],
                         *new_record_1]

                # Otherwise
                else:

                    # Create a new record for the first atom
                    record_1 = new_record_1

                # Update/add the record to the CONECT data,
                # making sure to sort the atoms in the record
                # in ascending order of serial number
                conect_data_updated[mod][atom_1] = sorted(record_1)

                # If we already have a CONECT record for the
                # second atom
                if atom_2 in conect_data_updated[mod]:

                    # Update it with the new record
                    record_2 = \
                        [*conect_data_updated[mod][atom_2],
                         *new_record_2]

                # Otherwise
                else:

                    # Create a new record for the second atom
                    record_2 = new_record_2
                
                # Update/add the record to the CONECT data,
                # making sure to sort the atoms in the record
                # in ascending order of serial number
                conect_data_updated[mod][atom_2] = sorted(record_2)

        # Sort the atoms in each model so that the CONECT data
        # are reported in the correct order (= in ascending
        # order of the atoms' serial numbers)
        conect_data_sorted = \
            {mod : \
                {k : conect_data_updated[mod][k] for k  \
                 in sorted(conect_data_updated[mod])} \
             for mod in conect_data_updated}

        # Update the CONECT data for the structure
        self._conect_data = conect_data_sorted


    def _get_atom_serial(self,
                         path):
        """Get the serial number of a specific atom, given
        the 'path' to the atom in the structure in terms
        of the model, chain, segment, and residue the atom
        belongs to, plus the atom's name.

        This method assumes that, within each residue, the
        atoms' names are unique.

        Parameters
        ----------
        path : ``tuple``
            The 'path' to the atom (the model, chain, segment,
            and residue it belongs to, plus the atom's name).

        Returns
        -------
        atom_serial : ``int```
            The serial number of the atom the 'path' points to.
        """
        
        # Get the model, segment, chain, residue, and name
        # of the atom whose serial number should be returned
        mod, ch, seg, res, atom_name = path

        # Get the sub-structure where the atom is stored
        sub_struct = self.atom_data[mod][ch][seg][res]["atoms"]

        # For each atom in the sub-structure
        for atom in sub_struct:

            # If the current atom's name matches the one
            # we are looking for
            if sub_struct[atom]["atom_name"] == atom_name:
                
                # Return the atom's serial number (there is no
                # need to check any further because we are
                # assuming that, for each residue, the atoms'
                # names are unique)
                return atom


    def _get_last(self,
                  items_type):
        """Get the last item of a specfic type in the structure
        by recursively traversing the 'atom_data' dictionary.

        Parameters
        ----------
        items_type : ``str``
            The type of items to get the last of (model, chain,
            segment, residue, or atom).

        Returns
        -------
        last_item : ``str`` or ``int`` or ``tuple``
            The last item of the specified type (model, chain,
            segment, residue, or atom).
        """


        #------------------- Define the recursion --------------------#


        def recursive_step(struct,
                           target_depth,
                           current_depth = 1):
            
            # If the current sub-structure is not a dictionary
            if not isinstance(struct, dict):

                # Return None
                return None

            # If the current depth is the target depth
            if current_depth == target_depth:

                # Return the last (= highest) key found in
                # the sub-structure
                return max(struct.keys(),
                           default = None)

            # Set the default last item to None
            last_item = None

            # For each sub-structure in the current dictionary
            for sub_struct in struct.values():

                # Get the last (= maximum) item in the sub-structure
                last_in_sub_struct = \
                    recursive_step(struct = sub_struct,
                                   target_depth = target_depth,
                                   current_depth = current_depth + 1)
                
                # If the recursion returned an item that is not None
                if last_in_sub_struct is not None:

                    # If the current last item is None or the
                    # item returned is higher than the current
                    # last item
                    if last_item is None \
                    or last_in_sub_struct > last_item:
                        
                        # Update the last item with the current one
                        last_item = last_in_sub_struct

            # Return the last item
            return last_item


        #--------------------- Get the last item ---------------------#


        # Get the depth corresponding to where the items we want
        # to take the last of are stored in the structure
        target_depth = self._ITEM2DEPTH[items_type]

        # Find the last item of the specified type (we can pass the
        # 'atom_data' attribute directly because it will not be
        # modified by the recursion)
        return recursive_step(struct = self.atom_data,
                              target_depth = target_depth)


    def _merge(self,
               other):
        """Merge two ``Structure`` instances and return
        a new ``Structure``.

        Parameters
        ----------
        other : ``pdbcraft.structure.Structure``
            The ``Structure`` to be merged with the current one.

        Returns
        -------
        merged : ``pdbcraft.structure.Structure``
            The merged ``Structure``.
        """


        #------------------- Define the recursion --------------------#


        def recursive_step(struct,
                           other_struct):

            # If none of the current sub-structures is
            # a dictionary
            if not isinstance(struct, dict) \
            or not isinstance(other_struct, dict):

                # Return the first one as it is
                return struct

            # Create a shallow copy of the first sub-structure,
            # representing the initial 'merged' sub-structure
            merged = dict(struct)

            # For each key and associated sub-structure in
            # the 'other' current sub-structure
            for key, sub_struct in other_struct.items():

                # If the key is already in the merged structure
                if key in merged:

                    # Add the values found associated with the same
                    # key in the 'other' sub-structure
                    merged[key] = \
                        recursive_step(struct = merged[key],
                                       other_struct = sub_struct)
                
                # Otherwise (= the key is only in the 'other'
                # sub-structure)
                else:

                    # Add the entire sub-structure the key is
                    # associated with to the merged structure
                    merged[key] = sub_struct

            # Return the merged structure
            return merged


        #------------------- Merge the structures --------------------#


        # If 'other' is not a Structure
        if not isinstance(other, self.__class__):

            # Raise an error
            errstr = \
                "'other' must be a Structure instance."
            raise TypeError(errstr)

        # If the two structures do not have the same number
        # of models
        if len(self.atom_data) != len(other.atom_data):

            # Raise an error
            errstr = \
                "The current structure and 'other' do not " \
                "have the same number of models. Therefore, " \
                "they cannot be merged."
            raise ValueError(errstr)

        # Get the serial number of the last atom in the current
        # structure
        last_atom = self._get_last(items_type = "atom")

        # Renumber the atoms in the 'other' structure starting
        # from the last atom in the current structure
        other._renumber_atoms(start = last_atom)

        # Merge the atomic coordinates of the two structures
        merged_atom_data = \
            recursive_step(struct = self.atom_data,
                           other_struct = other.atom_data)

        # Merge the CONECT data of the two structures
        merged_conect_data = \
            recursive_step(struct = self.conect_data,
                           other_struct = other.conect_data)

        # Create the merged structure
        merged = \
            self.__class__(atom_data = merged_atom_data,
                           conect_data = merged_conect_data)

        # Renumber the atoms of the merged structure starting
        # from 1
        merged._renumber_atoms()

        # Return the merged structure
        return merged


    def _relabel(self,
                 action,
                 item,
                 mapping = None,
                 start = None,
                 attribute = None,
                 models = None,
                 chains = None,
                 segments = None,
                 residues = None,
                 residues_types = None):
        """Re-label models, chains, segments, residues, or atoms,
        or assigm a new value to one of their attributes.

        Parameters
        ----------
        action : ``str``, {``"renumber"``, ``"rename"``}
            The re-labeling action performed (either renumbering
            or renaming).

        item : ``str``
            The type of item to be re-labeled.

        mapping : ``dict``, optional
            A dictionary mapping the old labels to the new labels.

        start : ``int``, optional
            The new starting number for items whose IDs are either
            integers (models) or contain integers (residues).

        attribute : ``str``, optional
            The name of the attribute that should be relabeled.

        models : ``list``, optional
            A list of the models whose items will be affected by
            the relabeling.

        chains : ``list``, optional
            A list of chains whose items will be affected by
            the relabeling.

        segments : ``list``, optional
            A list of segments whose items will be affected by
            the relabeling.
        
        residues : ``list``, optional
            A list of residues whose items will be affected by
            the relabeling.

        residues_types : ``list``, optional
            A list of residues' types. Items belonging to residues
            of the selected types will be affected by the relabeling.
        """


        #------------------- Define the recursion --------------------#


        def recursive_step(struct,
                           target_depth,
                           attribute,
                           mapping, 
                           start,
                           models,
                           chains,
                           segments,
                           residues,
                           residues_types,
                           current_depth = 1):

            # If we are at the target depth
            if current_depth == target_depth:

                # If we are assigning a new value to a specific
                # attribute
                if attribute is not None:

                    # Get the current value of the attribute
                    value = struct[attribute]

                    # If the value is included in the mapping
                    # provided
                    if value in mapping:

                        # Substitute the value with the corresponding
                        # one found in the mapping
                        struct[attribute] = mapping[value]

                    # Return the current sub-structure
                    return struct

                # If we are substituting the identifiers (and not
                # specific attributes)
                else:
            
                    # If a starting number was passed
                    if start is not None:
                    
                        # Renumber the keys in the current
                        # sub-structure, taking care of the
                        # identifiers that are tuples
                        # (= residue idenfitiers)
                        return {((i, old_key[1]) \
                                  if isinstance(old_key, tuple) \
                                  else i) : struct[old_key]  \
                                for i, old_key \
                                in enumerate(struct, start)}

                    # If a mapping was passed
                    elif mapping is not None:

                        # Create an empty dictionary to store
                        # the updated sub-structure
                        new_struct = {}

                        # For each key in the current sub-structure
                        for old_key in struct:

                            # By default, the key in the updated
                            # sub-structure will be the same as
                            # in the old sub-structure
                            new_key = old_key

                            # If the old key is a tuple (= residue
                            # identifier)
                            if isinstance(old_key, tuple):

                                # If the first element of the tuple
                                # (= the residue's sequence number)
                                # is in the mapping
                                if old_key[0] in mapping:

                                    # Create the key for the updated
                                    # sub-structure
                                    new_key = \
                                        (mapping[old_key[0]],
                                         old_key[1])

                            # Otherwise
                            else:

                                # If the entire old key is in the
                                # mapping
                                if old_key in mapping:

                                    # Create the key for the updated
                                    # sub-structure
                                    new_key = mapping[old_key]

                            # Add the key and the associated value
                            # to the updated struture
                            new_struct[new_key] = struct[old_key]

                        # Sort the updated structure
                        new_struct = \
                            {key : new_struct[key] for key in \
                             sorted(new_struct)}

                        # Return the updated structure
                        return new_struct  

            # If we have not reached the target depth yet, for
            # each key and associated sub-structure in the current
            # sub-structure
            for key, sub_struct in struct.items():

                # Set the flag to decide whether we should
                # recurse to True (= by default, recurse)
                should_recurse = True
                
                # If:
                # - a list of models was specified
                # - AND we have reached the models' level
                # - AND we are not at any of the selected models
                if models is not None and current_depth == 1 \
                and key not in models:

                    # Do not recurse
                    should_recurse = False
                
                # If:
                # - a list of chains was specified
                # - AND we have reached the chains' level
                # - AND we are not at any of the selected chains
                if chains is not None and current_depth == 2 \
                and key not in chains:

                    # Do not recurse
                    should_recurse = False
                
                # If:
                # - a set of segments was specified
                # - AND we have reached the segments' level
                # - AND we are not at any of the selected segments
                if segments is not None and current_depth == 3 \
                and key not in segments:

                    # Do not recurse
                    should_recurse = False
                
                # If:
                # - a list of residue was specified
                # - AND and we have reached the residues' level
                if residues is not None and current_depth == 4:

                    # If neither the full residue ID (sequence
                    # number and insertion code) or the residue's
                    # sequence number are in the set of selected
                    # residues
                    if key not in residues and key[0] not in residues:

                        # Do not recurse
                        should_recurse = False
                
                # If:
                # - a list of residue types was specified
                # - AND we have reached the inner residues'
                #   attributes' level
                # - AND and the current residue is not of a type
                #   of interest
                if residues_types is not None and current_depth == 5 \
                and struct["res_name"] in residues_types:

                    # Do not recurse
                    should_recurse = True

                # If we should recurse and the current sub-structure
                # is a dictionary
                if should_recurse and isinstance(sub_struct, dict):

                    # Recurse over the sub-structure, substituting
                    # the sub-structure associated to the current key
                    # with the result of the recursion
                    struct[key] = \
                        recursive_step(\
                            struct = sub_struct,
                            target_depth = target_depth,
                            attribute = attribute,
                            mapping = mapping,
                            start = start,
                            models = models,
                            chains = chains,
                            segments = segments,
                            residues = residues,
                            residues_types = residues_types,
                            current_depth = current_depth + 1)

            # Return the updated structure
            return struct


        #--------------------------- Check ---------------------------#


        # If both 'start' and 'mapping' were specified
        if start is not None and mapping is not None:

            # Raise an error
            errstr = \
                "'start' and 'mapping' are mutually exclusive " \
                "options."
            raise ValueError(errstr)

        # If a starting number was specified
        if start is not None:

            # If the starting number is not an integer
            if not isinstance(start, int):

                # Raise an error
                errstr = "'start' must be an integer."
                raise TypeError(errstr)

            # If a starting number was specified for something
            # other than models or residues
            if item not in ("model", "residue"):

                # Raise an error
                errstr = \
                    "'start' can only be used for re-numbering " \
                    "models and residues."
                raise ValueError(errstr)

        # If a mapping was specified
        if mapping is not None:

            # If the mapping is not a dictionary
            if not isinstance(mapping, dict):

                # Raise an error
                errstr = "'mapping' must be a dictionary."
                raise TypeError(errstr)

            # If we are renumbering with the mapping
            if action == "renumber":

                # Get the values in the mapping as a list and
                # as a set
                values = list(mapping.values())
                values_set = set(mapping.values())

                # If they are not unique
                if len(values_set) < len(values):

                    # Raise an error
                    errstr = \
                        "When renumbering, both keys and values " \
                        "in 'mapping' must be unique."
                    raise ValueError(errstr)


        #------------------------- Re-label --------------------------#


        # Get the target depth from the type of items to be re-labeled
        target_depth = self._ITEM2DEPTH[item]

        # Recurse and update the structure's atomic coordinates
        # (we need to pass a copy of the 'atom_data' dictionary
        # beecause it will be modified by the recursion)
        self._atom_data = \
            recursive_step(struct = copy.deepcopy(self.atom_data),
                           target_depth = target_depth,
                           start = start,
                           mapping = mapping,
                           attribute = attribute,
                           models = models,
                           chains = chains,
                           segments = segments,
                           residues = residues,
                           residues_types = residues_types)


    def _renumber_atoms(self,
                        start = 1):
        """Renumber the atoms in the structure so that each model's
        atom numbering starts from 'start', and each chain's atom
        numbering starts from the number of the atom ending the
        previous chains plus 2 (to accommodate the TER record at
        the end of a chain).

        Parameters
        ----------
        start : ``int``, ``1``
            The new starting number for the atoms' serial numbers.
        """


        #------------------- Define the recursion --------------------#


        def recursive_step(struct,
                           current_atom_count,
                           current_depth = 1,
                           current_path = (),
                           oldnum2newnum = None):

            # If no mapping has been created yet between the
            # old numbering and the new numbering, create an
            # empty dictionary that will contain the mapping (this
            # mapping will be created only once when we are
            # traversing the dictionary representing the first
            # model, and then updated, since the mapping will
            # be None only at the very beginning)
            if oldnum2newnum is None:

                # Create a mapping that will contain, for each
                # model, the old atom numbering mapped to
                # the new atom numbering
                oldnum2newnum = {mod : {} for mod in struct}
            
            # If we have reached the atoms' level
            if current_depth == 6:

                # Set the default starting number for the current
                # set of atoms
                start_num = current_atom_count + 1

                # Renumber the atoms
                atoms_renumbered = \
                    {atom_rel_index + start_num : atom_data \
                     for atom_rel_index, atom_data \
                     in enumerate(struct.values())}
                
                # Update the atoms' count for the current model
                # with the atoms of the current set
                current_atom_count = \
                    start_num + len(struct) - 1

                # Update the mapping between the old numbering
                # and the new numbering with the mapping for
                # the current set of atoms
                oldnum2newnum[current_path[0]].update(\
                    {atom_old_num : atom_rel_index + start_num \
                     for atom_rel_index, (atom_old_num, _) \
                     in enumerate(struct.items())})
                
                # Return the renumbered atoms, the current
                # atom count, and the updated mapping
                return atoms_renumbered, \
                       current_atom_count, \
                       oldnum2newnum

            # If we have not reached the atoms yet
            else:
                
                # Get the previous model
                prev_model = \
                    current_path[0] if len(current_path) > 0 else None
                
                # Get the previous chain
                prev_chain = \
                    current_path[1] if len(current_path) > 1 else None
                
                # For each key and associated sub-structure in the
                # current sub-structure
                for key, sub_struct in struct.items():
                    
                    # If the sub-structure is a dictionary
                    if isinstance(sub_struct, dict):
                        
                        # If we are in a new model
                        if current_depth == 1 and key != prev_model:
                            
                            # Reset the count for the atoms
                            current_atom_count = -1
                        
                        # If we are in a new chain
                        elif current_depth == 2 and key != prev_chain:
                            
                            # Carry over the previous count plus 1
                            current_atom_count += 1

                        # Recursively process the sub-structure
                        struct[key], current_atom_count, oldnum2newnum = \
                            recursive_step(\
                                struct = sub_struct,
                                current_atom_count = current_atom_count,
                                current_depth = current_depth + 1,
                                current_path = current_path + (key,),
                                oldnum2newnum = oldnum2newnum)
                
                # Return the structure, the current atom count,
                # and the mapping between the old and the new
                # numbering
                return struct, current_atom_count, oldnum2newnum


        #------------------------- Renumber --------------------------#


        # Get the updated data for the atoms and the mapping
        # between the old numbering and the new numbering (we need
        # to pass a copy of the 'atom_data' dictionary because it
        # will be modified by the recursion)
        atom_data_updated, _, oldnum2newnum = \
            recursive_step(struct = copy.deepcopy(self.atom_data),
                           current_atom_count = start - 2)

        # Update the data for the atoms
        self._atom_data = atom_data_updated


        #-------------------- Update CONECT data ---------------------#


        # Update the CONECT data associated with the structure
        self._update_conect_data(mapping = oldnum2newnum)


    def _keep_or_remove(self,
                        action,
                        items,
                        items_type,
                        models = None,
                        chains = None,
                        segments = None,
                        residues = None):
        """Keep or remove models, chains, segments, residues, or atoms.

        Parameters
        ----------
        action : ``str``, {``"keep"``, ``"remove"``}
            Whether to keep or remove the selected items.

        items : ``list``
            The list of items to be kept/removed.

        items_type : ``str``
            The type of items to be kept/removed.

        models : ``list``, optional
            The models to be considered when keeping/removing
            the items.

            If not provided, all models will be considered.

        chains : ``list``, optional
            The chains to be considered when keeping/removing
            the items.

            If not provided, all chains will be considered.

        segments : ``list``, optional
            The segments to be considered when keeping/removing
            the items.

            If not provided, all segments will be considered.

        residues : ``list``, optional
            The residues to be considered when keeping/removing
            the items.

            If not provided, all residues will be considered.
        """
        

        #------------------- Define the recursion --------------------#


        def recursive_step(struct,
                           action,
                           items,
                           target_depth,
                           models,
                           chains,
                           segments,
                           residues,
                           current_depth = 1):
            
            # If we are at the target depth
            if current_depth == target_depth:

                # If we have to keep only the selected items
                if action == "keep":

                    # Return the current sub-structure with the
                    # items to be kept (deal also with the case of
                    # residues, where you have a tuple representing
                    # the item but the user may have passed only the
                    # sequence number or a residue name)
                    return {k : v for k, v in struct.items() \
                            if (isinstance(k, tuple) \
                                and (k[0] in items \
                                or v["res_name"] in items)) \
                            or (not isinstance(k, tuple) \
                                and k in items)}

                # If we have to remove the selected items
                elif action == "remove":

                    # Return the current sub-structure without the
                    # items to be removed (deal also with the case of
                    # residues, where you have a tuple representing
                    # the item but the user may have passed only the
                    # sequence number or a residue name)
                    return {k : v for k, v in struct.items() \
                            if (isinstance(k, tuple) \
                                and (k[0] not in items \
                                and v["res_name"] not in items)) \
                            or (not isinstance(k, tuple) \
                                and k not in items)}

            # For each key and associated sub-structure in the
            # current sub-structure
            for key, sub_struct in list(struct.items()):
                
                # If the sub-structure is a dictionary
                if isinstance(sub_struct, dict):

                    # Evaluate the condition for whether this is
                    # the right model (we are at the right depth
                    # and either the user did not select any models
                    # or the current model is among the selected
                    # ones)
                    is_right_model = \
                        current_depth == 1 and \
                        (models is None or key in models)

                    # Evaluate the condition for whether this is
                    # the right chain (we are at the right depth
                    # and either the user did not select any chains
                    # or the current chain is among the selected
                    # ones)
                    is_right_chain = \
                        current_depth == 2 and \
                        (chains is None or key in chains)

                    # Evaluate the condition for whether this is
                    # the right segment (we are at the right depth
                    # and either the user did not select any segments
                    # or the current segment is among the selected
                    # ones)
                    is_right_segment = \
                        current_depth == 3 and \
                        (segments is None or key in segments)

                    # Evaluate the condition for whether this is
                    # the right residue (we are at the right depth
                    # and either the user did not select any residues
                    # or the current residue is among the selected
                    # ones)
                    is_right_residue = \
                        current_depth == 4 and \
                        (residues is None or key in residues)

                    # If the conditions for model, chain, segment,
                    # or residue are met
                    if is_right_model or is_right_chain \
                        or is_right_segment or is_right_residue:
                        
                        # Recurse through the sub-structure
                        struct[key] = \
                            recursive_step(\
                                struct = sub_struct,
                                action = action,
                                items = items,
                                target_depth = target_depth,
                                models = models,
                                chains = chains,
                                segments = segments,
                                residues = residues,
                                current_depth = current_depth + 1)

            # Return the updated structure
            return struct


        #------------------- Keep/remove the items -------------------#


        # Get the target depth corresponding to the items to
        # be kept or removed
        target_depth = self._ITEM2DEPTH[items_type]

        # Keep or remove the items (we need to pass a copy of the
        # 'atom_data' attribute because it will be modified by the
        # recursion)
        self._atom_data = \
            recursive_step(struct = copy.deepcopy(self.atom_data),
                           action = action,
                           items = items,
                           target_depth = target_depth,
                           models = models,
                           chains = chains,
                           segments = segments,
                           residues = residues)


        #-------------------- Renumber the atoms ---------------------#


        # Renumber the atoms (it also updates the CONECT data)
        self._renumber_atoms()


    def _move_atoms(self,
                    atoms,
                    create_new_struct,
                    new_chain,
                    new_segment,
                    new_residue,
                    in_place):
        """Move atoms to another part of the structure or to a new
        structure.

        Parameters
        ----------
        atoms : ``list``
            The atoms to be moved.

        create_new_struct : ``bool``
            Whether a new structure should be created to store the
            atoms.

        new_chain : ``str``
            The identifier of the chain where the atoms should
            be moved.

        new_segment : ``str``
            The identifier of the segment where the atoms should
            be moved.

        new_residue : ``tuple``
            The identifier of the residue where the atoms should
            be moved.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.structure.Structure`` if ``in_place = False``,
        ``None`` otherwise.
        """


        #------------------- Define the recursion --------------------#


        def recursive_step(struct_iter,
                           struct_add,
                           struct_remove,
                           atoms,
                           new_chain,
                           new_segment,
                           new_residue,
                           current_depth = 1,
                           current_path = ()):
            
            # If we have reached the atoms' level (it is at depth
            # 5 in this case because we 'skipped' a level when
            # defining the sub-structure to recurse through at
            # depth 4)
            if current_depth == 5:

                # Get the model, chain, segment, and residue
                # the atoms belong to
                model, chain, segment, residue = current_path

                # For each atom in the current residue
                for atom_serial in struct_iter:

                    # If the atom is one of those to be moved
                    if atom_serial in atoms:

                        # If no new chain was specified
                        if new_chain is None:

                            # The atom will stay on the current
                            # chain
                            new_chain = chain

                        # If no new segment was specified
                        if new_segment is None:

                            # The atom will stay on the current
                            # segment
                            new_segment = segment

                        # If no new residue was specified
                        if new_residue is None:

                            # The atom will stay on the current
                            # residue
                            new_residue_seq = residue

                            # The residue name will be the one
                            # of the current residue
                            new_residue_name = \
                                struct_remove[model][chain][\
                                    segment][residue]["res_name"]

                        # Otherwise
                        else:

                            # Separate the residue's sequence number
                            # and insertion code from the residue's
                            # name
                            new_residue_seq, new_residue_name = \
                                new_residue[:2], new_residue[2]
                        
                        # In the structure to which the atoms should
                        # be added, add a new model with the ID
                        # of current model, if necessary
                        struct_add.setdefault(\
                            model,
                            {})
                        
                        # In the structure to which the atoms should
                        # be added, add a new chain with the
                        # specified ID, if necessary
                        struct_add[model].setdefault(\
                            new_chain,
                            {})

                        # In the structure to which the atoms should
                        # be added, add a new segment with the
                        # specified ID, if necessary
                        struct_add[model][new_chain].setdefault(\
                            new_segment,
                            {})

                        # In the structure to which the atoms should
                        # be added, add a new residue with the
                        # specified ID and name, if necessary
                        struct_add[model][new_chain][\
                            new_segment].setdefault(\
                                new_residue_seq, 
                                {"res_name" : new_residue_name,
                                 "atoms": {}})
                        
                        # In the structure to which the atoms should
                        # be added, add the current atom
                        struct_add[model][new_chain][\
                            new_segment][new_residue_seq][\
                                "atoms"][atom_serial] = \
                                    struct_iter[atom_serial]

                        # Remove the atom from the structure from
                        # which the atoms should be removed
                        del struct_remove[model][chain][\
                                segment][residue]["atoms"][\
                                    atom_serial]

                # If removing the atoms resulted in an empty residue
                if not struct_remove[model][chain][\
                    segment][residue]["atoms"]:

                    # Remove the residue
                    del struct_remove[model][chain][segment][residue]

                # If removing the atoms resulted in an empty segment
                if not struct_remove[model][chain][segment]:

                    # Remove the segment
                    del struct_remove[model][chain][segment]

                # If removing the atoms resulted in an empty chain
                if not struct_remove[model][chain]:

                    # Remove the chain
                    del struct_remove[model][chain]

                # If removing the atoms resulted in an empty model
                if not struct_remove[model]:

                    # Remove the model
                    del struct_remove[model]
            
            # If we have not reached the atoms' level
            else:

                # For each key and associated sub-structure
                # in the current level of the structure
                for key, sub_dict in list(struct_iter.items()):
                    
                    # If we reached the residues' level
                    if current_depth == 4:

                        # Set the next level to be the atoms'
                        # level (since they are stored in the
                        # 'atoms' dictionary within each
                        # residue)
                        sub_dict = struct_iter[key]["atoms"]
                    
                    # If we are at the models' level
                    elif current_depth == 1:

                        # Reset the current path
                        current_path = ()

                    # Recursively traverse the sub-dictionary
                    recursive_step(\
                        struct_iter = sub_dict,
                        struct_add = struct_add,
                        struct_remove = struct_remove,
                        atoms = atoms,
                        new_chain = new_chain,
                        new_segment = new_segment,
                        new_residue = new_residue,
                        current_depth = current_depth + 1,
                        current_path = current_path + (key,))


        #---------------------- Move the atoms -----------------------#


        # If we need to add the atoms to a new structure
        if not in_place:

            # Create an empty structure
            new_struct = Structure()

            # Move the atoms from one structure to the other
            recursive_step(\
                struct_iter = copy.deepcopy(self.atom_data),
                struct_add = new_struct._atom_data,
                struct_remove = self._atom_data,
                atoms = atoms,
                new_chain = new_chain,
                new_segment = new_segment,
                new_residue = new_residue)

            # Copy the CONECT data from the old structure to the
            # new structure
            new_struct._conect_data = copy.deepcopy(self.conect_data)

            # Renumber the atoms in the current structure (it also
            # updates the CONECT data, removing the records for
            # atoms no longer in the structure)
            self._renumber_atoms()

            # Renumber the atoms in the new structure (it also
            # updated the CONECT data, removing the records for
            # atoms no longer in the structure)
            new_struct._renumber_atoms()
            
            # Return the new structure
            return new_struct

        # If we need to move the atoms within the same structure
        else:

            # Move the atoms within the structure
            recursive_step(\
                    struct_iter = copy.deepcopy(self.atom_data),
                    struct_add = self._atom_data,
                    struct_remove = self._atom_data,
                    atoms = atoms,
                    new_chain = new_chain,
                    new_segment = new_segment,
                    new_residue = new_residue)

            # Renumber the atoms in the current structure (it also
            # updates the CONECT data, removing the records for
            # atoms no longer in the structure)
            self._renumber_atoms()


    def _sort_atoms_residue_name(self,
                                 residue_name,
                                 atoms_names,
                                 other_atoms_position):
        """Sort the atoms in residues having a specific name.

        Parameters
        ----------
        residue_name : ``str``
            The name of the residues whose atoms will be sorted.

        atoms_names : ``list``
            The names of the atoms in residues named
            ``residue_name``, sorted.

        other_atoms_position : ``str``, {``"before"``, ``"after"``}
            Where to put the atoms not included in ``atoms_names``,
            if found in the residues named ``residue_name``:
            ``"before"`` or ``"after"`` the atoms specified in
            ``atoms_names``.
        """


        #------------------- Define the recursion --------------------#


        def recursive_step(struct,
                           residue_name,
                           atoms_names,
                           other_atoms_position,
                           current_depth = 1,
                           current_path = ()):

            # If we have reached the atoms' level (it is at depth
            # 5 in this case because we 'skipped' a level when
            # defining the sub-structure to recurse through at
            # depth 4)
            if current_depth == 5: 

                # Create a new dictionary to store the atoms
                # in the new order
                new_atoms_order = {}
                
                # Get the atoms matching the names of those provided
                expected_atoms = \
                    {atom_serial : atom_data \
                     for atom_serial, atom_data \
                     in struct.items() \
                     if atom_data["atom_name"] in atoms_names}
                
                # Get the extra atoms in the current residue
                # (atoms whose names are not included in the
                # list of atoms to sort)
                other_atoms = \
                    {atom_serial : atom_data \
                     for atom_serial, atom_data \
                     in struct.items() \
                     if atom_data["atom_name"] not in atoms_names}

                # Get the atoms than were found in the current residue
                found_atoms = \
                    [atom_data["atom_name"] for atom_data \
                    in expected_atoms.values()]
                
                # Get the atoms missing from the residues by comparing
                # the set of atoms' names in the list passed and
                # the names of the atoms found in the residue
                missing_atoms = \
                    set(atoms_names) - set(found_atoms)
                
                # If any atoms is missing
                if missing_atoms:

                    # Get the current model, chain, segment,
                    # and residue
                    model, chain, segment, residue = current_path

                    # Warn the user
                    warnstr = \
                        f"Missing atoms {', '.join(missing_atoms)} " \
                        f"in model {model}, chain {chain}, segment " \
                        f"{segment}, residue {residue}."
                    logger.warning(warnstr)

                # For each provided atom name
                for atom_name in atom_names:

                    # For each atom found in the residue
                    for atom_serial, atom_data \
                    in expected_atoms.items():

                        # If the current atom's name matches
                        # the current name
                        if atom_data["atom_name"] == atom_name:

                            # Add it to the dictionary storing
                            # the sorted atoms
                            new_atoms_order[atom_serial] = atom_data

                # If we want to add the extra atoms before the
                # sorted atoms
                if other_atoms_position == "before":

                    # Create a new dictionary with the extra
                    # atoms first
                    new_atoms_order = \
                        {**other_atoms, **new_atoms_order}
                
                # If we want to add the extra atoms after the
                # sorted atoms
                elif other_atoms_position == "after":

                    # Create a new dictionary with the
                    # extra atoms after
                    new_atoms_order = \
                        {**new_atoms_order, **other_atoms}
                
                # If an invalid value was passed
                else:

                    # Raise an error
                    errstr = \
                        "'other_atoms_position' must be 'before' " \
                        "or 'after'"
                    raise ValueError(errstr)

                # Return the atoms in the new order
                return new_atoms_order
            
            # If we are at the residues' level
            elif current_depth == 4:

                # For each residue and associated data
                for res, res_data in struct.items():

                    # If the current residue's name
                    # matches the residue type for which
                    # we need to sort the atoms
                    if res_data["res_name"] == residue_name:

                        # Recurse through the atoms belonging to
                        # the current residue
                        sorted_atoms = \
                            recursive_step(\
                                struct = res_data["atoms"],
                                residue_name = residue_name,
                                atoms_names = atoms_names,
                                other_atoms_position = \
                                    other_atoms_position, 
                                current_depth = current_depth + 1,
                                current_path = current_path + (res,))
                        
                        # Add the sorted atoms to the current
                        # residue
                        res_data["atoms"] = sorted_atoms
            
            # Otherwise
            else:

                # For each key and associated sub-structure
                # in the current sub-structure
                for key, sub_struct in struct.items():

                    # If the sub-structure is a dictionary
                    if isinstance(sub_struct, dict):

                        # Recurse through it
                        recursive_step(\
                            struct = sub_struct,
                            residue_name = residue_name,
                            atoms_names = atoms_names,
                            other_atoms_position = \
                                 other_atoms_position,
                            current_depth = current_depth + 1,
                            current_path = current_path + (key,))

            # Return the updated structure
            return struct


        #---------------------- Sort the atoms -----------------------#


        # Sort the atoms (we need to pass a copy of the 'atom_data'
        # attribute because it will be modified)
        self._atom_data = \
            recursive_step(\
                struct = copy.deepcopy(self.atom_data),
                residue_name = residue_name,
                atoms_names = atoms_names,
                other_atoms_position = other_atoms_position)


        #-------------------- Renumber the atoms ---------------------#


        # Renumber the atoms (it also updates the CONECT data)
        self._renumber_atoms()


    #----------------------------- Copy ------------------------------#


    def get_copy(self):
        """Return a copy of the current structure.

        Returns
        -------
        ``struct_copy`` : ``pdbcraft.structure.Structure``
            A copy of the current structure.
        """
        
        return Structure(atom_data = copy.deepcopy(self.atom_data),
                         conect_data = copy.deepcopy(self.conect_data),
                         name = copy.deepcopy(self.name))


    #---------------------- Update CONECT data -----------------------#


    def add_bonds(self,
                  bonds,
                  in_place = False):
        """Add bonds between pairs of atoms. This updates the
        CONECT data associated with the structure.

        Parameters
        ----------
        bonds : ``list``
            A list of tuples, each representing a bond between two
            atoms and whether the bond is a ``"single"``, 
            ``"double"``, or ``"triple"`` bond.

            An atom can be represented either by its unique
            serial number or by its "path", namely a tuple
            containing the model, chain, segment, and residue
            the atom belongs to, plus the atom's name. This is
            based on the assumption that, within each residue,
            atoms must have unique names.

            For instance, a single bond between atom 1 and atom 2
            could be represented as ``(1, 2, "single")`` or as
            ``((1, "A", "", (1, " "), "N"), 
            (1, "A", "", (1, " "), "CA"),
            "single")``, assuming that both atoms belong to
            model 1, chain A, an unnamed segment, residue 1 with
            insertion code " " (= no insertion code), and one is
            named N, while the other is named CA. The model's
            specification is necessary because atoms in the same
            chain, segment, and residue may have different serial
            numbers in different models.

            Each bond should be specified only once (for instance,
            we should not specify a bond between atom 2 and atom 1
            if we have already specified a bond between atom 1 and
            atom 2) since, for each bond specified, the CONECT records
            for both atoms involved in the bond will be updated.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = {"bonds" : bonds}

        # If the structure needs to be modified in place
        if in_place:

            # Add the bonds to the structure
            self._add_bonds(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Add the bonds to the new structure
            struct_new.add_bonds(**kwargs)

            # Return the new structure
            return struct_new


    #----------------------------- Merge -----------------------------#


    def merge(self,
              other):
        """Merge the current structure with another ``Structure``,
        and return the merged structure.

        Parameters
        ----------
        other : ``pdbcraft.structure.Structure``
            The other structure.

        Returns
        -------
        merged : ``pdbcraft.structure.Structure``
            The merged structure.
        """

        # Merge the structures and return the resulting structure
        return self._merge(other = other)


    #--------------------- Rename/renumber items ---------------------#


    def renumber_models(self,
                        mapping = None,
                        start = None,
                        models = None,
                        in_place = False):
        """Renumber models in the structure.

        You can either provide a dictionary mapping the old
        numbering to the new numbering (``mapping``) or a
        starting number for the new numbering (``start``).

        You can choose to renumber only some models by using
        the ``models`` option.

        Parameters
        ----------
        mapping : ``dict``, optional
            The mapping between the old numbering
            and the new numbering.

        start : ``int``, optional
            The new starting number for the models.

        models : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of models to be
            renumbered. If not provided, all models
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "renumber",
             "item" : "model",
             "mapping" : "mapping",
             "start" : "start",
             "models" : models}

        # If the structure needs to be modified in place
        if in_place:

            # Renumber the models
            self._relabel(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Renumber the models in the new structure
            struct_new._relabel(**kwargs)

            # Return the new structure
            return struct_new


    def rename_chains(self,
                      mapping,
                      models = None,
                      in_place = False):
        """Rename chains in the structure.

        You must provide a dictionary mapping the old
        chain identifiers to the new identifiers
        (``mapping``). Only chains whose identifiers
        are in the dictionary will be renamed.

        You can choose to rename chains in only some
        models by using the ``models`` option.

        Parameters
        ----------
        mapping : ``dict``
            The mapping between the old identifiers
            and the new identifiers.

        models : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of models whose chains
            will be renamed. If not provided, all models
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "rename",
             "item" : "chain",
             "mapping" : mapping,
             "models" : models}

        # If the structure needs to be modified in place
        if in_place:

            # Rename the chains
            self._relabel(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Rename the chains in the new structure
            struct_new._relabel(**kwargs)

            # Return the new structure
            return struct_new


    def rename_segments(self,
                        mapping,
                        models = None,
                        chains = None,
                        in_place = False):
        """Rename the segments in the structure.

        You must provide a dictionary mapping the old
        segment identifiers to the new identifiers
        (``mapping``). Only segments whose identifiers
        are in the dictionary will be renamed.

        You can choose to rename segments in only some
        models or chains by using the ``models`` and
        ``chains`` options.

        Parameters
        ----------
        mapping : ``dict``
            The mapping between the old identifiers
            and the new identifiers.

        models : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of models whose segments
            will be renamed. If not provided, all models
            will be considered.

        chains : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of chains whose segments
            will be renamed. If not provided, all chains
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "rename",
             "item" : "segment",
             "mapping" : mapping,
             "models" : models,
             "chains" : chains}

        # If the structure needs to be modified in place
        if in_place:

            # Rename the segments
            self._relabel(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Rename the segments in the new structure
            struct_new._relabel(**kwargs)

            # Return the new structure
            return struct_new


    def renumber_residues(self,
                          mapping = None,
                          start = None,
                          models = None,
                          chains = None,
                          segments = None,
                          in_place = False):
        """Renumber residues in the structure.

        You can either provide a dictionary mapping the old
        numbering to the new numbering (``mapping``) or a
        starting number for the new numbering (``start``).

        The ``mapping`` must indicate residues with their
        sequence number (no insertion code), and all matching
        residues will be considered by the renumbering,
        regardless of their insertion code.

        You can choose to renumber only residues in specific
        models, chains, and segments by using the ``models``,
        ``chains``, and ``segments`` options.

        Parameters
        ----------
        mapping : ``dict``, optional
            The mapping between the old numbering
            and the new numbering.

        start : ``int``, optional
            The new starting number for the models.

        models : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of models whose residues
            will be renumbered. If not provided, all models
            will be considered.

        chains : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of chains whose residues
            will be renumbered. If not provided, all chains
            will be considered.

        segments : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of segments whose residues
            will be renumbered. If not provided, all segments
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "renumber",
             "item" : "residue",
             "mapping" : mapping,
             "start" : start,
             "models" : models,
             "chains" : chains,
             "segments" : segments}

        # If the structure needs to be modified in place
        if in_place:

            # Renumber the residues
            self._relabel(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Renumber the residues in the new structure
            struct_new._relabel(**kwargs)

            # Return the new structure
            return struct_new


    def rename_residues(self,
                        mapping,
                        models = None,
                        chains = None,
                        segments = None,
                        in_place = False):
        """Change the names (type) of residues.

        You must provide a dictionary mapping the
        residues' old names to the residues' new names
        (``mapping``).

        Parameters
        ----------
        mapping : ``dict``
            The mapping between the residues' old names
            to the residues' new names.

        models : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of models whose residues
            will be renamed. If not provided, all models
            will be considered.

        chains : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of chains whose residues
            will be renamed. If not provided, all chains
            will be considered.

        segments : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of segments whose residues
            will be renamed. If not provided, all segments
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "rename",
             "item" : "residue_attribute",
             "mapping" : mapping,
             "attribute" : "res_name",
             "models" : models,
             "chains" : chains,
             "segments" : segments}

        # If the structure needs to be modified in place
        if in_place:

            # Rename the residues
            self._relabel(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Rename the residues in the new structure
            struct_new._relabel(**kwargs)

            # Return the new structure
            return struct_new


    def rename_atoms(self,
                     mapping,
                     models = None,
                     chains = None,
                     segments = None,
                     residues = None,
                     residues_types = None,
                     in_place = False):
        """Change the names of atoms. The atoms' types
        remain invaried. This is useful especially to
        convert between nomenclatures for hydrogens
        used by different programs generating
        PDB files.

        You must provide a dictionary mapping the
        atoms' old names to the atoms' new names
        (``mapping``).

        Parameters
        ----------
        mapping : ``dict``
            The mapping between the atoms' new names
            and the atoms' new names.

        models : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of models whose atoms
            will be renamed. If not provided, all models
            will be considered.

        chains : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of chains whose atoms
            will be renamed. If not provided, all chains
            will be considered.

        segments : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of segments whose atoms
            will be renamed. If not provided, all segments
            will be considered.

        residues : ``list``, ``tuple``, or ``set``, optional
            A list, tuple, or set of residues whose atoms
            will be renamed. If not provided, all residues
            will be considered.

        residues_types : ``list``, ``tuple``, or ``set``, \
            optional
            A list, tuple, or set of residue types whose atoms
            will be renamed. If not provided, residues of all
            types will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "rename",
             "item" : "atom_attribute",
             "mapping" : mapping,
             "attribute" : "atom_name",
             "models" : models,
             "chains" : chains,
             "segments" : segments,
             "residues" : residues,
             "residues_types" : residues_types}

        # If the structure needs to be modified in place
        if in_place:

            # Rename the atoms
            self._relabel(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Rename the atoms in the new structure
            struct_new._relabel(**kwargs)

            # Return the new structure
            return struct_new


    #-------------------------- Keep items ---------------------------#


    def keep_models(self,
                    models,
                    in_place = False):
        """Keep only selected models.

        Parameters
        ----------
        models : ``list``
            The models to be kept.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "keep",
             "items" : models,
             "items_type" : "model"}

        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected models
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected models in the new structure
            self._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new

        # There is no need to renumber the atoms since
        # the atom numbering re-starts for each model


    def keep_chains(self,
                    chains,
                    models = None,
                    in_place = False):
        """Keep only selected chains.

        Parameters
        ----------
        chains : ``list``
            The chains to be kept.

        models : ``list``, optional
            The models considered when selecting the
            chains to be kept.

            The chains whose identifiers match the ones
            provided in ``chains`` will be removed only
            from the selected models.

            If no ``models`` are provided, all models
            will be considered.
        
        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "keep",
             "items" : chains,
             "items_type" : "chain",
             "models" : models}

        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected chains
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected chains in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def keep_segments(self,
                      segments,
                      models = None,
                      chains = None,
                      in_place = False):
        """Keep only selected segments.

        Parameters
        ----------
        segments : ``list``
            The segments to be kept.

        models : ``list``, optional
            The models considered when selecting the
            segments to be kept.

            The segments whose identifiers match the ones
            provided in ``segments`` will be removed only
            from the selected models.

            If no ``models`` are provided, all models
            will be considered.

        chains : ``list``, optional
            The chains considered when selecting the
            segments to be kept.
            
            The segments whose identifiers match the ones
            provided in ``segments`` will be removed only
            from the selected models.

            If no ``chains`` are provided, all chains
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "keep",
             "items" : segments,
             "items_type" : "segment",
             "models" : models,
             "chains" : chains}

        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected segments
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected segments in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def keep_residues(self,
                      residues,
                      models = None,
                      chains = None,
                      segments = None):
        """Keep only selected residues.

        Parameters
        ----------
        residues : ``list``
            The residues to be kept.

            It can be eiter a list of integers representing
            the sequence numbers of the residues to be removed,
            a list of tuples representing the sequence
            number/insertion code pairs of the residues to be
            removed, a list of strings representing the types
            of the residues to be removed, or a list with
            mixed integers, string, and tuples.

            For the integers in the list, residues whose
            sequence numbers match the integers will be
            kept regardless of their insertion codes.

        models : ``list``, optional
            The models considered when selecting the
            residues to be kept.

            The residues whose identifiers, sequence
            numbers, or types match the ones provided
            in ``residues`` will be removed only
            from the selected models.

            If no ``models`` are provided, all models
            will be considered.

        chains : ``list``, optional
            The chains considered when selecting the
            residues to be kept.
            
            The residues whose identifiers, sequence
            numbers, or types match the ones provided
            in ``residues`` will be removed only
            from the selected chains.

            If no ``chains`` are provided, all chains
            will be considered.

        segments : ``list``, optional
            The segments considered when selecting the
            residues to be kept.

            The residues whose identifiers, sequence
            numbers, or types match the ones provided
            in ``residues`` will be removed only
            from the selected segments.

            If no ``segments`` are provided, all segments
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "keep",
             "items" : residues,
             "items_type" : "residue",
             "models" : models,
             "chains" : chains,
             "segments" : segments}

        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected residues
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected residues in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def keep_atoms(self,
                   atoms,
                   models = None,
                   chains = None,
                   segments = None,
                   residues = None,
                   in_place = False):
        """Keep only selected atoms.

        Parameters
        ----------
        atoms : ``list``
            The atoms to be kept.

        models : ``list``, optional
            The models considered when selecting the
            atoms to be kept.

            The atoms whose identifiers match the ones
            provided in ``atoms`` will be kept only
            in the selected models.

            If no ``models`` are provided, all models
            will be considered.

        chains : ``list``, optional
            The chains considered when selecting the
            atoms to be kept.

            The atoms whose identifiers match the ones
            provided in ``atoms`` will be kept only
            in the selected chains.

            If no ``chains`` are provided, all chains
            will be considered.

        segments : ``list``, optional
            The segments considered when selecting the
            atoms to be kept.

            The atoms whose identifiers match the ones
            provided in ``atoms`` will be kept only
            in the selected segments.

            If no ``segments`` are provided, all segments
            will be considered.

        residues : ``list``, optional
            The residues considered when selecting
            the atoms to be kept.

            It can be eiter a list of integers representing
            the sequence numbers of the residues to be removed,
            a list of tuples representing the sequence
            number/insertion code pairs of the residues to be
            removed, a list of strings representing the types
            of the residues to be removed, or a list with
            mixed integers, string, and tuples.

            For the integers in the list, residues whose
            sequence numbers match the integers will be
            considered regardless of their insertion codes.

            If no ``residues`` are provided, all residues
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "keep",
             "items" : atoms,
             "items_type" : "atom",
             "models" : models,
             "chains" : chains,
             "segments" : segments,
             "residues" : residues}

        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected atoms
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected atoms in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def keep_protein(self,
                     extra_protein_residues = None,
                     in_place = False):
        """Keep only protein residues.

        Parameters
        ----------
        extra_protein_residues : ``list``, optional
            A list of residues that are not part of the canonical
            set but need to be considered protein residues.
        
        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Get the canonical set of protein residues
        protein_residues = self._RESNAMES_PROTEIN

        # If a list of extra residues was passed
        if extra_protein_residues is not None:

            # Add the residues to the set
            protein_residues.extend(extra_protein_residues)

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "keep",
             "items" : protein_residues,
             "items_type" : "residue"}

        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected residues
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected residues in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def keep_dna(self,
                 extra_dna_residues = None,
                 in_place = False):
        """Keep only DNA residues.

        Parameters
        ----------
        extra_dna_residues : ``list``, optional
            A list of residues that are not part of the canonical
            set but need to be considered DNA residues.
        
        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Get the canonical set of DNA residues
        dna_residues = self._RESNAMES_DNA

        # If a list of extra residues was passed
        if extra_dna_residues is not None:

            # Add the residues to the set
            dna_residues.extend(extra_dna_residues)

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "keep",
             "items" : dna_residues,
             "items_type" : "residue"}

        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected residues
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected residues in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def keep_rna(self,
                 extra_rna_residues = None,
                 in_place = False):
        """Keep only RNA residues.

        Parameters
        ----------
        extra_rna_residues : ``list``, optional
            A list of residues that are not part of the canonical
            set but need to be considered RNA residues.
        
        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Get the canonical set of RNA residues
        rna_residues = self._RESNAMES_RNA

        # If a list of extra residues was passed
        if extra_rna_residues is not None:

            # Add the residues to the set
            rna_residues.extend(extra_rna_residues)

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "keep",
             "items" : rna_residues,
             "items_type" : "residue"}

        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected residues
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected residues in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    #------------------------- Remove items --------------------------#


    def remove_models(self,
                      models,
                      in_place = False):
        """Remove selected models.

        Parameters
        ----------
        models : ``list``
            The models to be removed.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "remove",
             "items" : models,
             "items_type" : "model"}

        # If the structure needs to be modified in place
        if in_place:
            
            # Remove the selected models
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Remove the selected models in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new

        # There is no need to renumber the atoms since
        # the atom numbering re-starts for each model


    def remove_chains(self,
                      chains,
                      models = None,
                      in_place = False):
        """Remove selected models.

        Parameters
        ----------
        chains : ``list``
            The chains to be removed.

        models : ``list``, optional
            The models considered when selecting the
            chains to remove.

            The chains whose identifiers match the ones
            provided in ``chains`` will be removed only
            from the selected models.

            If no ``models`` are provided, all models
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "remove",
             "items" : chains,
             "items_type" : "chain",
             "models" : models}

        # If the structure needs to be modified in place
        if in_place:
            
            # Remove the selected chains
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Remove the selected chains in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def remove_segments(self,
                        segments,
                        models = None,
                        chains = None,
                        in_place = False):
        """Remove segments from the structure.

        Parameters
        ----------
        segments : ``list``
            The segments to be removed.

        models : ``list``, optional
            The models considered when selecting the
            segments to remove.

            The segments whose identifiers match the ones
            provided in ``segments`` will be removed only
            from the selected models.

            If no ``models`` are provided, all models
            will be considered.

        chains : ``list``, optional
            The chains considered when selecting the
            segments to remove.

            The segments whose identifiers match the ones
            provided in ``segments`` will be removed only
            from the selected chains.

            If no ``chains`` are provided, all chains
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "remove",
             "items" : segments,
             "items_type" : "segment",
             "models" : models,
             "chains" : chains}

        # If the structure needs to be modified in place
        if in_place:
            
            # Remove the selected segments
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Remove the selected segments in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def remove_residues(self,
                        residues,
                        models = None,
                        chains = None,
                        segments = None,
                        in_place = False):
        """Remove residues from the structure.

        Parameters
        ----------
        residues : ``list``
            The residues to be removed.

            It can be eiter a list of integers representing
            the sequence numbers of the residues to be removed,
            a list of tuples representing the sequence
            number/insertion code pairs of the residues to be
            removed, a list of strings representing the types
            of the residues to be removed, or a list with
            mixed integers, string, and tuples.

            For the integers in the list, residues whose
            sequence numbers match the integers will be
            removed regardless of their insertion codes.

        models : ``list``, optional
            The models considered when selecting the
            residues to remove.

            The residues whose identifiers, sequence
            numbers, or types match the ones provided
            in ``residues`` will be removed only
            from the selected models.

            If no ``models`` are provided, all models
            will be considered.

        chains : ``list``, optional
            The chains considered when selecting the
            residues to remove.

            The residues whose identifiers, sequence
            numbers, or types match the ones provided
            in ``residues`` will be removed only
            from the selected chains.

            If no ``chains`` are provided, all chains
            will be considered.

        segments : ``list``, optional
            The segments considered when selecting the
            residues to remove.

            The residues whose identifiers, sequence
            numbers, or types match the ones provided
            in ``residues`` will be removed only
            from the selected segments.

            If no ``segments`` are provided, all segments
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "remove",
             "items" : residues,
             "items_type" : "residue",
             "models" : models,
             "chains" : chains,
             "segments" : segments}

        # If the structure needs to be modified in place
        if in_place:
            
            # Remove the selected residues
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Remove the selected residues in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def remove_atoms(self,
                     atoms,
                     models = None,
                     chains = None,
                     segments = None,
                     residues = None,
                     in_place = False):
        """Remove atoms from the structure.

        Parameters
        ----------
        atoms : ``list``
            The atoms to be removed.

        models : ``list``, optional
            The models considered when selecting
            the atoms to remove.

            The atoms whose identifiers match the ones
            provided in ``atoms`` will be removed only
            from the selected models.

            If no ``models`` are provided, all models
            will be considered.

        chains : ``list``, optional
            The chains considered when selecting the
            atoms to remove.

            The atoms whose identifiers match the ones
            provided in ``atoms`` will be removed only
            from the selected chains.

            If no ``chains`` are provided, all chains
            will be considered.

        segments : ``list``, optional
            The segments considered when selecting
            the atoms to remove.

            The atoms whose identifiers match the ones
            provided in ``atoms`` will be removed only
            from the selected segments.

            If no ``segments`` are provided, all segments
            will be considered.

        residues : ``list``, optional
            The residues considered when selecting
            the atoms to remove.

            It can be eiter a list of integers representing
            the sequence numbers of the residues to be removed,
            a list of tuples representing the sequence
            number/insertion code pairs of the residues to be
            removed, a list of strings representing the types
            of the residues to be removed, or a list with
            mixed integers, string, and tuples.

            For the integers in the list, residues whose
            sequence numbers match the integers will be
            considered regardless of their insertion codes.

            If no ``residues`` are provided, all residues
            will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "remove",
             "items" : atoms,
             "items_type" : "atom",
             "models" : models,
             "chains" : chains,
             "segments" : segments,
             "residues" : residues}

        # If the structure needs to be modified in place
        if in_place:
            
            # Remove the selected atoms
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Remove the selected atoms in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def remove_protein(self,
                       extra_protein_residues = None,
                       in_place = False):
        """Remove protein residues.

        Parameters
        ----------
        extra_protein_residues : ``list``, optional
            A list of residues that are not part of the canonical
            set but need to be considered protein residues.
        
        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Get the canonical set of protein residues
        protein_residues = self._RESNAMES_PROTEIN

        # If a list of extra residues was passed
        if extra_protein_residues is not None:

            # Add the residues to the set
            protein_residues.extend(extra_protein_residues)

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "remove",
             "items" : protein_residues,
             "items_type" : "residue"}

        # If the structure needs to be modified in place
        if in_place:
            
            # Remove the selected residues
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Remove the selected residues in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def remove_dna(self,
                   extra_dna_residues = None,
                   in_place = False):
        """Remove DNA residues.

        Parameters
        ----------
        extra_dna_residues : ``list``, optional
            A list of residues that are not part of the canonical
            set but need to be considered DNA residues.
        
        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Get the canonical set of DNA residues
        dna_residues = self._RESNAMES_DNA

        # If a list of extra residues was passed
        if extra_dna_residues is not None:

            # Add the residues to the set
            dna_residues.extend(extra_dna_residues)

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "remove",
             "items" : dna_residues,
             "items_type" : "residue"}

        # If the structure needs to be modified in place
        if in_place:
            
            # Remove the selected residues
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Remove the selected residues in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    def remove_rna(self,
                   extra_rna_residues = None,
                   in_place = False):
        """Remove RNA residues.

        Parameters
        ----------
        extra_rna_residues : ``list``, optional
            A list of residues that are not part of the canonical
            set but need to be considered RNA residues.
        
        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Get the canonical set of RNA residues
        rna_residues = self._RESNAMES_RNA

        # If a list of extra residues was passed
        if extra_rna_residues is not None:

            # Add the residues to the set
            rna_residues.extend(extra_rna_residues)

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "remove",
             "items" : rna_residues,
             "items_type" : "residue"}

        # If the structure needs to be modified in place
        if in_place:
            
            # Remove the selected residues
            self._keep_or_remove(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Remove the selected residues in the new structure
            struct_new._keep_or_remove(**kwargs)

            # Return the new structure
            return struct_new


    #-------------------------- Move atoms ---------------------------#


    def move_atoms(self,
                   atoms,
                   new_chain = None,
                   new_segment = None,
                   new_residue = None,
                   in_place = False):
        """Move a set of atoms to another (possibly new) chain,
        segment, and/or residue.

        Parameters
        ----------
        atoms : ``list``
            A list containing the serial numbers of the atoms
            to be moved.

        new_chain : ``str``, optional
            The identifier of a new chain to move the atoms to.

            If no chain identifier is passed, the atoms are
            kept on the same chain they were before.

        new_segment : ``str``, optional
            The identifier of a new segment to move the atoms to.

            If no segment identifier is passed, the atoms are
            kept on the same segment they were before.

        new_residue : ``tuple``, optional
            A tuple containing the sequence number, the
            iCode, and the name of a new residue to move
            the atoms to.

            If no residue is passed, the atoms are kept on
            the same residue they were before.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Move the atoms
        return self._move_atoms(atoms = atoms,
                                create_new_struct = not in_place,
                                new_chain = new_chain,
                                new_segment = new_segment,
                                new_residue = new_residue,
                                in_place = in_place)


    #-------------------------- Sort atoms ---------------------------#


    def sort_atoms_residue_name(self,
                                atoms_names,
                                residue_name,
                                other_atoms_position = "after",
                                in_place = False):
        """Sort the atoms of all residues having a specific name
        according to a list of atoms' names.

        Parameters
        ----------
        atoms_names : ``list``
            A list of atoms' names defining the sorting order
            for the atoms.

        residue_name : ``str``
            The name of residue whose atoms will be sorted.

        other_atoms_position : ``str``, \
            {``"before"``, ``"after"``}, ``"after"``
            In case more atoms that the ones specified in the
            list are found in a residue, where to place these
            atoms relative to the sorted ones.

            ``"before"`` means that, in each residue whose
            atoms are sorted, the atoms whose names do not
            match those provided in ``atoms_names`` will
            be placed before the sorted atoms.

            ``"after"`` means that, in each residue whose
            atoms are sorted, the atoms whose names do not
            match those provided in ``atoms_names`` will
            be placed after the sorted atoms.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        ``pdbcraft.Structure`` if ``in_place = False``, ``None``
        otherwise
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"atoms_names" : atoms_names,
             "residue_name" : residue_name,
             "other_atoms_position" : other_atoms_position}

        # If the structure needs to be modified in place
        if in_place:
            
            # Sort the atoms
            self._sort_atoms_residue_name(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Sort the atoms in the new structure
            struct_new._sort_atoms_residue_name(**kwargs)

            # Return the new structure
            return struct_new


    #-------------------------- Write files --------------------------#


    def write_pdb_file(self,
                       pdb_file,
                       write_conect_data = True,
                       write_models_records = True):
        """Write out the structure to a PDB file.

        Parameters
        ----------
        pdb_file : ``str``
            The PDB file where the structure will be written.

        write_conect_data : ``bool``, ``True``
            Write the CONECT data associated with the structure,
            if any are present.

        write_models_records : ``bool``, ``True``
            Write the MODEL and ENDMDL records for each model found
            in the structure.

            You can set it to ``False`` only if you have one model
            in the structure. 
        """


        #----------------- Sections' format strings ------------------#


        # Get the format string for ATOM records       
        fmt_atom = \
            "ATOM  {:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}   " \
            "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          " \
            "{:>2s}{:2s}\n"

        # Get the format string for HETATM records
        fmt_hetatom = \
            "HETATM{:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}   " \
            "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          " \
            "{:>2s}{:2s}\n"

        # Get the format string for MODEL records
        fmt_model = "MODEL     {:>4}\n"

        # Get the format string for TER records
        fmt_ter = "TER   {:5d}      {:3s} {:1s}{:4d}\n"

        # Get the format string for ENDMDL records
        fmt_endmdl = "ENDMDL\n"


        #-------------------- Write the PDB file ---------------------#


        # Open the output PDB file
        with open(pdb_file, "w") as out:


            #------------- Write the atomic coordinates --------------#


            # For each model and associated data
            for i, (mod, mod_d) in enumerate(self.atom_data.items()):

                # If the models' records need to be written
                if write_models_records:

                    # Write out the corresponding MODEL record
                    out.write(fmt_model.format(mod))

                # Otherwise
                else:

                    # If there are multiple models in the structure
                    if len(self.atom_data) > 1:

                        # If we are at the first model
                        if i == 0:

                            # Warn the user
                            warnstr = \
                                "Even though 'write_models_header' " \
                                "is set to 'False', we need to " \
                                "write the models' records (MODEL " \
                                "and ENDMDL) because there are " \
                                "multiple models in the structure."
                            logger.warning(warnstr)

                        # Write out the records anyway
                        out.write(fmt_model.format(mod))

                # For each chain and associated data
                for ch, ch_d in mod_d.items():

                    # For each segment and associated data
                    for seg, seg_d in ch_d.items():

                        # For each residue and associated data
                        for res, res_d in seg_d.items():

                            # For each atom and associated data
                            for atom, atom_d in res_d["atoms"].items():

                                # If we need to write a HETATM record
                                if atom_d["is_hetatm"]:
                                    
                                    # Select the corresponding format
                                    # string
                                    fmt_str = fmt_hetatom

                                # If we need to write an ATOM record
                                else:

                                    # Select the corresponding format
                                    # string
                                    fmt_str = fmt_atom

                                # Write out the ATOM/HETATM record
                                out.write(\
                                    fmt_str.format(\
                                        atom,
                                        atom_d["atom_name"],
                                        atom_d["alt_loc"],
                                        res_d["res_name"],
                                        ch,
                                        res[0],
                                        res[1],
                                        atom_d["x"],
                                        atom_d["y"],
                                        atom_d["z"],
                                        atom_d["occupancy"],
                                        atom_d["temp_factor"],
                                        atom_d["element"],
                                        atom_d["charge"]))

                                # Save the current atom serial to write
                                # the TER record for the current chain
                                # (plus one since the TER record has a
                                # serial one unit higher than the last
                                # atom of the chain)
                                ter_atom_serial = atom + 1

                            # Save the current residue number to
                            # write the TER record for the current
                            # chain
                            ter_res_seq = res[0]

                            # Save the current residue name to
                            # write the TER record for the current
                            # chain
                            ter_res_name = res_d["res_name"]

                    # Save the current chain ID to write the
                    # TER record for the current chain
                    ter_chain_id = ch

                    # Write out the TER record for the current
                    # chain
                    out.write(fmt_ter.format(ter_atom_serial,
                                             ter_res_name,
                                             ter_chain_id,
                                             ter_res_seq))

                # If the models' records need to be written
                if write_models_records:
                    
                    # Write out the ENDMDL record for the current
                    # model
                    out.write(fmt_endmdl)

                # Otherwise
                else:

                    # If there are multiple models in the structure
                    if len(self.atom_data) > 1:

                        # Write them out anyway (do not warn the
                        # user here since they have already been
                        # warned before)
                        out.write(fmt_model.format(mod))


            #----------------- Write the CONECT data -----------------#


            # If the user requested the writing of CONECT data
            if write_conect_data:

                # If there are conect data
                if self.conect_data is not None:

                    # For each model
                    for mod in self.conect_data:

                        # For each CONECT record in the model
                        for atom, rec in self.conect_data[mod].items():

                            # Write out the header
                            out.write("CONECT")

                            # For each atom in the CONECT record
                            for a in [atom, *rec]:

                                # Write out the atom
                                out.write("{:5d}".format(a))

                            # Write out a newline character at the
                            # end of the record
                            out.write("\n")

                # Otherwise
                else:

                    # Warn the user that there are no CONECT data
                    warnstr = \
                        "No CONECT data found for the structure."
                    logger.warning(warnstr)


    def write_fasta_file(self,
                         fasta_file,
                         split_models = False,
                         split_chains = False,
                         disc_chains_mode = "join_with_gaps",
                         res_i_code = " ",
                         gap_char = "-",
                         wrap_at = None,
                         resnames_3to1 = None):
        """Write a FASTA file containing the sequence corresponding
        to the structure.

        Parameters
        ----------
        fasta_file : ``str``
            The FASTA file to be written.

            If ``split_models``, or ``split_chains`` are
            ``True``, the name of the file is used as a
            prefix for the multiple FASTA files that will
            be written.

        split_models : ``bool``, ``True``
            If ``True``, one FASTA file per model will be
            written.

            Please refer to the Notes section below for a
            more detailed explanation of how the different
            combinations of ``split_models`` and ``split_chains``
            affect the content of the output FASTA files.

        split_chains : ``bool``, ``True``
            If ``True``, one FASTA file per chain will be
            written.

            Please refer to the Notes section below for a
            more detailed explanation of how the different
            combinations of ``split_models`` and ``split_chains``
            affect the content of the output FASTA files.

        disc_chains_mode : ``str``, {``"join"``, \
            ``"join_with_caps"``, ``"split"``}, ``"join_with_caps"``
            How to represent discontinuous chains in the FASTA
            sequences.

            A discontinuous chain is defined as a chain
            with a discontinuous numbering of its residues.

            If ``"join"``, the discontinuous portions of the
            chain will be part of the same FASTA sequence,
            without any indication of the gaps between them.

            If ``"join_with_gaps"``, the discontinuous portions
            of the chain will be part of the same FASTA sequence,
            with as many instances of a character indicating
            a gap (``gap_char``) as there are residues missing
            between two portions of the chain.

            If ``"split"``, the discontinuous portions of
            the chain will be written as separate FASTA
            sequences.

        res_i_code : ``str``, ``" "``
            The insertion code of the residues that will be
            included in the FASTA sequences.

            If not provided, the function will consider only
            residues without an insertion code.

        gap_char : ``str``, ``"-"``
            The character used to represent gaps in the FASTA
            sequences, if ``disc_chains_mode = "join_with_caps"``.

        wrap_at : ``int``, optional
            Wrap each FASTA sequence at ``wrap_at`` characters.
            
            If it is not provided, do not wrap the sequences.

        resnames_3to1 : ``dict``, optional
            A dictionary containing the mapping between the
            residues' three-letter names and their one-letter
            name.

            It is used when writing FASTA files from the
            structure.

            There is no need to pass it if your structure only
            contains the 20 canonical protein residues, since
            their mapping is hard-coded. However, if an
            updated mapping for these residues is provided,
            the new one is used when writing FASTA files.

        Notes
        -----
        If ``split_models = False`` and ``split_chains = False``,
        the FASTA sequences of all chains of all models will
        be written in one FASTA file.
        The file will have the name provided in the ``fasta_file``
        option.

        If ``split_models = False`` and ``split_chains = True``,
        one FASTA file per chain will be written, containing
        the FASTA sequence of that chain in all models.
        Assuming the ``fasta_file`` option has the form 
        ``{file_name}.{extension}``, each file will be named
        ``{file_name}_{chain}.{extension}``.

        If ``split_models = True`` and ``split_chains = False``,
        one FASTA file per model will be written, containing
        the FASTA sequences of the chains in that model.
        Assuming the ``fasta_file`` option has the form 
        ``{file_name}.{extension}``, each file will be named
        ``{file_name}_{model}.{extension}``.

        If ``split_models = True`` and ``split_chains = True``,
        one FASTA file per chain in each model will be written,
        containing the FASTA sequence of that chain in that
        model.
        Assuming the ``fasta_file`` option has the form 
        ``{file_name}.{extension}``, each file will be named
        ``{file_name}_{model}_{chain}.{extension}``.
        """


        # Create an empty set to store the files that have
        # been created
        files_created = set()


        #--------------------- FASTA file's name ---------------------#


        # Get the file's name (including the path leading to
        # it, if any) and the file's extension, separately
        fasta_name, fasta_ext = \
            os.path.splitext(fasta_file)


        #------------------- FASTA headers' prefix -------------------#


        # If the structure has no name
        if self.name is None:

            # Set the prefix for the FASTA headers to
            # an empty string
            prefix = ""

        # Otherwise
        else:

            # The prefix will contain the name of the structure
            prefix = f"{self.name}_"


        #------------------ Residues' names mapping ------------------#


        # If no mapping between the residues' three-letter
        # names and their one-letter names was passed
        if resnames_3to1 is None:

            # Set it to the default one
            resnames_3to1 = self._RESNAMES_3TO1


        #----------------- Open file - do not split ------------------#


        # If we do not need to split neither models nor chains
        if not split_models and not split_chains:

            # Set the name of the only FASTA file that will be
            # written
            file_name = fasta_file

            # If the file has beem already created
            if file_name in files_created:

                # Open the file in 'append' mode
                file_mode = "a"

            # If the file has not been already created
            else:

                # Open the file in 'write' mode
                file_mode = "w"

                # Add the name of the file to the set of
                # files created
                files_created.add(file_name)

            # Open the file in the selected mode
            out = open(file_name, file_mode)


        #-------------------------- Models ---------------------------#


        # For each model and associated data in the
        # structure
        for mod, mod_data in self.atom_data.items():


            #--------------- Open file - split models ----------------#


            # If we need to split models but not chains
            if split_models and not split_chains:

                # Set the name of the FASTA file for the
                # current model
                file_name = f"{fasta_name}_{mod}{fasta_ext}"

                # If the file has beem already created
                if file_name in files_created:

                    # Open the file in 'append' mode
                    file_mode = "a"

                # If the file has not been already created
                else:

                    # Open the file in 'write' mode
                    file_mode = "w"

                    # Add the name of the file to the set of
                    # files created
                    files_created.add(file_name)

                # Open the file in the selected mode
                out = open(file_name, file_mode)


            #------------------------ Chains -------------------------#

            
            # For each chain and associated data
            # in the current model
            for ch, ch_data in mod_data.items():

                # Inizialite a list to store the FASTA
                # sequences (and corresponding headers)
                # for the current chain in the current
                # model
                fasta_seqs = \
                    [[f">{prefix}{mod}_{ch}_", ""]]


                #--------------------- Segments ----------------------#


                # For each chain segment and associated data
                for seg, seg_data in ch_data.items():

                    # Initialize the variable storing the
                    # previously visited residue's sequence
                    # number to None
                    prev_res_num = None


                    #------------------- Residues --------------------#


                    # For each residue in the segment
                    for i, res in enumerate(seg_data):

                        # Get the residue's sequence number and
                        # insertion code
                        res_num, i_code = res

                        # If the insertion code is not the one
                        # of the residues to be considered
                        if i_code != res_i_code:

                            # Skip the current residue
                            continue

                        # Get the residue's three-letter name
                        res_name_3 = seg_data[res]["res_name"]

                        # Try to get the residue's one-letter name
                        # from its three-letter name
                        try:

                            res_name_1 = resnames_3to1[res_name_3]

                        # If no one-letter name was found for the
                        # three-letter name
                        except KeyError:

                            # Re-raise the error with a more
                            # verbose message
                            errstr = \
                                "No one-letter name found for " \
                                f"residues of type '{res_name_3}'."
                            raise KeyError(errstr)

                        # If we are at the first residue of the
                        # segment
                        if i == 0:

                            # Update the current chain's sequence
                            fasta_seqs[-1][1] += \
                                res_name_1

                            # Update the current chain's header
                            fasta_seqs[-1][0] += \
                                f"{res_num}"
                            
                            # Set the previous residue to the
                            # current residue
                            prev_res_num = res_num
                            
                            # Go to the next residue
                            continue

                        # If we reached a discontinuous point of the
                        # chain
                        if res_num > (prev_res_num + 1):

                            # If we should join discontinuous portions
                            # of the chain without gaps
                            if disc_chains_mode == "join":
                                
                                # Update the current chain's sequence
                                fasta_seqs[-1][1] += \
                                    res_name_1
                                
                                # Update the current chain's header
                                # by adding the end of the previous
                                # chain and the beginning of the
                                # current one
                                fasta_seqs[-1][0] += \
                                    f"-{prev_res_num}_{res_num}"

                            # If we should join discontinuous portions
                            # of the chain using gaps
                            elif disc_chains_mode == "join_with_gaps":

                                # Get the appropriate number of gaps
                                gaps = \
                                    gap_char * (res_num - prev_res_num)

                                # Update the current chain's sequence
                                fasta_seqs[-1][1] += \
                                    gaps + res_name_1

                                # Update the current chain's header
                                # by adding the end of the previous
                                # chain and the beginning of the
                                # current one
                                fasta_seqs[-1][0] += \
                                    f"-{prev_res_num}_{res_num}"

                            # If we should split discontinuous portions
                            # of the chain
                            elif disc_chains_mode == "split":

                                # Update the current chain's sequence
                                fasta_seqs.append(\
                                    [f">{self.name}_{res_num}",
                                     res_name_1])

                                # Update the previous chain's header
                                # by adding the end of the previous
                                # chain
                                fasta_seqs[-2][0] += \
                                    f"-{prev_res_num}"

                        # If we are at a continuous point in the chain
                        else:

                            # Update the current chain's sequence
                            fasta_seqs[-1][1] += \
                                res_name_1

                        # If we are at the last residue of the
                        # segment
                        if i == len(seg_data) - 1:

                            # Update the current chain's header
                            # by adding the end of the current
                            # chain
                            fasta_seqs[-1][0] += \
                                f"_{res_num}"

                        # Update the 'previous residue number'
                        # variable
                        prev_res_num = res_num


                #--------- Open file - split models + chains ---------#


                # If we need to split chains and not models
                if split_models and split_chains:

                    file_name = f"{fasta_name}_{mod}_{ch}{fasta_ext}"

                    # If the file has beem already created
                    if file_name in files_created:

                        # Open the file in 'append' mode
                        file_mode = "a"

                    # If the file has not been already created
                    else:

                        # Open the file in 'write' mode
                        file_mode = "w"

                        # Add the name of the file to the set of
                        # files created
                        files_created.add(file_name)

                    # Open the file in the selected mode
                    out = open(file_name, file_mode)


                #------------- Open file - split chains --------------#


                # If we need to split chains but not models
                elif not split_models and split_chains:

                    file_name = f"{fasta_name}_{ch}{fasta_ext}"

                    # If the file has beem already created
                    if file_name in files_created:

                        # Open the file in 'append' mode
                        file_mode = "a"

                    # If the file has not been already created
                    else:

                        # Open the file in 'write' mode
                        file_mode = "w"

                        # Add the name of the file to the set of
                        # files created
                        files_created.add(file_name)

                    # Open the file in the selected mode
                    out = open(file_name, file_mode)


                #------------- Write header and sequence -------------#


                # For each header and associated sequence
                for fasta_head, fasta_seq in fasta_seqs:

                    # Write out the header
                    out.write(f"{fasta_head}\n")

                    # If no line wrapping was requested
                    if wrap_at is None:
                        
                        # Write out the entire sequence
                        out.write(f"{fasta_seq}\n")

                    # Otherwise
                    else:

                        # Get 'wrap_at'-sized chunks of the
                        # sequence
                        fasta_seq_chunks = \
                            [fasta_seq[i:i + wrap_at] \
                            for i \
                            in range(0, len(fasta_seq), wrap_at)]

                        # For each chunk
                        for fasta_seq_chunk in fasta_seq_chunks:

                            # Write out the chunk
                            out.write(f"{fasta_seq_chunk}\n")
