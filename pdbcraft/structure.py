#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    structure.py
#
#    Utilities to manipulate atomic structures.
#
#    Copyright (C) 2024 Valentina Sora 
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


#######################################################################


# Set the module's description.
__doc__ = "Utilities to manipulate atomic structures."


#######################################################################


# Import from the standard library.
from collections import defaultdict
import copy
import itertools
import logging as log
import os
import re
# Import from 'pdbcraft'.
from . import _defaults


#######################################################################


# Get the module's logger.
logger = log.getLogger(__name__)


#######################################################################


class Structure:

    """
    Class implementing a molecular structure, usually parsed from
    either a PDB or a mmCIF file.
    """


    ###################################################################


    # Set the ['_items', '_attributes'] depth for each level of the
    # hierarchy.
    _ITEMSATTRS_DEPTH_TO_LEVELS = \
        {# Models' level - ['_items', '_attributes']
         2: "models",
         # Chains' level - ['_items', '_attributes']
         4: "chains",
         # Segments' level - ['_items', '_attributes']
         6: "segments",
         # Residues' level - ['_items', '_attributes']
         8: "residues",
         # Atoms' level - ['_items', '_attributes']
         10: "atoms"}

    # Set the identifiers' depth for each level of the hierarchy.
    _IDS_DEPTH_TO_LEVELS = \
        {# Models' level
         1 : "models",
         # Chains' level
         3 : "chains",
         # Segments' level
         5 : "segments", 
         # Residues' level
         7 : "residues",
         # Atoms' level
         9 : "atoms"}


    ###################################################################
    

    def __init__(self,
                 atom_data = None,
                 conect_data = None,
                 name = None):
        """Create a structure.

        Parameters
        ----------
        atom_data : ``dict``, optional
            A nested dictionary containing the atomic coordinates
            of the structure.

        conect_data : ``dict``, optional
            A nested dictionary containing the connectivity data
            for the structure.

        name : ``str``, optional
            The name of the structure.

        Notes
        -----
        To create an empty structure, leave both ``atom_data`` and
        ``conect_data`` blank.
        """

        # If no atomic coordinates were provided
        if atom_data is None:

            # Set the atomic coordinates to an empty dictionary.
            self._atom_data = {}

        # Otherwise
        else:
            
            # Set the atomic coordinates to the dictionary provided.
            self._atom_data = atom_data

        #-------------------------------------------------------------#

        # If no connectivity data were provided
        if conect_data is None:

            # If no atoms' coordinates were provided
            if atom_data is None:

                # Set the connectivity data to an empty dictionary.
                self._conect_data = {}

            # Otherwise
            else:

                # Set the connectivity data to a dictionary containing
                # an empty dictionary for each of the models found in
                # the atomic coordinates' dictionary.
                self._conect_data = \
                    {mod : {} for mod in self.atom_data}

        # Otherwise
        else:

            # Set the connectivity data to the dictionary provided.
            self._conect_data = conect_data

        #-------------------------------------------------------------#

        # Set the name of the structure.
        self._name = name

    
    ###################################################################


    @property
    def atom_data(self):
        """A dictionary containing the atomic coordinates.
        """
        
        return self._atom_data


    @atom_data.setter
    def atom_data(self,
                  value):
        """Raise an error if the user tries to modify the 'atom_data'
        attribute.
        """

        errstr = \
            "The 'atom_data' attribute cannot be directly modified."
        raise ValueError(errstr)


    @property
    def conect_data(self):
        """A dictionary containing the connectivity data associated
        with the structure.
        """
        
        return self._conect_data


    @conect_data.setter
    def conect_data(self,
                    value):
        """Raise an error if the user tries to modify the 'conect_data'
        attribute.
        """

        errstr = \
            "The 'conect_data' attribute cannot be directly modified."
        raise ValueError(errstr)


    @property
    def name(self):
        """The name of the structure.
        """
        
        return self._name

    
    ###################################################################


    # Define a helper function for the '__getitem__' method.
    def _recursive_getitem(self,
                           struct,
                           path):

        # Define the recursion as an inner function.
        def recurse(struct,
                    path):

            # If we reached the end of the path
            if len(path) == 1:

                # Try to return the last key in the path and the
                # associated value.
                try:
                    
                    return {path[0][0] : struct[path[0][0]]}

                # If there is no such key in the structure
                except KeyError:

                    # Raise an error.
                    errstr = \
                        f"No {path[0][1]} {path[0][0]} " \
                        "found in the structure."
                    raise KeyError(errstr)

            #---------------------------------------------------------#

            # Otherwise
            else:

                # For each key and associated sub-structure in the
                # current structure
                for key, sub_struct in list(struct.items()):

                    # If we are at an ['_items', '_attributes'] level
                    if key in ("_items", "_attributes"):

                        # If we are at the '_items'
                        if key == "_items":

                            # Recursively get items through the
                            # structure.
                            struct[key] = \
                                recurse(struct = sub_struct,
                                        path = path)

                    # Otherwise
                    else:

                        # If the current key is in the path
                        if key == path[0][0]:

                            # Recursively get items through the
                            # structure.
                            struct[key]["_items"] = \
                                recurse(struct = sub_struct["_items"],
                                        path = path[1:])

                        # Otherwise
                        else:

                            # Try to remove the key and associated
                            # sub-structure from the structure.
                            try:
                                
                                struct.pop(key)
                            
                            # If the key was not found
                            except KeyError:

                                # Raise an error.
                                errstr = \
                                    f"No {path[0][1]} {path[0][0]} " \
                                    "found in the structure."
                                raise KeyError(errstr)

            # Return the structure.
            return struct

        #-------------------------------------------------------------#

        # Return the result of the recursion.
        return recurse(struct = struct,
                       path = path)

    
    def __eq__(self,
               other):
        """Return ``True`` if all attributes of ``self`` are equal
        to the attributes or ``other`` (the structures' names can be
        different).

        Parameters
        ----------
        other: :class:`pdbcraft.structure.Structure`
            Another structure.
        """

        # If 'other' is not a 'Structure'
        if not isinstance(other, self.__class__):

            # Raise an error.
            errstr = \
                "Equality and inequality operations are " \
                "supported only between " \
                f"'{self.__class__.__name__}' instances."
            raise NotImplementedError(errstr)

        # Check the equality between the two structures.
        return (self.atom_data == other.atom_data \
                and self.conect_data == other.conect_data)


    def __ne__(self,
               other):
        """Return ``True`` if any of the data attributes of ``self``
        is different from the corresponding attribute of ``other``
        (the structures' names can be different).

        Parameters
        ----------
        other: :class:`pdbcraft.structure.Structure`
            Another structure.
        """

        # Check the inequality and return the result.
        return not self.__eq__(other = other)


    def __getitem__(self,
                    path):
        """Return a new structure with the items in a 'path',
        where a 'path' is a tuple containing the identifiers of
        the items.
        """

        # Set a tuple containing the different items' types
        # in hierarchical order (it is needed to provided useful
        # error messages if some item is not found in the
        # structure)
        item_types = ("model", "chain", "segment", "residue", "atom")

        #-------------------------------------------------------------#

        # Pair the items' types with the corresponding items in the
        # path.
        path = \
            tuple([(item, item_type) for item, item_type \
                   in zip(path, item_types[:len(path)])])

        #-------------------------------------------------------------#

        # Get the 'atom_data' for the new structure by keeping only the
        # items in the path.
        new_atom_data = \
            self._recursive_getitem(\
                struct = copy.deepcopy(self.atom_data),
                path = path)

        #-------------------------------------------------------------#

        # Build a new structure (for now, just copy the old structure's
        # connectivity data).
        new_struct = \
            self.__class__(atom_data = new_atom_data,
                           conect_data = \
                                copy.deepcopy(self.conect_data))

        #-------------------------------------------------------------#

        # Update the atom numbering and the connectivity data for the
        # new structure.
        new_struct._update_atom_and_conect()

        #-------------------------------------------------------------#

        # Return the new structure.
        return new_struct


    ###################################################################


    def _check_items(self,
                     model = None,
                     chain = None,
                     segment = None,
                     residue = None):
        """Check the items passed when using the following methods:

        * :meth:`pdbcraft.structure.Structure.chains`
        * :meth:`pdbcraft.structure.Structure.segments`
        * :meth:`pdbcraft.structure.Structure.residues`
        * :meth:`pdbcraft.structure.Structure.atoms`

        Parameters
        ----------
        model : ``int``, optional
            The model.

        chain : ``str``, optional
            The chain.

        segment : ``str``, optional
            The segment.

        residue : ``tuple``, optional
            The residue.
        """

        # If a model was passed
        if model is not None:

            # If the model is not an integer
            if not isinstance(model, int):

                # Raise an error.
                errstr = "'model' must be an integer."
                raise ValueError(errstr)

        #-------------------------------------------------------------#

        # If a chain was passed
        if chain is not None:

            # If the chain is not a string
            if not isinstance(chain, str):

                # Raise an error.
                errstr = "'chain' must be a string."
                raise ValueError(errstr)

        #-------------------------------------------------------------#

        # If a segment was passed
        if segment is not None:

            # If the segment is not a string
            if not isinstance(segment, str):

                # Raise an error.
                errstr = "'segment' must be a string."
                raise ValueError(errstr)

        #-------------------------------------------------------------#

        # If a residue was passed
        if residue is not None:

            # If the residue is not a tuple
            if not isinstance(residue, tuple):

                # Raise an error.
                errstr = "'residue' must be a tuple."
                raise ValueError(errstr)

            # If the tuple is not three elements long
            if len(residue) != 3:

                # Raise an error.
                errstr = \
                    "The 'residue' tuple must contain three " \
                    "elements: the residue's sequence number " \
                    "(int), the residue's insertion code (str), " \
                    "and the residue's name (str)."
                raise ValueError(errstr)

            # If the tuple's first element is not an integer
            if not isinstance(residue[0], int):

                # Raise an error.
                errstr = \
                    "The first element of the 'residue' tuple " \
                    "must be an integer (the residue's sequence " \
                    "number)."
                raise ValueError(errstr)

            # If the tuple's second element is not an integer
            if not isinstance(residue[1], str):

                # Raise an error.
                errstr = \
                    "The second element of the 'residue' tuple " \
                    "must be a string (the residue's insertion " \
                    "code)."
                raise ValueError(errstr)

            # If the tuple's third element is not an integer
            if not isinstance(residue[2], str):

                # Raise an error.
                errstr = \
                    "The third element of the 'residue' tuple " \
                    "must be a string (the residue's name)."
                raise ValueError(errstr)


    def _should_recurse(self,
                        key,
                        value,
                        selected,
                        current_depth):
        """Return whether we should keep on recursing in the
        'atom_data' dictionary by checking whether we are at the
        right model, chain, segment, or residue.

        Parameters
        ----------
        key : ``int``, ``str``, ``tuple``
            The key we have to check to see if we are at the correct
            model, chain, segment, or residue.

        value : ``dict``
            The dictionary we have to check to see if the attributes
            of the model, chain, segment, or residue have any
            of the accepted values.

        selected : ``dict``
            A dictionary containing the models, chains, segments, and
            residues through which we should recurse.

            It also contains the attributes of the selected models,
            chains, segments, and residues should have (if they exist)
            for us to recurse through them.

        current_depth : ``int``
            The depth we are currently at in the recursion.

        Returns
        -------
        should_recurse : ``bool``
            Whether we should recurse through the current sub-structure
            of the 'atom_data' dictionary.
        """

        # If we are at a level where we should check if we should
        # recurse
        if current_depth in self._IDS_DEPTH_TO_LEVELS:

            #---------------------------------------------------------#

            # Get the type of items found at this depth from the
            # depth itself.
            items_type = self._IDS_DEPTH_TO_LEVELS[current_depth]
            
            # Get the user-selected items for this item type.
            selected_items = selected[items_type]["_items"]

            # Get the user-selected attributes for this item type.
            selected_attrs = selected[items_type]["_attributes"]

            #---------------------------------------------------------#

            # If there are user-selected items and the current one is
            # not one of them
            if selected_items:

                # If the key is a tuple (it is a residue)
                if isinstance(key, tuple):

                    # If neither the residue's sequence number nor the
                    # residue's name nor the full residue's identifier
                    # are among the user-selected items
                    if key[0] not in selected_items \
                    and key[2] not in selected_items \
                    and key not in selected_items:

                        # Do not recurse.
                        return False
                
                # Otherwise 
                else:

                    # If the entire key is not among the selected items
                    if key not in selected_items:

                        # Do not recurse.
                        return False

            #---------------------------------------------------------#

            # If there are user-selected attributes
            if selected_attrs and isinstance(selected_attrs, dict):

                # For each attribute and associated values
                for attr_name, attr_value in selected_attrs.items():
                    
                    # If the attribute exists in the structure and its
                    # value does not match any of the accepted ones            
                    if attr_name in value["_attributes"] \
                    and value["_attributes"][attr_name] \
                        not in selected_attrs[attr_name]:

                        # Do not recurse.
                        return False

            #---------------------------------------------------------#

            # If we arrived here without returning, recurse.
            return True

        #-------------------------------------------------------------#

        # If we are not at a level where we can check, recurse by
        # default.
        return True


    def _get_selected(self,
                      models = None,
                      models_attributes = None,
                      chains = None,
                      chains_attributes = None,
                      segments = None,
                      segments_attributes = None,
                      residues = None,
                      residues_attributes = None,
                      atoms = None,
                      atoms_attributes = None,
                      extra_residues = None):
        """Get the items to be considered in an operation.

        Parameters
        ----------
        models : ``list``, optional
            The models to be considered.

        models_attributes : ``dict``, optional
            A dictionary of attributes and their accepted values.
            
            Only models whose values for the specified attributes
            (if they exist) are among the accepted values will be
            considered.

        chains : ``list``, optional
            The chains to be considered.

        chains_attributes : ``dict``, optional
            A dictionary of attributes and their accepted values.
            
            Only chains whose values for the specified attributes
            (if they exist) are among the accepted values will be
            considered.

        segments : ``list``, optional
            The segments to be considered.

        segments_attributes : ``dict``, optional
            A dictionary of attributes and their accepted values.
            
            Only segments whose values for the specified attributes
            (if they exist) are among the accepted values will be
            considered.

        residues : ``list``, optional
            The residues to be considered.

        residues_attributes : ``dict``, optional
            A dictionary of attributes and their accepted values.
            
            Only residues whose values for the specified attributes
            (if they exist) are among the accepted values will be
            considered.

        atoms : ``list``, optional
            The atoms to be considered.

        atoms_attributes : ``dict``, optional
            A dictionary of attributes and their accepted values.
            
            Only atoms whose values for the specified attributes
            (if they exist) are among the accepted values will be
            considered.

        extra_residues : ``dict``, optional
            A dictionary containing ``"protein"``, ``"dna"``, 
            ``"rna"``, or ``"het"`` residues to be considered part
            of the canonical sets for these entities, despite not
            being part of the sets.
        """

        # If a list of residues was specified
        if residues is not None:

            # For each possible group of residues
            for group in _defaults.STRUCT_RESNAMES:

                # If the user selected the group
                if group in residues:

                    # Get the canonical set of residues' names
                    # for that group and add it to the set of
                    # accepted residues.
                    residues.extend(\
                        _defaults.STRUCT_RESNAMES[group])

                    # If a list of extra residues was passed
                    if group in extra_residues:

                        # Add the residues to the set.
                        residues.extend(\
                            extra_residues[group])

                    # Remove the name of the group from the list of
                    # residues.
                    residues.remove(group)

            #---------------------------------------------------------#

            # If the user requested heteoresidues
            if "het" in residues:

                # If there are no residues' attributes yet
                if residues_attributes is None:
                    
                    # Create the residues' attributes' dictionary.
                    residues_attributes = {"is_het" : [True]}

                # Otherwise
                else:

                    # If the attribute representing atoms' names is
                    # already in the dictionary of attributes
                    if "is_het" in residues_attributes:

                        # Warn the user.
                        warnstr = \
                            "A list of accepted values for whether " \
                            "a residue is a heteroresidue or not " \
                            "was passed in the dictionary of " \
                            "attributes ('is_het'). Since the 'het' " \
                            "keyword was also passed in the list of " \
                            "accepted residues, 'True' will be " \
                            "added to the list of accepted values " \
                            "for the 'is_het' attribute."
                        logger.warning(warnstr)

                        # Extend the accepted values for the attribute.
                        residues_attributes["is_het"].extend([True])

                    # Otherwise
                    else:

                        # Add the element to the residues' attributes.
                        residues_attributes["is_het"] = [True]

                # Remove the element from the list of residues.
                residues.remove("het")

        #-------------------------------------------------------------#

        # If a list of atoms was specified
        if atoms is not None:

            # Assume that every string element in the list is the name
            # of an atom.
            atoms_names = \
                [elem for elem in atoms if isinstance(elem, str)]

            # Assume that every integer in the list is the serial
            # number of an atom.
            atoms = \
                [elem for elem in atoms if isinstance(elem, int)]

            #---------------------------------------------------------#

            # If any names were found
            if atoms_names:

                # If there are no atoms' attributes yet
                if atoms_attributes is None:
                    
                    # Create a dictionary containing the attributes
                    # to be considered and add the attribute
                    # representing the atoms' names to it.
                    #
                    # The accepted values for this attribute will be
                    # the user-provided names.
                    atoms_attributes = {"label_atom_id" : atoms_names}

                # Otherwise
                else:

                    # If the attribute representing atoms' names is
                    # already in the dictionary of attributes
                    if "label_atom_id" in atoms_attributes:

                        # Warn the user.
                        warnstr = \
                            "A list of accepted values for atoms' " \
                            "names was passed in the dictionary " \
                            "of attributes ('label_atom_id'). It " \
                            "will be extended with the names " \
                            "provided in the 'atoms' list."
                        logger.warning(warnstr)

                        # Extend the accepted values for the attribute.
                        atoms_attributes["label_atom_id"].extend(\
                            atoms_names)

                    # Otherwise
                    else:

                        # Add the attribute representing atoms' names
                        # to the dictionary of attributes.
                        atoms_attributes["label_atom_id"] = atoms_names

        #-------------------------------------------------------------#

        # Return the dictionary of selected items and corresponding
        # attributes - add an empty list if there are no selected
        # items for a given level of the hierarchy, and an empty
        # dictionary if there are no selected attributes.
        return \
            {"models" : \
              {"_items" : \
                 models if models else [],
               "_attributes" : \
                 models_attributes if models_attributes else {}},
             "chains" : \
               {"_items" : \
                  chains if chains else [],
                "_attributes" : \
                  chains_attributes if chains_attributes else {}},
             "segments" : \
               {"_items" : \
                  segments if segments else [],
                "_attributes" : \
                  segments_attributes if segments_attributes else {}},
             "residues" : \
               {"_items" : \
                  residues if residues else [],
                "_attributes" : \
                  residues_attributes if residues_attributes else {}},
             "atoms" : \
               {"_items" : \
                  atoms if atoms else [],
                "_attributes" : \
                  atoms_attributes if atoms_attributes else {}}}


    def _add_bonds(self,
                   bonds,
                   bonds_attributes):
        """Add bonds between atoms (it modifies the connectivity
        data, while the atomic coordinates remain unchanged).

        Parameters
        ----------
        bonds : ``list``
            A list of tuples, each representing a bond between two
            atoms. The tuple must contain the first atom, the
            second atom, the bond's type and the bond's order.

            An atom must be represented by its "path", namely a
            tuple containing the model, chain, segment, and residue
            the atom belongs to, plus the atom's name. This is
            based on the assumption that, within each residue,
            atoms must have unique names.

            For instance, a single covalent bond between
            the N atom of MET 1 in chain A, segment '' and
            the CA atom of the same residue must be represented as
            ``(("A", "", (1, "", "MET"), "N"), 
            ("A", "", (1, "", "MET"), "CA"), "covale", "sing")``.

            Supported bond types are:
            - ``"covale"`` for covalent bonds.
            - ``"disulf"`` for disulfide bridges.
            - ``"hydrog"`` for hydrogen bonds.
            - ``"metalc"`` for metal coordination.

            Supported bond orders are:
            - ``"sing"`` for single bonds.
            - ``"doub"`` for double bonds.
            - ``"trip"`` for triple bonds.
            - ``"quad"`` for quadruple bonds.

            Each bond should be specified only once since,
            for each bond specified, the connectivity data for
            both atoms involved in the bond will be updated.

        bonds_attributes : ``dict``, optional
            A dictionary mapping each bond in ``bonds`` to a dictionary
            of attributes that will be assigned to the bond.
        """

        # Create a copy of the current connectivity data that will be
        # updated with the new bonds.
        conect_data_updated = copy.deepcopy(self.conect_data)

        # If there are no connectivity data
        if not conect_data_updated:

            # Create an empty record for each model found in the
            # structure.
            conect_data_updated = {mod : {} for mod in self.atom_data}

        #-------------------------------------------------------------#

        # For each bond (first atom, second atom, bond type, and bond
        # order)
        for atom_1, atom_2, bond_type, bond_order in bonds:

            # Get the bond's unique identifier.
            bond_id = (atom_1, atom_2, bond_type, bond_order)

            #---------------------------------------------------------#

            # Get the supported bond types.
            supported_bond_types = _defaults.STRUCT_BOND_TYPES

            # If the bond type is invalid
            if bond_type not in supported_bond_types:

                # Get a string representing the supported bond types.
                supported_bond_types_str = \
                    ", ".join([f"'{bt}'" for bt \
                               in supported_bond_types])

                # Raise an error.
                errstr = \
                    f"Unrecognized bond type '{bond_type}' for the " \
                    f"bond between atom {atom_1} and atom {atom_2}. " \
                    "Supported bond types are: " \
                    f"{supported_bond_types_str}."
                raise ValueError(errstr)

            #---------------------------------------------------------#

            # Get the supported bond orders.
            supported_bond_orders = _defaults.STRUCT_BOND_ORDERS

            # If the bond order is invalid
            if bond_order not in supported_bond_orders:

                # Get a string representing the supported bond orders.
                supported_bond_orders_str = \
                    ", ".join([f"'{bo}'" for bo \
                               in supported_bond_orders])

                # Raise an error.
                errstr = \
                    f"Unrecognized bond order '{bond_order}' for " \
                    f"the bond between atom {atom_1} and atom " \
                    f"{atom_2}. Supported bond orders are: " \
                    f"{supported_bond_orders_str}."
                raise ValueError(errstr)

            #---------------------------------------------------------#

            # Create the connectivity record representing the new bond
            # for the first atom.
            new_record_1 = \
                {atom_2 : \
                    {"conn_type_id" : bond_type,
                     "pdbx_value_order" : bond_order}}

            # Create the connectivity record representing the new bond
            # for the second atom.
            new_record_2 = \
                {atom_1 : \
                    {"conn_type_id" : bond_type,
                     "pdbx_value_order" : bond_order}}

            #---------------------------------------------------------#

            # If there are attributes for the bond
            if bonds_attributes is not None \
            and bond_id in bonds_attributes:

                # Add them to the new record for the first atom.
                new_record_1[atom_2].update(bonds_attributes[bond_id])

                # Add them to the new record for the second atom.
                new_record_2[atom_1].update(bonds_attributes[bond_id])

            #---------------------------------------------------------#

            # For each model's connectivity records
            for mod in list(conect_data_updated):

                #-----------------------------------------------------#

                # If we already have a connectivity record for the
                # first atom
                if atom_1 in conect_data_updated[mod]:

                    # Update it with the new record.
                    record_1 = \
                        {**conect_data_updated[mod][atom_1],
                         **new_record_1}

                # Otherwise
                else:

                    # Create a new record for the first atom.
                    record_1 = new_record_1

                # Add the record to the connectivity data.
                conect_data_updated[mod][atom_1] = record_1

                #-----------------------------------------------------#

                # If we already have a connectivity record for the
                # second atom
                if atom_2 in conect_data_updated[mod]:

                    # Update it with the new record.
                    record_2 = \
                        {**conect_data_updated[mod][atom_2],
                         **new_record_2}

                # Otherwise
                else:

                    # Create a new record for the second atom.
                    record_2 = new_record_2
                
                # Add the record to the connectivity data.
                conect_data_updated[mod][atom_2] = record_2

        #-------------------------------------------------------------#

        # Update the connectivity data for the structure.
        self._conect_data = conect_data_updated


    def _get_atom_serial(self,
                         path):
        """Get the serial number of a specific atom, given the 'path'
        to the atom in the structure in terms of the model, chain,
        segment, and residue the atom belongs to, plus the atom's name.

        This method assumes that, within each residue, the atoms'
        names are unique.

        Parameters
        ----------
        path : ``tuple``
            The 'path' to the atom (the model, chain, segment, and
            residue it belongs to, plus the atom's name).

        Returns
        -------
        atom_serial : ``int```
            The serial number of the atom the 'path' points to.
        """
        
        # Get the model, segment, chain, residue, and name of the atom
        # whose serial number should be returned.
        mod, ch, seg, res, atom_name = path

        #-------------------------------------------------------------#

        # Try to get the sub-structure where the atom is stored.
        try:
            
            sub_struct = \
                self.atom_data[mod]["_items"][ch]["_items"][\
                    seg]["_items"][res]["_items"]

        # If the atom was not found
        except KeyError:

            # Return None.
            return

        #-------------------------------------------------------------#

        # For each atom in the sub-structure
        for atom in sub_struct:

            # If the current atom's name matches the one we are looking
            # for.
            if sub_struct[atom]["_attributes"]["label_atom_id"] == \
                atom_name:
                
                # Return the atom's serial number.
                # 
                # There is no need to check any further because we are
                # assuming that, for each residue, the atoms' names
                # are unique).
                return atom


    def _update_atom_and_conect(self,
                                oldids2newids = None,
                                start = 1):
        """Update both atomic coordinates and connectivity data after
        a modification.

        This method does two things:

        - It renumbers the atoms in the structure so that each model's
          atom numbering starts from 'start', and each chain's atom
          numbering starts from the number of the atom ending the
          previous chain plus 2 (to accommodate the TER record at the
          end of a chain, which is needed when we want to write out the
          structure to a PDB file).

        - It modifies the identifiers of the atoms in the connectivity
          data if the models, chains, segments, residues, or atoms were
          renamed.

        Parameters
        ----------
        oldids2newids : ``dict``, optional
            A dictionary mapping the old identifiers for models, chains
            segments, residues, or atoms, to the new ones after
            renaming or renumbering the corresponding items.

            Each identifier is provided as the unique 'path' to the
            item (plus the item's identifier) before the modification
            and after.

        start : ``int``, ``1``
            The new starting number for the atoms' serial numbers.
        """

        # We use iteration instead of recursion because, in this
        # case, it's much simpler and more readable.
        
        # Create a copy of the 'atom_data' dictionary to modify.
        struct_copy = copy.deepcopy(self.atom_data)

        #-------------------------------------------------------------#

        # For each model in the structure
        for mod, mod_d in self.atom_data.items():
            
            # Re-set the current atom count (since, for each model,
            # serial numbers start from 1).
            current_atom_count = 1

            #---------------------------------------------------------#
            
            # For each chain in the current model
            for i, (ch, ch_d) in enumerate(mod_d["_items"].items()):
                
                # If we are not in the first chain of the model
                if i != 0:

                    # The starting number will be the current one plus
                    # 1 to accommodate the previous chain's TER record.
                    current_atom_count = current_atom_count + 1

                #-----------------------------------------------------#
                
                # For each segment in the current chain
                for seg, seg_d in ch_d["_items"].items():

                    #-------------------------------------------------#
                    
                    # For each residue in the current segment
                    for res, res_d in seg_d["_items"].items():

                        # Renumber all atoms in the residue.
                        new_res_d = \
                            {rel_index + current_atom_count : at_d \
                             for rel_index, (_, at_d) \
                             in enumerate(res_d["_items"].items())}

                        # Add the new atoms to the residues,
                        # overwriting the old ones.
                        struct_copy[mod]["_items"][ch][\
                            "_items"][seg]["_items"][res][\
                                "_items"] = new_res_d
                            
                        # Update the current count for the atoms.
                        current_atom_count += len(new_res_d)

        #-------------------------------------------------------------#

        # Update the 'atom_data' dictionary.
        self._atom_data = struct_copy

        #-------------------------------------------------------------#

        # If there is a mapping between the old identifiers and the new
        # identifiers.
        if oldids2newids:

            # Get the first model's identifier.
            first_mod = list(oldids2newids.keys())[0]

            # Try to get the depth of the identifiers.
            try:
                
                depth = \
                    len(list(oldids2newids[first_mod].values())[0])

            # If we are at the models' level
            except AttributeError:

                # Set the depth to 0.
                depth = 0

        #-------------------------------------------------------------#

        # If the structure has associated connectivity data
        if self.conect_data:

            # Initialize an empty dictionary to store the updated
            # connectivity data.
            conect_data_updated = {}

            # For each model in the existent connectivity data
            for mod in self.conect_data:

                #-----------------------------------------------------#

                # For each first atom in the connectivity data for the
                # current model
                for atom_1, atoms_bonded \
                    in self.conect_data[mod].items():

                    # If we have a mapping between old identifiers and
                    # new identifiers
                    if oldids2newids:

                        # If we are substituting the identifiers of
                        # items other than models
                        if depth > 0:

                            # Get the portion of the identifier to
                            # change, if needed, and the one to keep
                            # unchanged.
                            id_to_change, id_to_keep = \
                                atom_1[:depth], atom_1[depth:]

                            # If the portion to change is in the
                            # mapping
                            if id_to_change in oldids2newids[mod]:

                                # Get the new identifier for the atom.
                                atom_1 = \
                                    oldids2newids[mod][\
                                        id_to_change] + id_to_keep

                        # If we are substituting the models'
                        # identifiers
                        else:

                            # If the current model is in the mapping
                            if (mod,) in oldids2newids:

                                # Susbtitute it with the new one.
                                mod = oldids2newids[(mod,)][0]

                    #-------------------------------------------------#

                    # Get the first atom's serial number.
                    atom_1_serial = \
                        self._get_atom_serial((mod, *atom_1))

                    # If the atom is not found (because, for instance,
                    # it was removed from the structure).
                    if atom_1_serial is None:

                        # Continue.
                        continue

                    #-------------------------------------------------#

                    # Initialize an empty dictionary to store the new
                    # record.
                    atoms_bonded_updated = {}

                    #-------------------------------------------------#

                    # For each atom in the record and associated data
                    for atom_2, bond_data in atoms_bonded.items():

                        # If we have a mapping between old identifiers
                        # and new identifiers
                        if oldids2newids:

                            # If we are substituting the identifiers of
                            # items other than models
                            if depth > 0:

                                # Get the portion of the identifier to
                                # change, if needed, and the one to
                                # keep unchanged.
                                id_to_change, id_to_keep = \
                                    atom_2[:depth], atom_2[depth:]

                                # If the portion to change is in the
                                # mapping
                                if id_to_change in oldids2newids[mod]:

                                    # Get the new identifier for the
                                    # atom.
                                    atom_2 = \
                                        oldids2newids[mod][\
                                            id_to_change] + id_to_keep

                        #---------------------------------------------#

                        # Get the second atom's serial number.
                        atom_2_serial = \
                            self._get_atom_serial((mod, *atom_2))

                        # If the atom is not found (because, for
                        # instance, it was removed from the structure).
                        if atom_2_serial is None:

                            # Continue.
                            continue

                        # Add the second atom to the new record.
                        atoms_bonded_updated[atom_2] = bond_data

                    #------------------------------------------------#

                    # If we have an empty record
                    if not atoms_bonded_updated:

                        # Ignore it.
                        continue

                    # If the model is already in the dictionary
                    if mod in conect_data_updated:

                        # Add the record for the current first atom.
                        conect_data_updated[mod][atom_1] = \
                            atoms_bonded_updated

                    # Otherwise
                    else:

                        # Create an entry for the model and add a
                        # record for the current atom.
                        conect_data_updated[mod] = \
                            {atom_1 : atoms_bonded_updated}

            #---------------------------------------------------------#

            # Update the connectivity data.
            self._conect_data = conect_data_updated


    def _add_atom(self,
                  atom_path,
                  atom_attributes):
        """Add a single atom to the structure, given the 'path' to the
        point in the structure where the atom should be added.

        Parameters
        ----------
        atom_path : ``tuple``
            The 'path' to the point in the structure where the atom
            should be added.

        atom_attributes : ``dict``
            A dictionary of attributes for the atom.
        """

        # Define the recursion to add the atom at the correct point
        # in the structure.
        def recurse(struct,
                    atom_attributes,
                    current_depth = 1,
                    current_path = ()):

            # If we are at the residues' depth
            if current_depth == 8:

                # Get the updated atoms for the residue, adding the one
                # of interest and the corresponding attributes.
                new_items = \
                    {**struct["_items"],
                     **{current_path[0] : \
                        {"_attributes" : atom_attributes}}}

                # Update the residue's atoms.
                struct["_items"] = new_items

                # Return the residue.
                return struct

            #---------------------------------------------------------#

            # If we are at the models' depth
            elif current_depth == 1:

                # For each key and associated sub-structure in the
                # current structure
                for key, sub_struct in struct.items():

                    # Recurse through the sub-structure without
                    # updating the current path.
                    struct[key] = \
                        recurse(struct = sub_struct,
                                atom_attributes = atom_attributes,
                                current_depth = current_depth + 1,
                                current_path = current_path)

            #---------------------------------------------------------#

            # Otherwise
            else:

                # For each key and associated sub-structure in the
                # current structure
                for key, sub_struct in struct.items():

                    # If we are at an '_items' level.
                    if key == "_items":

                        # Recurse through the structure without
                        # updating the current path.
                        struct[key] = \
                            recurse(\
                                struct = sub_struct,
                                atom_attributes = atom_attributes,
                                current_depth = current_depth + 1,
                                current_path = current_path)

                    # Otherwise
                    elif key == current_path[0]:

                        # Recurse through the sub-structure.
                        struct[key] = \
                            recurse(\
                                struct = sub_struct,
                                atom_attributes = atom_attributes,
                                current_depth = current_depth + 1,
                                current_path = current_path[1:])

            #---------------------------------------------------------#

            # Return the structure.
            return struct

        #-------------------------------------------------------------#

        # For each model in the structure
        for mod in self.atom_data:

            # Try to access the path to the atom.
            try:

                self[(mod, *atom_path)]

            # If the path does not exist
            except KeyError:

                # Raise an error.
                errstr = \
                    f"Chain '{atom_path[0]}', " \
                    f"segment '{atom_path[1]}', " \
                    f"residue '{atom_path[2]}' does not " \
                    f"exist in model {mod}."

        #-------------------------------------------------------------#

        # If the user did not specify the atom's name
        if "label_atom_id" not in atom_attributes:

            # Raise an error.
            errstr = \
                "The atom's name must be specified ('label_atom_id')."
            raise KeyError(errstr)

        # If the user did not specify the atom's alternative location
        if "label_alt_id" not in atom_attributes:

            # Set it to an empty string.
            atom_attributes["label_alt_id"] = ""

        # If the user did not specify the atom's entity ID
        if "label_entity_id" not in atom_attributes:

            # Set it to an empty string.
            atom_attributes["label_entity_id"] = ""

        # If the user did not specify the atom's type
        if "type_symbol" not in atom_attributes:

            # Raise an error.
            errstr = \
                "The atom's type must be specified ('type_symbol')."
            raise KeyError(errstr)

        # If the user did not specify the X coordinate of the atom's
        # position
        if "Cartn_x" not in atom_attributes:

            # Raise an error.
            errstr = \
                "The X coordinate of the atom's position must be " \
                "specified ('Cartn_x')."
            raise KeyError(errstr)

        # If the user did not specify the Y coordinate of the atom's
        # position
        if "Cartn_y" not in atom_attributes:

            # Raise an error.
            errstr = \
                "The Y coordinate of the atom's position must be " \
                "specified ('Cartn_y')."
            raise KeyError(errstr)

        # If the user did not specify the Z coordinate of the atom's
        # position
        if "Cartn_z" not in atom_attributes:

            # Raise an error.
            errstr = \
                "The Z coordinate of the atom's position must be " \
                "specified ('Cartn_z')."
            raise KeyError(errstr)

        # If the user did not specify the atom's formal charge
        if "pdbx_formal_charge" not in atom_attributes:

            # Set it to an empty string.
            atom_attributes["pdbx_formal_charge"] = ""

        #-------------------------------------------------------------#

        # Get the serial number of the last atom.
        last_atom = self.get_items(action = "get_max",
                                   level = "atoms",
                                   squeeze = "everything")

        #-------------------------------------------------------------#

        # Recursively modify the structure.
        self._atom_data = \
            recurse(struct = copy.deepcopy(self.atom_data),
                    atom_attributes = atom_attributes,
                    current_path = atom_path + (last_atom+1,))

        #-------------------------------------------------------------#

        # Renumber the atoms and update the connectivity data.
        self._update_atom_and_conect()


    def _get(self,
             action,
             level,
             squeeze,
             elements_type,
             selected,
             attribute = None):
        """Generic 'get' method to retrieve elements from the
        'atom_data' dictionary.

        Parameters
        ----------
        action : ``str``, {``"get_unique"``, ``"get_max"``, \
            ``"get_min"``}
            The name of the 'get' action to perform:

            - ``"get_unique"`` retrieves all unique elements of a
               specific type or all unique values for a specific
               attribute.

            - ``"get_min"`` retrieves the minimum value found for all
              elements of a specific type or the minimum value found
              for a specific attribute.

            - ``"get_max"`` retrieves the maximum value found for all
              elements of a specific type or the maximum value found
              for a specific attribute.

        level : ``str``, {``"models"``, ``"chains"``, \
            ``"segments"``, ``"residues"``, ``"atoms"``}
            At what level of the hierarchy we are getting the elements.

        squeeze : ``str``, {``"models"``, ``"chains"``, \
            ``"segments"``, ``"residues"``, ``"atoms"``}
            At which level of the hierarchy items should be 'squeezed',
            meaning that the 'get' action will involve all values found
             at deeper levels and consider them all together.

            For instance, if ``action = "get_unique"``, 
            ``level = "residues"``, and ``squeeze = "chains"``, the
            result will be a dictionary mapping each model to the
            chains it contains, and each chain will be mapped to a list
            containing all unique residues found in the chain. This
            will ignore the segments they residues belong to, meaning
            that if residues with the same ID are found in different
            segments, they will be reported only once.

            If ``squeeze`` is not specified, it will automatically be
            set to the level above ``level``.

            The hierarchy goes: models > chains > segments > residues >
            atoms.

        elements_type : ``str``, {``"items"``, ``"attributes"``}
            The type of elements we are retrieving.

        selected : ``dict``
            The dictionary containing items (selected by either their
            IDs or by their attributes) that will be considered in
            the 'get' action.

        attribute : ``str``, optional
            The name of the attribute whose values should be retrieved.

        Returns
        -------
        result : ``dict``
            A dictionary containing the result of the 'get' action
            performed.
        """

        # Define the inner recursion to get the desired values below
        # 'merge_depth'
        def _inner_recurse(struct,
                           selected,
                           elements_type,
                           current_depth):

            # If we reached the target depth
            if current_depth == target_depth:

                # If the target depth is the models' depth
                if current_depth == 1:
                    
                    # The depth of the items we are checking is only
                    # one level down in the '_ACTION2DEPTH2ITEMS'
                    # dictionary.
                    items_depth = current_depth + 1
                
                # Otherwise
                else:
                    
                    # The depth of the items we are checking is the
                    # level down in the '_ACTION2DEPTH2items'
                    # dictionary.
                    items_depth = current_depth + 2

                #-----------------------------------------------------#

                # Get the user-selected items.
                selected_items = \
                    selected[self._ITEMSATTRS_DEPTH_TO_LEVELS[\
                        items_depth]]["_items"]

                # Get the user-selected attributes.
                selected_attrs = \
                    selected[self._ITEMSATTRS_DEPTH_TO_LEVELS[\
                        items_depth]]["_attributes"]

                #-----------------------------------------------------#

                # If we are checking models
                if current_depth == 1:

                    # If no selected models were provided
                    if not selected_items:

                        # All models will be included in the selection.
                        selected_items = list(struct.keys())

                    # The items to check will be the entire current
                    # structure.
                    items_to_check = struct

                    # The attributes of the previous items in the
                    # hierarchy will be an empty dictionary (models
                    # are the first items in the hierarchy).
                    attrs_prev = {}

                    # There are no user-selected attributes for the
                    # previous items in the hierarchy (models are the
                    # first items in the hierarchy).
                    selected_attrs_prev = {}

                #-----------------------------------------------------#

                # Otherwise
                else:

                    # If no selected items were provided
                    if not selected_items:

                        # All items will be included.
                        selected_items = list(struct["_items"].keys())
                    
                    # The items to check will be the ones in the
                    # '_items' dictionary.
                    items_to_check = struct["_items"]

                    # The attributes of the previous items in the
                    # hierarchy will be the attributes associated with
                    # the current structure.
                    attrs_prev = struct["_attributes"]

                    # Get the user-selected attributes for the previous
                    # items in the hierarchy (they sit exactly at the
                    # current depth).
                    selected_attrs_prev = \
                        selected[self._ITEMSATTRS_DEPTH_TO_LEVELS[\
                            current_depth]]["_attributes"]

                #-----------------------------------------------------#

                # Initialize an empty dictionary where the items that
                # will be considered will be stored.
                considered_items = {}

                #-----------------------------------------------------#
                
                # For each item and associated dictionary in the
                # current structure.
                for item, item_dict in items_to_check.items():

                    #-------------------------------------------------#

                    # Set a flag for whether to consider the item and
                    # initialize it to False.
                    to_consider = False

                    #-------------------------------------------------#
                    
                    # If the item is a tuple (it is a residue)
                    if isinstance(item, tuple):

                        # If the first element is in the selected
                        # items
                        if item[0] in selected_items:

                            # Consider the current item.
                            to_consider = True

                        # If the third element is in the selected
                        # items
                        if item[2] in selected_items:

                            # Consider the current item.
                            to_consider = True

                    #-------------------------------------------------#
                    
                    # If the item is in the selected ones
                    if item in selected_items:

                        # Do not consider the current item.
                        to_consider = True

                    #-------------------------------------------------#

                    # For each attribute of the previous item in the
                    # hierarchy
                    for attr_name, attr_value in attrs_prev.items():

                        # If the attribute is one of those we have to
                        # check, and its value is not in the list of
                        # accepted values.
                        if attr_name in selected_attrs_prev \
                        and attr_value \
                            not in selected_attrs_prev[attr_name]:

                            # Do not consider the current item.
                            to_consider = False

                    #-------------------------------------------------#

                    # If we are collecting items (rather than
                    # attributes)
                    if elements_type == "items":

                        # For each item's attribute and the
                        # attribute's associated value
                        for attr_name, attr_value \
                        in item_dict["_attributes"].items():
 
                            # If the attribute is one of those we have
                            # to check, and its value is not in the
                            # list of accepted values
                            if attr_name in selected_attrs \
                            and attr_value \
                                not in selected_attrs[attr_name]:
                                
                                # Do not consider the item.
                                to_consider = False

                    #-------------------------------------------------#
                    
                    # If wee need to consider the item
                    if to_consider:

                        # Add the item and associted dictionary to the
                        # items to be considered.
                        considered_items[item] = item_dict

                #-----------------------------------------------------#

                # If we are collecting items
                if elements_type == "items":

                    # Return the unique items found at this level (it
                    # should be all of them since dictionary keys are
                    # unique).
                    return set(considered_items.keys())

                # If we are collecting attributes
                elif elements_type == "attributes":

                    # Return the unique values for the attribute of
                    # interest in the considered items.
                    return \
                        set([considered_items[k]["_attributes"][\
                                selected_attrs[0]] \
                             for k in considered_items])

            #---------------------------------------------------------#

            # Otherwise
            else:

                # Create an empty set to store the collected values.
                collected_values = set()

                # For each key and associated sub-structure in the
                # current structure
                for key, sub_struct in struct.items():

                    # Get whether we should recurse.
                    should_recurse = \
                        self._should_recurse(\
                            key = key,
                            value = sub_struct,
                            selected = selected,
                            current_depth = current_depth)

                    # If we should recurse
                    if should_recurse:

                        # Update the collected values.
                        collected_values.update(\
                            _inner_recurse(\
                                struct = sub_struct,
                                selected = selected,
                                elements_type = elements_type,
                                current_depth = current_depth + 1))

                # Return the collected values.
                return collected_values

        #-------------------------------------------------------------#

        # Define the outer recursion to get the desired values above
        # and at 'merge_depth'.
        def _outer_recurse(struct,
                           action,
                           selected,
                           target_depth,
                           merge_depth,
                           elements_type,
                           current_depth = 1):

            # If the current depth is the depth at which values should
            # be merged
            if current_depth == merge_depth:

                # Initialize an empty dictionary to store the current
                # result.
                result = {}
                
                # For each key and associated sub-structure in the
                # current structure
                for key, sub_struct in struct.items():

                    # Get whether we should recurse.
                    should_recurse = \
                        self._should_recurse(\
                            key = key,
                            value = sub_struct,
                            selected = selected,
                            current_depth = current_depth)

                    # If we should recurse
                    if should_recurse:

                        # Get the collected values.
                        collected_values = \
                            _inner_recurse(\
                                struct = sub_struct,
                                selected = selected,
                                elements_type = elements_type,
                                current_depth = current_depth + 1)

                        # If we need to get the unique values
                        if action == "get_unique":

                            # Convert the set into a list, sort it,
                            # and return it.
                            result[key] = sorted(collected_values)

                        # If we need to get the minimum value
                        elif action == "get_min":

                            # Find the minimum value and return it.
                            result[key] = min(collected_values)

                        # If we need to get the maximum
                        elif action == "get_max":

                            # Find the maximum value and return it.
                            result[key] = max(collected_values)

                # Return the result.
                return result

            #---------------------------------------------------------#

            # Otherwise
            else:

                # Create an empty dictionary to store the result.
                result = {}

                # For each key and associated sub-structure in the
                # current structure.
                for key, sub_struct in struct.items():

                    # Get whether we should recurse.
                    should_recurse = \
                        self._should_recurse(\
                            key = key,
                            value = sub_struct,
                            selected = selected,
                            current_depth = current_depth)

                    # If we should recurse
                    if should_recurse:

                        # If we are at a ['_items', '_attributes']
                        # level
                        if key in ("_items", "_attributes"):

                            # If we are at the '_attributes'
                            if key == "_attributes":

                                # Ignore them.
                                continue
                            
                            # The sub-structure to recurse through
                            # will be the current sub-structure.
                            result[key] = \
                                _outer_recurse(\
                                    struct = sub_struct,
                                    action = action,
                                    selected = selected,
                                    target_depth = target_depth,
                                    merge_depth = merge_depth,
                                    elements_type = elements_type,
                                    current_depth = current_depth + 1)

                        #---------------------------------------------#

                        # Otherwise
                        else:

                            # The sub-structure to recurse through
                            # will be the '_items' dictionary in
                            # the current sub-structure (in this
                            # case, 'current_depth' becomes
                            # 'current_depth + 2' for the next
                            # level since we are skipping the
                            # ['_items', '_attribute'] level just
                            # below.
                            result[key] = \
                                _outer_recurse(\
                                    struct = sub_struct["_items"],
                                    action = action,
                                    selected = selected,
                                    target_depth = target_depth,
                                    merge_depth = merge_depth,
                                    elements_type = elements_type,
                                    current_depth = current_depth + 2)

                #-----------------------------------------------------#

                # Return the result.
                return result

        #-------------------------------------------------------------#

        # Get the target depth from the 'level' keword.
        target_depth = \
            {v : k for k, v \
             in self._ITEMSATTRS_DEPTH_TO_LEVELS.items()}[level] - 2

        #-------------------------------------------------------------#

        # If 'squeeze' was not specified
        if squeeze is None:

            # The 'merge_depth' will be just one level above (for
            # instance, 'segments' for 'residues' and 'chains' for
            # 'segments').
            merge_depth = target_depth - 1

        # If we need to 'squeeze' everything into a list/value in the
        # end
        elif squeeze == "everything":

            # The merge depth will be that of models.
            merge_depth = \
                {v : k for k, v \
                 in self._IDS_DEPTH_TO_LEVELS.items()}["models"]

        # Otherwise
        else:

            # The 'merge_depth' will be the specified level.
            merge_depth = \
                {v : k for k, v \
                 in self._IDS_DEPTH_TO_LEVELS.items()}[squeeze]

        #-------------------------------------------------------------#

        # If an attribute was specified
        if attribute is not None:

            # If attributes for the given level were specified
            if selected[level]["_attributes"]:

                # Warn the user.
                warnstr = \
                    "When an 'attribute' is specified for level " \
                    f"'{level}', the '{level}_attributes' option " \
                    "is ignored."
                logger.warning(warnstr)

            # Set the attribute for the given level.
            selected[level]["_attributes"] = [attribute]

        #-------------------------------------------------------------#
        
        # Recursively traverse the structure and get the result.
        result = \
            _outer_recurse(\
                struct = self.atom_data,
                action = action,
                selected = selected,
                target_depth = target_depth,
                merge_depth = merge_depth,
                elements_type = elements_type)

        #-------------------------------------------------------------#

        # If we need to 'squeeze' everything into one list/value
        if squeeze == "everything":

            # Get all values from the dictionary of results.
            values = list(result.values())

            # We are sure that we either have one value or a
            # dictionary of depth 1 because of how we recursed
            # through the structure if 'squeeze' was 'everything'.

            # If we have only one value (= we had only one model)
            if len(values) == 1:

                # If the values are a list of NOT lists
                if not isinstance(values[0], list):

                    # The 'unique' values will be just the values
                    # we have.
                    unique_values = values

                # If the values are a list of lists
                else:

                    # The 'unique' values will be the ones in the
                    # first position of the list.
                    unique_values = values[0]

            # If we are multiple values (= we had multiple models)
            else:

                # If the values are a list of NOT lists
                if not isinstance(values[0], list):

                    # Get the unique values.
                    unique_values = set(values)

                # Otherwise.
                else:

                    # Get the unique values from the list of lists.
                    unique_values = \
                        set(itertools.chain.from_iterable(values))

            # If we need to get unique values
            if action == "get_unique":
                
                # Convert the set of unique values into a list,
                # sort it, and return it.
                return sorted(list(unique_values))

            # If we need to get the minimum value
            elif action == "get_min":

                # Return the minimum value found in the set.
                return min(unique_values)

            # If we need to get the maximum value
            elif action == "get_max":

                # Return the maximum value found in the set.
                return max(unique_values)

        #-------------------------------------------------------------#

        # Otherwise
        else:

            # Return the result as it is
            return result


    def _merge(self,
               other):
        """Merge two structures.

        Parameters
        ----------
        other : :class:`pdbcraft.structure.Structure`
            The structure to be merged with the current one.

        Returns
        -------
        merged : :class:`pdbcraft.structure.Structure`
            The merged structure.
        """

        # Define the recursion to merge the two structures.
        def recurse(struct,
                    other_struct):

            # If none of the current sub-structures is a dictionary
            if not isinstance(struct, dict) \
            or not isinstance(other_struct, dict):

                # Return the first one as it is.
                return struct

            # Create a copy of the first sub-structure, representing
            # the initial 'merged' sub-structure.
            merged = dict(struct)

            # For each key and associated sub-structure in the
            # 'other' current sub-structure
            for key, sub_struct in other_struct.items():

                # If the key is already in the merged structure
                # (= they key is already in the 'self' structure)
                if key in merged:

                    # Add the values found associated with the same
                    # key in the 'other' sub-structure.
                    merged[key] = \
                        recurse(struct = merged[key],
                                other_struct = sub_struct)
                
                # Otherwise (= the key is only in the 'other'
                # sub-structure)
                else:

                    # Add the entire sub-structure the key is
                    # associated with to the merged structure.
                    merged[key] = sub_struct

            # Return the merged structure.
            return merged

        #-------------------------------------------------------------#

        # If 'other' is not a 'Structure'
        if not isinstance(other, self.__class__):

            # Raise an error.
            errstr = \
                f"'other' must be a {self.__class__.__name__} " \
                "instance."
            raise TypeError(errstr)

        # If the two structures do not have the same number of models
        if len(self.atom_data) != len(other.atom_data):

            # Raise an error.
            errstr = \
                f"The current structure has {len(self.atom_data)} " \
                f"model, while 'other' has {len(other.atom_data)} " \
                "models. To be merged, two structures must have " \
                "the same number of models."
            raise ValueError(errstr)

        #-------------------------------------------------------------#

        # Get the serial number of the last atom.
        last_atom = self.get_items(action = "get_max",
                                   level = "atoms",
                                   squeeze = "everything")

        #-------------------------------------------------------------#

        # Renumber the atoms in the 'other' structure starting from the
        # last atom in the current structure - this is needed so that
        # no atoms get overwritten while merging the two structures.
        other._update_atom_and_conect(start = last_atom + 1)

        #-------------------------------------------------------------#

        # Merge the atomic coordinates of the two structures.
        merged_atom_data = \
            recurse(struct = self.atom_data,
                    other_struct = other.atom_data)

        # Merge the connectivity data of the two structures.
        merged_conect_data = \
            recurse(struct = self.conect_data,
                    other_struct = other.conect_data)

        # Create the merged structure
        merged = \
            self.__class__(atom_data = merged_atom_data,
                           conect_data = merged_conect_data)

        # Renumber the atoms of the merged structure starting from 1
        # and update the connectivity data.
        merged._update_atom_and_conect()

        #-------------------------------------------------------------#

        # Return the merged structure.
        return merged


    def _keep(self,
              struct,
              selected_items,
              selected_attrs,
              current_depth):
        """Keep only selected items in the structure or in a portion
        of the structure.

        It is called by the '_modify' function when needed.

        Parameters
        ----------
        struct : ``dict``
            Either the full 'atom_data' dictionary or a portion of it.

        selected_items : ``list``
            A list of items to be kept.

        selected_attrs : ``dict``
            A dictionary mapping each attribute of interest to a list
            of values for the attribute.

            Only items for which the attribute's value is among the
            ones in the list will be kept.

        current_depth : ``int``
            An integer representing the depth we are currently at in
            the 'atom_data' dictionary.

        Returns
        -------
        kept_items : ``dict``
            A dictionary containing the items to be kept and their
            attributes.
        """

        # If we are keeping models
        if current_depth == 1:

            # If no selected models were provided
            if not selected_items:

                # Add models will be included.
                selected_items = list(struct.keys())

            # The items to check will be the entire structure.
            items_to_check = struct

        # Otherwise
        else:

            # If no selected items were provided
            if not selected_items:

                # All items will be included.
                selected_items = list(struct["_items"].keys())
            
            # The items to check will be the ones in the '_items'
            # dictioanary.
            items_to_check = struct["_items"]

        #-------------------------------------------------------------#

        # Initialize an empty dictionary where the items that will
        # be kept will be stored.
        kept_items = {}

        #-------------------------------------------------------------#
        
        # For each item and associated dictionary in the structure.
        for item, item_dict in items_to_check.items():

            #---------------------------------------------------------#

            # Initialize a flag for whether to keep the item to
            # False (= by default, we do not keep the item).
            to_keep = False

            #---------------------------------------------------------#
            
            # If the item is a tuple AND we are at the residues'
            # level (we use both conditions in case we will
            # introduce tuple identifiers for other levels in
            # the future)
            if isinstance(item, tuple) and current_depth == 6:

                # If the first element is in the selected items
                if item[0] in selected_items:

                    # We will keep the item.
                    to_keep = True

                # If the third element is in the selected items
                if item[2] in selected_items:

                    # We will keep the item.
                    to_keep = True

            #---------------------------------------------------------#
            
            # If the item is in the selected ones
            if item in selected_items:

                # We will keep the item.
                to_keep = True

            #---------------------------------------------------------#
            
            # For each item's attribute and the attribute's
            # associated value
            for attr_name, attr_value in \
                item_dict["_attributes"].items():

                # If the attribute is one of those selected,
                # and its value is not among the accepted values
                if attr_name in selected_attrs \
                and attr_value not in selected_attrs[attr_name]:
                    
                    # We will not keep the item.
                    to_keep = False

            #---------------------------------------------------------#
            
            # If wee need to keep the item
            if to_keep:

                # Add the item and associted dictionary to the
                # items to be kept.
                kept_items[item] = item_dict

        #-------------------------------------------------------------#

        # If we kept models
        if current_depth == 1:

            # Just return the dictionary.
            return kept_items

        # Otherwise
        else:
            
            # Return the items to be kept and the current
            # structure's attributes.
            return {"_items" : kept_items,
                    "_attributes" : \
                        {k : v for k, v \
                         in struct["_attributes"].items()}}

    def _remove(self,
                struct,
                selected_items,
                selected_attrs,
                current_depth):
        """Remove selected items from the full structure or from
        a portion of the structure.

        It is called by the '_modify' function when needed.

        Parameters
        ----------
        struct : ``dict``
            Either the full 'atom_data' dictionary or a portion of it.

        selected_items : ``list``
            A list of items to be removed.

        selected_attrs : ``dict``
            A dictionary mapping each attribute of interest to a list
            of values for the attribute.

            Only items for which the attribute's value is not among the
            ones in the list will be kept.

        current_depth : ``int``
            An integer representing the depth we are currently at in
            the 'atom_data' dictionary.

        Returns
        -------
        kept_items : ``dict``
            A dictionary containing the items to be kept and their
            attributes.
        """

        # If we are removing models
        if current_depth == 1:

            # The items to check for removal will be the entire
            # structure.
            items_to_check = struct

        # Otherwise
        else:
            
            # The items to check for removal will be the ones in the
            # '_items' dictioanary.
            items_to_check = struct["_items"]

        # If no selected items were provided
        if not selected_items:

            # No items will be selected for removal.
            selected_items = {}

        #-------------------------------------------------------------#

        # Initialize an empty dictionary where the items that will be
        # kept will be stored.
        kept_items = {}

        #-------------------------------------------------------------#
        
        # For each item and associated dictionary in the structure
        for item, item_dict in items_to_check.items():

            # Initialize a flag for whether to keep the item to True
            # (= by default, we keep the item).
            to_keep = True

            #---------------------------------------------------------#

            # If the item is a tuple AND we are at the residues' level
            # (we use both conditions in case we will introduce tuple
            # identifiers for other levels in the future).
            if isinstance(item, tuple) and current_depth == 6:

                # If the first element is in the selected items
                if item[0] in selected_items:

                    # We will remove the item.
                    to_keep = False

                # If the third element is in the selected items
                if item[2] in selected_items:

                    # We will remove the item.
                    to_keep = False

            #---------------------------------------------------------#
            
            # If the item is in the selected ones
            if item in selected_items:

                # We will remove the item.
                to_keep = False

            #---------------------------------------------------------#
            
            # For each item's attribute and the attribute's associated
            # value
            for attr_name, attr_value \
            in item_dict["_attributes"].items():

                # If the attribute is one of those selected, and its
                # value is among the values in the list
                if attr_name in selected_attrs \
                and attr_value in selected_attrs[attr_name]:

                    # We will remove the item.
                    to_keep = False

            #---------------------------------------------------------#
            
            # If wee need to keep the item
            if to_keep:

                # Add the item and associated dictionary to the items
                # to be kept.
                kept_items[item] = item_dict

        #-------------------------------------------------------------#

        # If we removed models
        if current_depth == 1:

            # Just return the kept items.
            return kept_items

        # Otherwise
        else:

            # Return the items to be kept and the current structure's
            # attributes.
            return {"_items" : kept_items,
                    "_attributes" : \
                        {k : v for k, v \
                         in struct["_attributes"].items()}}


    def _rename_or_assign(self,
                          struct,
                          selected_items,
                          selected_attrs,
                          mapping,
                          elements_type,
                          attribute,
                          new_value,
                          current_depth,
                          current_path,
                          oldids2newids):
        """Rename selected items or assign a new value to one of their
        attributes for items in the the full structure or in a portion
        of the structure.

        It is called by the '_modify' function when needed.

        Parameters
        ----------
        struct : ``dict``
            Either the full 'atom_data' dictionary or a portion of it.

        selected_items : ``list``
            A list of items to be renamed or for which an attribute's
            value will be reassigned.

        selected_attrs : ``dict``
            A dictionary mapping each attribute of interest to a list
            of values for the attribute.

            Only items for which the attribute's value is among the
            ones in the list will be renamed.

        mapping : ``dict``
            A dictionary mapping the old names to the new names.

        elements_type : ``str``, {``"items"``, ``"attributes"``}
            Whether the renaming needs to be performed on items or
            attributes.

        attribute : ``str``, optional
            The name of the attribute to be renamed or to which
            a new value shoud be assigned, if we are doing either
            thing.

        new_value : ``str``, ``int``, ``float``, ``tuple``
            The new value for the attribute, if we are assigning a new
            value to an attribute.    

        current_depth : ``int``
            An integer representing the depth we are currently at in
            the 'atom_data' dictionary.

        current_path : ``tuple``
            The 'path' to the point in the 'atom_data' dictionary we
            are currently at (in terms of keys in the 'atom_data'
            dictionary that we need to traverse in order to reach the
            current point).

        oldids2newids : ``dict``
            A dictionary mapping the old identifiers to the new
            identifiers, in case we are renaming items.

        Returns
        -------
        updated_items : ``dict``
            A dictionary containing the updated items and their
            attributes.
        """

        # If we are renaming models
        if current_depth == 1:

            # If no selected models were provided
            if not selected_items:

                # Add models will be included in the renaming.
                selected_items = list(struct.keys())

            # The items to check will be the entire structure.
            items_to_check = struct

        # If we are renaming chains, segments, residues, or atoms
        else:

            # If no selected items were provided
            if not selected_items:

                # All items in the current structure will be included
                # in the renaming.
                selected_items = list(struct["_items"].keys())
            
            # The items to check will be the ones in the '_items'
            # dictioanary.
            items_to_check = struct["_items"]

        #-------------------------------------------------------------#

        # Initialize an empty dictionary where the updated items
        # will be stored - we do this instead of initializing a
        # dictionary with all the items as they are and then
        # renaming them on the go because by adding them (renamed or
        # in their original form) one by one we ensure that we keep
        # the original order.
        updated_items = {}

        #-------------------------------------------------------------#
        
        # For each item and associated dictionary in the current
        # structure
        for item, item_dict in items_to_check.items():

            # Set a flag for whether to rename the item and initialize
            # it to False.
            to_rename = False

            #---------------------------------------------------------#
            
            # If the item is a tuple (it is a residue)
            if isinstance(item, tuple):

                # If the first element is in the selected items
                if item[0] in selected_items:

                    # Set the flag to True.
                    to_rename = True

                # If the third element is in the selected items
                if item[2] in selected_items:

                    # Set the flag to True.
                    to_rename = True

            #---------------------------------------------------------#
            
            # If the item is in the selected ones
            if item in selected_items:

                # Set the flag to True.
                to_rename = True

            #---------------------------------------------------------#
            
            # For each item's attribute and the attribute's associated
            # value
            for attr_name, attr_value \
            in item_dict["_attributes"].items():

                # If the attribute is one of those selected, and its
                # value is among the accepted values
                if attr_name in selected_attrs \
                and attr_value in selected_attrs[attr_name]:
                    
                    # Set the flag to True.
                    to_rename = True

            #---------------------------------------------------------#

            # If the 'to_rename' flag is False
            if not to_rename:

                # Add the item as it is.
                updated_items[item] = item_dict

                # Go to the next item.
                continue

            #---------------------------------------------------------#

            # If we are renaming an attribute
            if elements_type == "attributes":

                # Initialize the dictionary of updated attributes
                # as the current dictionary of attributes
                # (= not renamed).
                updated_attrs = item_dict["_attributes"]

                # If we have the attribute in the current item
                if attribute in updated_attrs:

                    # If we are renaming an atom
                    if attribute == "label_atom_id" \
                    and current_depth == 8:
                        
                        # Get the old name of the atom.
                        old_value = updated_attrs[attribute]
                        
                        # Update the dictionary mapping the old
                        # identifiers to the new identifiers with
                        # the old and new name of the atomn (we are
                        # sure the path is at least two items long,
                        # because we are at the depth where the atoms
                        # are if we encounter this attribute).
                        oldids2newids[current_path[0]][\
                            current_path[1:] + (old_value,)] = \
                                current_path[1:] + (new_value,)

                    # Update the attribute's value.
                    updated_attrs[attribute] = new_value

                # Update the attribute in the item's dictionary.
                item_dict["_attributes"] = updated_attrs

                # Add the item to the list of updated items.
                updated_items[item] = item_dict
                
                # Go to the next item.
                continue

            #---------------------------------------------------------#
            
            # If the item is a tuple AND we are at the residues'
            # level (we use both conditions in case we introduce
            # tuple identifiers for other levels in he future)
            if isinstance(item, tuple) and current_depth == 6:

                # If the first element of the identifier (= the
                # residue's sequence number) is in the mapping
                # provided
                if item[0] in mapping:

                    # Assemble the new identifier by substituting the
                    # residue's original sequence number with the
                    # corresponding sequence number provided in the
                    # mapping, and leaving the residue's insertion
                    # code and name unchanged.
                    new_item = (mapping[item[0]], item[1], item[2])

                    # Add the new identifier to the dictionary of
                    # updated items.
                    updated_items[new_item] = item_dict

                    # Update the mapping between the old identifiers
                    # and the new identifiers (we are sure the path
                    # is at least two items long because we are at
                    # the residues' level).
                    oldids2newids[current_path[0]][\
                        current_path[1:] + (item,)] = \
                            current_path[1:] + (new_item,)

                    # Go to the next item (we need to continue
                    # because, if we do not do it and the residue's
                    # name or the entire ID of the residue matches a
                    # key in the mapping provided, we are going to add
                    # the residue multiple times).
                    continue

                # If the third element of the identifier (= the
                # residue's name) is in the mapping provided
                if item[2] in mapping:

                    # Assemble the new identifier by substituting the
                    # residue's original name with the corresponging
                    # name found in the mapping, and leaving the
                    # residue's sequence number and insertion code
                    # unchanged.
                    new_item = (item[0], item[1], mapping[item[2]])

                    # Add the new identifier to the dictionary of
                    # updated items.
                    updated_items[new_item] = item_dict

                    # Update the mapping between the old identifiers
                    # and the new identifiers (we are sure the path
                    # is at least two items long because we are at
                    # the residues' level).
                    oldids2newids[current_path[0]][\
                        current_path[1:] + (item,)] = \
                            current_path[1:] + (new_item,)

                    # Go to the next item (we need to continue
                    # because, if we do not do it and the residue's
                    # sequence number or the entire ID of the residue
                    # matches a key in the mapping provided, we are
                    # going to add the same residue multiple times).
                    continue

            #---------------------------------------------------------#

            # If the entire item is in the mapping
            if item in mapping:
                
                # Add the new identifier to the dictionary of updated
                # items.
                updated_items[mapping[item]] = item_dict

                # If the path is longer than one item
                if len(current_path) > 1:

                    # Update the mapping between the old identifiers
                    # and the new identifiers, assuming the path is
                    # at least two items long.
                    oldids2newids[current_path[0]][\
                        current_path[1:] + (item,)] = \
                            current_path[1:] + (mapping[item],)

                # If the path is 1 item long
                elif len(current_path) == 1:

                    # Update the mapping between the old identifiers
                    # and the new identifiers, assuming the path is
                    # one item long.
                    oldids2newids[current_path[0]][(item,)] = \
                        (mapping[item],)

                # If the path is 0 items long (= we are at the models'
                # level)
                else:

                    # Update the mapping between the old models'
                    # identifiers and the new models' idenfitiers.
                    oldids2newids[(item,)] = (mapping[item],)

                # Go to the next item
                continue

            #---------------------------------------------------------#

            # Add the item as it is to the dictionary of updated items.
            updated_items[item] = item_dict

        #-------------------------------------------------------------#

        # If we updated models
        if current_depth == 1:

            # Just return the updated items.
            return updated_items
        
        # Otherwise
        else:
            
            # Return the updated items and the current structure's
            # attributes.
            return {"_items" : updated_items,
                    "_attributes" : \
                        {k : v for k, v \
                         in struct["_attributes"].items()}}


    def _renumber(self,
                  struct,
                  start,
                  current_depth,
                  current_path,
                  oldids2newids):
        """Renumber selected items.

        It is called by the '_modify' function when needed.

        Parameters
        ----------
        struct : ``dict``
            Either the full 'atom_data' dictionary or a portion of it.

        start : ``imt``
            The new starting point for the renumbering. 

        current_depth : ``int``
            An integer representing the depth we are currently at in
            the 'atom_data' dictionary.

        current_path : ``tuple``
            The 'path' to the point in the 'atom_data' dictionary we
            are currently at (in terms of keys in the 'atom_data'
            dictionary that we need to traverse in order to reach the
            current point).

        oldids2newids : ``dict``
            A dictionary mapping the old identifiers to the new
            identifiers, in case we are renaming items.

        Returns
        -------
        updated_items : ``dict``
            A dictionary containing the updated items and their
            attributes.
        """

        # If we are checking models
        if current_depth == 1:

            # The items to check will be the entire structure.
            items_to_check = struct

        # Otherwise
        else:
            
            # The items to check will be the ones in the '_items'
            # dictioanary in the current structure.
            items_to_check = struct["_items"]

        #-------------------------------------------------------------#

        # Initialize an empty dictionary where the updated items will
        # be stored.
        updated_items = {}

        # Initialize the new starting number from the one provided.
        new_start = start

        #-------------------------------------------------------------#
        
        # For each item and associated dictionary in the structure to
        # be checked
        for item, item_dict in items_to_check.items():

            # If the item is a tuple AND we are at the residues'
            # level (we use both conditions in case we introduce
            # tuple identifiers for other levels in he future)
            if isinstance(item, tuple) and current_depth == 6:

                # Assemble the new item by substituting the residue's
                # sequence number with the one resulting from the
                # renumbering, and leaving the residue's insertion
                # code and name unchanged.
                new_item = (new_start, item[1], item[2])

                # Add the new identifier to the dictionary of updated
                # items.
                updated_items[new_item] = item_dict

                # Update the mapping between the old identifiers and
                # the new identifiers (we are sure the path is at least
                # two items long because we are at the residues'
                # level).
                oldids2newids[current_path[0]][\
                    current_path[1:] + (item,)] = \
                        current_path[1:] + (new_item,)

            #---------------------------------------------------------#

            # Otherwise
            else:

                # Add the new identifier to the dictionary of updated
                # items.
                updated_items[new_start] = item_dict

                # If the item is a tuple AND we are at the residues'
                # level (we use both conditions in case we introduce
                # tuple identifiers for other levels in he future)
                if len(current_path) > 2 and current_depth == 6:

                    # Update the mapping between the old identifiers
                    # and the new identifiers (we are sure the path is
                    # at least two items long because we are at the
                    # residues' level).
                    oldids2newids[current_path[0]][\
                        current_path[1:] + (item,)] = \
                            current_path[1:] + (new_start,)

                # If the path is 0 items long AND we are at the models'
                # level (we use both conditions just to be sure, but
                # either should be sufficient)
                elif len(current_path) == 0 and current_depth == 1:

                    # Update the mapping between the old identifiers
                    # and the new identifiers (we are sure the path is
                    # at zero items long because we are at the models'
                    # level).
                    oldids2newids[(item,)] = (new_start,)  

            # Update the starting number.
            new_start += 1

        #-------------------------------------------------------------#

        # If we renumbered models
        if current_depth == 1:

            # Just return the dictionary of updated items.
            return updated_items

        # Otherwise
        else:

            # Return the renumbered items and the current structure's
            # attributes.
            return {"_items" : updated_items,
                    "_attributes" : \
                        {k : v for k, v \
                         in struct["_attributes"].items()}}


    def _modify(self,
                action,
                elements_type,
                selected,
                level = None,
                mapping = None,
                start = None,
                attribute = None,
                new_value = None):
        """Perform various modifications on the 'atom_data' dictionary.

        Parameters
        ----------
        action : ``str``, {``"keep"``, ``"remove"``, ``"rename"``, \
            ``"renumber"``, ``"assign"``}
            The modification to be performed on the 'atom_data'
            dictionary:

            - ``"keep"`` keeps only selected items.
            - ``"remove"`` removes selected items.
            - ``"rename"`` renames selected items, given a mapping
              between the old names and the new names.
            - ``"renumber"`` renumbers selected items, given a new
              starting number.

        elements_type : ``str``, {``"items"``, ``"attributes"``}
            Whether the ``action`` at ``level`` needs to
            be performed on items or attributes.

        selected : ``dict``
            A dictionary of models, chains, segments, residues, and
            atoms considered during the modification, and which
            values specific model, chain, segment, residue, or atom
            attributes the items should have for them to be
            considered.

        level : ``str``, {``"models"``, ``"chains"``, ``"segments"``, \
            ``"residues"``, ``"atoms"``}, optional
            At what level of the hierarchy the ``action`` should be
            performed, if ``action`` is ``"rename"``, ``"renumber"``,
            or ``"assign"``.

        mapping : ``dict``, optional
            A dictionary mapping the old names to the new names,
            when ``action = "rename"``.

        start : ``int``, optional
            The new starting number, when ``action = "renumber"``.

        attribute : ``str``, optional
            The name of the attribute to be renamed, if
            ``action = "rename"`` and ``elements_type = "attributes"``.

        new_value : ``str``, ``int``, ``float``, ``tuple``, \
            optional
            The new value for the attribute, if ``action = "assign"``
            and ``elements_type = "attributes"``. 
        """

        # Define the recursion through the structure.
        def recurse(struct,
                    action,
                    elements_type,
                    selected,
                    target_depth,
                    mapping = None,
                    start = None,
                    attribute = None,
                    new_value = None,
                    current_depth = 1,
                    current_path = (),
                    oldids2newids = None,
                    should_recurse = True):

            # If we are at the target depth
            if current_depth == target_depth:

                # If the target depth is the models' depth
                if current_depth == 1:
                    
                    # The depth of the items we are checking is only
                    # one level down.
                    items_depth = current_depth + 1
                
                # Otherwise
                else:
                    
                    # The depth of the items we are checking is two
                    # levels down.
                    items_depth = current_depth + 2

                # Get the user-selected items.
                selected_items = \
                    selected[self._ITEMSATTRS_DEPTH_TO_LEVELS[\
                        items_depth]]["_items"]

                # Get the user-selected attributes.
                selected_attrs = \
                    selected[self._ITEMSATTRS_DEPTH_TO_LEVELS[\
                        items_depth]]["_attributes"]

                #-----------------------------------------------------#

                # If we are keeping the selected items
                if action == "keep":

                    # Keep the items and return the result.
                    return self._keep(\
                                struct = struct,
                                selected_items = selected_items,
                                selected_attrs = selected_attrs,
                                current_depth = current_depth)
                
                #-----------------------------------------------------#

                # If we are removing the selected items
                elif action == "remove":

                    # Remove the items and return the result.
                    return self._remove(\
                                struct = struct,
                                selected_items = selected_items,
                                selected_attrs = selected_attrs,
                                current_depth = current_depth)

                #-----------------------------------------------------#

                # If we are renaming the selected items or assigning
                # a new value to one of their attributes
                elif action in ("rename", "assign"):

                    # Rename the items and return the result.
                    return self._rename_or_assign(\
                                struct = struct,
                                selected_items = selected_items,
                                selected_attrs = selected_attrs,
                                mapping = mapping,
                                elements_type = elements_type,
                                attribute = attribute,
                                new_value = new_value,
                                current_depth = current_depth,
                                current_path = current_path,
                                oldids2newids = oldids2newids)

                #-----------------------------------------------------#

                # If we are renumbering the selected items
                elif action == "renumber":

                    # Renumber the items and return the result.
                    return self._renumber(\
                                struct = struct,
                                start = start,
                                current_depth = current_depth,
                                current_path = current_path,
                                oldids2newids = oldids2newids)

            #---------------------------------------------------------#

            # If we are not at the target depth yet
            else:

                # If we are at a level containing items
                if current_depth in self._IDS_DEPTH_TO_LEVELS:

                    # For each key and associated sub-structure in the
                    # current structure
                    for key, sub_struct in list(struct.items()):

                        # Get whether we should recurse through the
                        # structure or not.
                        should_recurse = \
                            self._should_recurse(\
                                key = key,
                                value = sub_struct,
                                selected = selected,
                                current_depth = current_depth)

                        # If we should recurse
                        if should_recurse:
                            
                            # Recursively traverse the sub-structure
                            # and get the updated sub-structure.
                            new_sub_struct = \
                                recurse(\
                                    struct = sub_struct,
                                    action = action,
                                    elements_type = elements_type,
                                    selected = selected,
                                    target_depth = target_depth,
                                    mapping = mapping,
                                    start = start,
                                    attribute = attribute,
                                    new_value = new_value,
                                    current_depth = current_depth + 1,
                                    current_path = \
                                        current_path + (key,),
                                    should_recurse = should_recurse,
                                    oldids2newids = oldids2newids)

                            # If the sub-structure does not contain any
                            # items
                            if new_sub_struct["_items"] == {}:

                                # Remove the sub-structure from the
                                # dictionary.
                                struct.pop(key)

                            # Otherwise
                            else:

                                # Substitute the original sub-structure
                                # with the updated one.
                                struct[key] = new_sub_struct

                        # If we should not recurse and the action to be
                        # performed is 'keep'
                        elif not should_recurse and action == "keep":

                            # Remove the sub-structure from the
                            # dictionary (= ignored sub-structures are
                            # sub-structures to be removed if the
                            # action is 'keep').
                            struct.pop(key)

                #-----------------------------------------------------#

                # If we are at a ['_items', '_attributes'] depth
                else:

                    # For each key and associated sub-structure in the
                    # current structure
                    for key, sub_struct in struct.items():

                        # If we are at the '_items'
                        if key == "_items":

                            # Recursively traverse the sub-structure
                            # and substitute it in the structure with
                            # with the resulting sub-structure.
                            struct[key] = \
                                recurse(\
                                    struct = sub_struct,
                                    action = action,
                                    elements_type = elements_type,
                                    selected = selected,
                                    target_depth = target_depth,
                                    mapping = mapping,
                                    start = start,
                                    attribute = attribute,
                                    new_value = new_value,
                                    current_depth = current_depth + 1,
                                    current_path = current_path,
                                    should_recurse = should_recurse,
                                    oldids2newids = oldids2newids)

                        # We ignore the '_attributes' key.

            #---------------------------------------------------------#

            # Return the updated structure.
            return struct

        #-------------------------------------------------------------#
   
        # If no 'level' was passed
        if level is None:

            # Get the hierarchy of the the selected items/attributes,
            # setting 'None' for levels where no items/attributes
            # were selected.
            selected_list = \
                [k if v["_items"] or v["_attributes"] \
                 else None for k, v in selected.items()]

            # Get the target level from the last not-None element of
            # the reversed hierarchy of levels.
            level = \
                next((i for i in reversed(selected_list) \
                      if i is not None))

        # Get the target depth for the recursion from the mapping
        # between depths and corresponding items and the index of
        # our target in the reversed-hierarchy list.
        target_depth = \
            {v : k for k, v \
             in self._ITEMSATTRS_DEPTH_TO_LEVELS.items()}[level] - 2

        # If we ended up with a target depth of 0
        if target_depth == 0:

            # Set it to 1 since we are operating on models.
            target_depth = 1

        #-------------------------------------------------------------#
        
        # If we have to renumber the items
        if action == "renumber":

            # If the starting number is not an integer
            if not isinstance(start, int):

                # Raise an error.
                errstr = "'start' must be an integer."
                raise TypeError(errstr)

            # If a starting number was specified for something other
            # than models or residues
            if level not in ("models", "residues"):

                # Raise an error.
                errstr = \
                    "'renumber' can only be used for models and " \
                    "residues."
                raise ValueError(errstr)

        # If we have to rename the items
        elif action == "rename":

            # If the mapping is not a dictionary
            if not isinstance(mapping, dict):

                # Raise an error.
                errstr = "'mapping' must be a dictionary."
                raise TypeError(errstr)

        # If we have to assign a new value to an attribute
        elif action == "assign":

            # If the attribute's name is not a string
            if not isinstance(attribute, str):

                # Raise an error.
                errstr = "'attribute' must be a string."
                raise TypeError(errstr)

            # If the attribute's new value is not an integer, floating-
            # point number, string, or tuple
            if not isinstance(new_value, (int, float, str, tuple)):

                # Raise an error.
                errstr = \
                    "'new_value' must be an integer, floating-point " \
                    "number, string, or tuple."
                raise TypeError(errstr)

        #-------------------------------------------------------------#

        # Create an empty dictionary to store the mapping between the
        # old identifiers and the new identifiers, to be used if we are
        # renaming or renumbering items.
        oldids2newids = defaultdict(lambda: defaultdict(tuple))

        # Recursively modify the structure.
        self._atom_data = \
            recurse(struct = copy.deepcopy(self.atom_data),
                    action = action,
                    elements_type = elements_type,
                    selected = selected,
                    target_depth = target_depth,
                    mapping = mapping,
                    start = start,
                    attribute = attribute,
                    new_value = new_value,
                    oldids2newids = oldids2newids)

        #-------------------------------------------------------------#

        # Renumber the atoms and update the connectivity data.
        self._update_atom_and_conect(oldids2newids = oldids2newids)


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
            The names of the atoms in residues named ``residue_names``
            in the order they should be sorted.

        other_atoms_position : ``str``, {``"before"``, ``"after"``}
            Where to put the atoms not included in ``atoms_names``, if
            found in the residues named ``residue_name``.

            They can be put ``"before"`` or ``"after"`` the atoms
            specified in ``atoms_names``.
        """

        # Define the recursion through the structure.
        def recurse(struct,
                    residue_name,
                    atoms_names,
                    other_atoms_position,
                    current_depth = 1,
                    current_path = ()):

            # If we have reached the atoms' level
            if current_depth == 7: 

                # Create a new dictionary to store the atoms after the
                # sorting.
                new_atoms_order = {}

                #-----------------------------------------------------#

                # Get the atoms matching the names of those provided.
                expected_atoms = \
                    {atom_serial : atom_data \
                     for atom_serial, atom_data \
                     in struct.items() \
                     if atom_data["_attributes"]["label_atom_id"] \
                        in atoms_names}
                
                # Get the extra atoms in the current residue (= atoms
                # whose names are not included in the list of atoms to
                # sort).
                other_atoms = \
                    {atom_serial : atom_data \
                     for atom_serial, atom_data \
                     in struct.items() \
                     if atom_data["_attributes"]["label_atom_id"] \
                         not in atoms_names}

                #-----------------------------------------------------#

                # Get the atoms that were found in the current residue.
                found_atoms = \
                    [atom_data["_attributes"]["label_atom_id"] \
                     for atom_data in expected_atoms.values()]
                
                # Get the atoms missing from the residues by comparing
                # the set of atoms' names in the list passed and the
                # names of the atoms found in the residue.
                missing_atoms = \
                    set(atoms_names) - set(found_atoms)
                
                # If any atoms are missing
                if missing_atoms:

                    # Get the current model, chain, segment, and
                    # residue.
                    model, chain, segment, residue = \
                        current_path[0], current_path[2], \
                        current_path[4], current_path[5]

                    # Warn the user.
                    warnstr = \
                        f"Missing atoms {', '.join(missing_atoms)} " \
                        f"in model {model}, chain '{chain}', " \
                        f"segment '{segment}', residue '{residue}'."
                    logger.warning(warnstr)

                #-----------------------------------------------------#

                # For each provided atom name
                for atom_name in atoms_names:

                    # For each atom found in the residue
                    for atom_serial, atom_data \
                    in expected_atoms.items():

                        # If the current atom's name matches the
                        # current name
                        if atom_data["_attributes"][\
                            "label_atom_id"] == atom_name:

                            # Add the atom to the dictionary storing
                            # the sorted atoms.
                            new_atoms_order[atom_serial] = atom_data

                #-----------------------------------------------------#

                # If we want to add the extra atoms before the sorted
                # atoms
                if other_atoms_position == "before":

                    # Create a new dictionary with the extra atoms
                    # before the sorted ones.
                    new_atoms_order = \
                        {**other_atoms, **new_atoms_order}
                
                # If we want to add the extra atoms after the sorted
                # atoms.
                elif other_atoms_position == "after":

                    # Create a new dictionary with the extra atoms
                    # after the sorted ones.
                    new_atoms_order = \
                        {**new_atoms_order, **other_atoms}
                
                # If an invalid value was passed
                else:

                    # Raise an error.
                    errstr = \
                        "'other_atoms_position' must be 'before' " \
                        "or 'after'"
                    raise ValueError(errstr)

                # Return the atoms, sorted.
                return new_atoms_order

            #---------------------------------------------------------#
            
            # If we are at the residues' level
            elif current_depth == 6:

                # For each residue and associated data
                for res, res_data in struct["_items"].items():

                    # If the current residue's name matches the type
                    # of residue for which we need to sort the atoms
                    if res[2] == residue_name:

                        # Recurse through the atoms belonging to the
                        # current residue and return them sorted.
                        sorted_atoms = \
                            recurse(\
                                struct = res_data["_items"],
                                residue_name = residue_name,
                                atoms_names = atoms_names,
                                other_atoms_position = \
                                    other_atoms_position, 
                                current_depth = current_depth + 1,
                                current_path = current_path + (res,))
                        
                        # Add the sorted atoms to the current residue,
                        # substituting the ones that were there before.
                        res_data["_items"] = sorted_atoms

            #---------------------------------------------------------#
            
            # Otherwise
            else:

                # For each key and associated sub-structure in the
                # current structure
                for key, sub_struct in struct.items():

                    # If the sub-structure is a dictionary
                    if isinstance(sub_struct, dict):

                        # Recurse through it.
                        recurse(\
                            struct = sub_struct,
                            residue_name = residue_name,
                            atoms_names = atoms_names,
                            other_atoms_position = \
                                 other_atoms_position,
                            current_depth = current_depth + 1,
                            current_path = current_path + (key,))

            #---------------------------------------------------------#

            # Return the updated structure.
            return struct

        #-------------------------------------------------------------#

        # Sort the atoms (we need to pass a copy of the 'atom_data'
        # attribute because it will be modified).
        self._atom_data = \
            recurse(struct = copy.deepcopy(self.atom_data),
                    residue_name = residue_name,
                    atoms_names = atoms_names,
                    other_atoms_position = other_atoms_position)

        #-------------------------------------------------------------#

        # Renumber the atoms and update the connectivity data.
        self._update_atom_and_conect()


    ###################################################################


    def models(self):
        """Get all models in the structure.

        Returns
        -------
        models : ``generator``
            A generator containing the models.

        Examples
        --------
        .. code-block:: python

            # Get all models in a multi-model structure.

            # Import 'pdbcraft'.
            >>> import pdbcraft
            
            # Load the example '2LPC.pdb' - it is a multi-model
            # structure obtained using NMR spectroscopy.
            >>> struct = pdbcraft.load_example("2LPC.pdb")

            # Get a list of all models.
            >>> models = list(struct.models())
        
            # Inspect the list of models.
            >>> models
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
            17, 18, 19, 20, 21]

        See also
        --------
        chains : Get all chains in a model.
        segments : Get all segments in a chain.
        residues : Get all residues in a segment.
        atoms : Get all atoms in a residue.
        """
        
        # Return the models.
        return (item for item \
                    in self.atom_data.keys())


    def chains(self,
               model):
        """Get all chains in a model.

        Parameters
        ----------
        model : ``int``
            The identifier of the model the chains belong to.

        Returns
        -------
        chains : ``generator``
            A generator containing all chains in the selected model.

        Examples
        --------
        .. code-block:: python

            # Get all chains in a specific model.

            # Import 'pdbcraft'.
            >>> import pdbcraft
            
            # Load the example '2LPC.pdb' - it is a multi-model
            # structure obtained using NMR spectroscopy.
            >>> struct = pdbcraft.load_example("2LPC.pdb")

            # Get a list of all chains in model 2.
            >>> chains = list(struct.chains(model = 2))
        
            # Inspect the list of chains.
            >>> chains
            ['A']

        See also
        --------
        models : Get all models in a structure.
        segments : Get all segments in a chain.
        residues : Get all residues in a segment.
        atoms : Get all atoms in a residue.
        """

        # Check the model.
        self._check_items(model = model)
        
        # Return the chains.
        return (item for item \
                    in self.atom_data[model]["_items"].keys())


    def segments(self,
                 model,
                 chain):
        """Get all segments in a chain.

        Parameters
        ----------
        model : ``int``
            The identifier of the model the segments belongs to.

        chain : ``str``
            The identifier of the chain the segments belong to.

        Returns
        -------
        segments : generator
            A generator containing all segments in the selected model
            and chain.

        Examples
        --------
        .. code-block:: python

            # Get all segments in a specific chain.

            # Import 'pdbcraft'.
            >>> import pdbcraft
            
            # Load the example '2LPC.pdb' - it is a multi-model
            # structure obtained using NMR spectroscopy.
            >>> struct = pdbcraft.load_example("2LPC.pdb")

            # Get a list of all segments in chain 'A', model 2.
            >>> segments = list(struct.segments(model = 2,
                                                chain = "A"))
        
            # Inspect the list of segments.
            >>> segments
            ['']

        See also
        --------
        models : Get all models in a structure.
        chains : Get all chains in a model.
        residues : Get all residues in a segment.
        atoms : Get all atoms in a residue.
        """

        # Check the model and the chain.
        self._check_items(model = model,
                          chain = chain)

        # Return the segments.
        return (item for item \
                    in self.atom_data[model]["_items"][\
                        chain]["_items"].keys())

    def residues(self,
                 model,
                 chain,
                 segment):
        """Get all residues in a segment.

        Parameters
        ----------
        model : ``int``
            The identifier of the model the residues belongs to.

        chain : ``str``
            The identifier of the chain the residues belong to.

        segment : ``str``
            The identifier of the segment the residues belong to.

        Returns
        -------
        residues : generator
            A generator containing all residues in the selected model,
            chain, and segment.

        Examples
        --------
        .. code-block:: python

            # Get all residues in a specific segment.

            # Import 'pdbcraft'.
            >>> import pdbcraft
            
            # Load the example '2LPC.pdb' - it is a multi-model
            # structure obtained using NMR spectroscopy.
            >>> struct = pdbcraft.load_example("2LPC.pdb")

            # Get a list of all residues in segment '', chain 'A',
            # model 2.
            >>> residues = list(struct.residues(model = 2,
                                                chain = "A",
                                                segment = ""))
        
            # Inspect the list of residues.
            >>> residues
            [(1, '', 'MET'), (2, '', 'SER'), (3, '', 'GLN'),
            (4, '', 'SER'), (5, '', 'ASN'), (6, '', 'ARG'),
            (7, '', 'GLU'), (8, '', 'LEU'), (9, '', 'VAL'),
            (10, '', 'VAL'), (11, '', 'ASP'), (12, '', 'PHE'),
            (13, '', 'LEU'), ..., ]

        See also
        --------
        models : Get all models in a structure.
        chains : Get all chains in a model.
        segments : Get all segments in a chain.
        atoms : Get all atoms in a residue.
        """

        # Check the model, the chain, and the segment.
        self._check_items(model = model,
                          chain = chain,
                          segment = segment)

        # Return the residues.
        return (item for item \
                    in self.atom_data[model]["_items"][\
                        chain]["_items"][segment]["_items"].keys())

    def atoms(self,
              model,
              chain,
              segment,
              residue):
        """Get all atoms in a residue.

        Parameters
        ----------
        model : ``int``
            The identifier of the model the atoms belong to.

        chain : ``str``
            The identifier of the chain the atoms belong to.

        segment : ``str``
            The identifier of the segment the atoms belong to.

        residue : ``tuple``
            The identifier of the residue the atoms belong to.

            The identifier is a tuple containing the residue's
            number, insertion code, and name.

        Returns
        -------
        atoms : generator
            A generator containing all atoms in the selected model,
            chain, segment, and residue.

        Examples
        --------
        .. code-block:: python

            # Get all atoms in a specific residue.

            # Import 'pdbcraft'.
            >>> import pdbcraft
            
            # Load the example '2LPC.pdb' - it is a multi-model
            # structure obtained using NMR spectroscopy.
            >>> struct = pdbcraft.load_example("2LPC.pdb")

            # Get a list of all atoms in residue (3, '', 'GLN'),
            # segment '', chain 'A', model 2. The atoms' numeric
            # identifiers will be returned.
            >>> atoms = list(struct.atoms(model = 2,
                                          chain = "A",
                                          segment = "",
                                          residue = (3, "", "GLN")))
        
            # Inspect the list of atoms.
            >>> atoms
            [31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
            45, 46, 47]

        See also
        --------
        models : Get all models in a structure.
        chains : Get all chains in a model.
        segments : Get all segments in a chain.
        residues : Get all residues in a segment.
        """

        # Check the model, the chain, the segment, and the residue.
        self._check_items(model = model,
                          chain = chain,
                          segment = segment,
                          residue = residue)

        # Return the atoms.
        return (item for item \
                    in self.atom_data[model]["_items"][\
                        chain]["_items"][segment]["_items"][\
                            residue]["_items"].keys())


    def atom_serial(self,
                    path):
        """Get the serial number of a specific atom, given the 'path'
        to the atom in the structure in terms of the model, chain,
        segment, and residue the atom belongs to, plus the atom's name.

        This method assumes that, within each residue, the names of
        the atoms are unique.

        Parameters
        ----------
        path : ``tuple``
            The 'path' to the atom (the model, chain, segment, and
            residue it belongs to, plus the atom's name).

            In the 'path':
            - The model is identified by its number (an integer).
            - The chain is identified by its ID (a string).
            - The segment is identified by its ID (a string)
            - The residue is identified by a tuple containing the
            residue's sequence number, insertion code, and name.

        Returns
        -------
        atom_serial : ``int```
            The serial number of the atom the 'path' points to.

        Examples
        --------
        .. code-block:: python

            # Get the serial number of a specific atom, given its
            # 'path' in the structure.

            # Import 'pdbcraft'.
            >>> import pdbcraft
            
            # Load the example '1L9Z.cif' - it is a single-model
            # structure obtained using X-ray crystallography.
            >>> struct = pdbcraft.load_example("1L9Z.cif")

            # Set the 'path' to the CA atom of residue
            # (38, '', 'ASN'), segment '', chain 'B'.
            >>> atom_path = (1, "B", "", (38, "", "ASN"), "CA")

            # Get the serial number of the atom from its 'path'.
            >>> atom_serial = struct.atom_serial(path = atom_path)
        
            # Inspect the atom's serial number.
            >>> atom_serial
            1482

        See also
        --------
        models : Get all models in a structure.
        chains : Get all chains in a model.
        segments : Get all segments in a chain.
        residues : Get all residues in a segment.
        atoms : Get all atoms in a residue.
        """

        # Return the serial number of the atom the 'path' points to.
        return self._get_atom_serial(path = path)


    ###################################################################


    def get_copy(self):
        """Get a copy of the current structure.

        Returns
        -------
        ``struct_copy`` : :class:`pdbcraft.structure.Structure`
            A copy of the current structure.

        Examples
        --------
        .. code-block:: python

            # Get a copy of a structure.

            # Import 'pdbcraft'.
            >>> import pdbcraft
            
            # Load the example '1L9Z.cif' - it is a single-model
            # structure obtained using X-ray crystallography.
            >>> struct = pdbcraft.load_example("1L9Z.cif")

            # Create a copy of the structure.
            >>> struct_copy = struct.get_copy()

            # Verify the equality of the two structures.
            >>> struct == struct_copy
            True

        See also
        --------
        merge : Merge two structures.
        """
        
        # We 'deepcopy' each element of the structure and return a new
        # new 'Structure' because just 'deepcopying' the 'Structure'
        # raises errors.
        return self.__class__(\
                atom_data = copy.deepcopy(self.atom_data),
                conect_data = copy.deepcopy(self.conect_data),
                name = copy.deepcopy(self.name))


    ###################################################################


    def merge(self,
              other):
        """Merge the current structure with another structure, and
        return the merged structure.

        Parameters
        ----------
        other : :class:`pdbcraft.structure.Structure`
            The other structure.

        Returns
        -------
        merged : :class:`pdbcraft.structure.Structure`
            The merged structure.

        Examples
        --------
        .. code-block:: python

            # Merge two structures.

            # Import 'pdbcraft'.
            >>> import pdbcraft

            # Load the example '6X8O_AB.pdb' - it is a single-model
            # structure obtained using X-ray crystallography.
            >>> struct1 = pdbcraft.load_example("6X8O_AB.pdb")

            # Load the example '6X8O_C.pdb' - it is a single-model
            # structure obtained using X-ray crystallography.
            >>> struct2 = pdbcraft.load_example("6X8O_C.pdb")

            # Create a new structure by merging the two structures.
            >>> merged = struct1.merge(other = struct2)

            # Get the chains in the new structure.
            >>> list(merged.chains(model = 1))
            ['A', 'B', 'C']

        See also
        --------
        get_copy : Get a copy of a structure.
        """

        # Merge the structures and return the resulting structure.
        return self._merge(other = other)


    ###################################################################


    def add_bonds(self,
                  bonds,
                  bonds_attributes = None,
                  in_place = False):
        """Add bonds between pairs of atoms, updating the connectivity
        data of the structure while preserving atomic coordinates.

        Parameters
        ----------
        bonds : ``list``
            A list of tuples, each representing a bond between two
            atoms. The tuple must contain the first atom, the
            second atom, the bond's type and the bond's order.

            An atom must be represented by its "path", namely a
            tuple containing the model, chain, segment, and residue
            the atom belongs to, plus the atom's name. This is
            based on the assumption that, within each residue,
            atoms must have unique names.

            For instance, a single covalent bond between
            the N atom of MET 1 in chain A, segment '' and
            the CA atom of the same residue must be represented as
            ``(("A", "", (1, "", "MET"), "N"), 
            ("A", "", (1, "", "MET"), "CA"), "covale", "sing")``.

            Supported bond types are:
            - ``"covale"`` for covalent bonds.
            - ``"disulf"`` for disulfide bridges.
            - ``"hydrog"`` for hydrogen bonds.
            - ``"metalc"`` for metal coordination.

            Supported bond orders are:
            - ``"sing"`` for single bonds.
            - ``"doub"`` for double bonds.
            - ``"trip"`` for triple bonds.
            - ``"quad"`` for quadruple bonds.

            Each bond should be specified only once since,
            for each bond specified, the connectivity data for
            both atoms involved in the bond will be updated.

        bonds_attributes : ``dict``, optional
            A dictionary mapping each bond in ``bonds`` to a dictionary
            of attributes that will be assigned to the bond.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new ``Structure`` will be returned.

        Returns
        -------
        :class:`pdbcraft.structure.Structure` if ``in_place`` is ``False``;
        otherwise, ``None``.

        Examples
        --------
        .. code-block:: python

            # Add bonds between pairs of atoms.

            # Import 'pdbcraft'.
            >>> import pdbcraft

            # Load the example '6X8O_C_noconect.pdb' - it is a
            # single-model structure obtained using X-ray
            # crystallography.
            >>> struct = pdbcraft.load_example("6X8O_C_noconect.pdb")

            # Get the connectivity data of the structure.
            >>> struct.conect_data
            {}

            # Define the bond(s) to be added.
            >>> bonds = \
                    [(("C", "", (166, "", "ARG"), "C"),
                     ("C", "", (167, "", "NH2"), "N"),
                     "covale",
                     "sing")]

            # Add the bond(s) to the structure.
            >>> struct.add_bonds(bonds = bonds,
                                 in_place = True)

            # Get the connectivity data of the structure.
            >>> struct.conect_data
            {1: {('C', '', (166, '', 'ARG'), 'C'):
                    {('C', '', (167, '', 'NH2'), 'N'):
                        {'conn_type_id': 'covale',
                         'pdbx_value_order': 'sing'}},
                 ('C', '', (167, '', 'NH2'), 'N'):
                    {('C', '', (166, '', 'ARG'), 'C'):
                        {'conn_type_id': 'covale',
                         'pdbx_value_order': 'sing'}}}}

        See also
        --------
        add_atom : Add an atom to the structure.
        """

        # Set the keyword arguments to pass to the internal method.
        kwargs = \
            {"bonds" : bonds,
             "bonds_attributes" : bonds_attributes}

        # If the structure needs to be modified in place
        if in_place:

            # Add the bonds to the structure.
            self._add_bonds(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure.
            struct_new = self.get_copy()

            # Add the bonds to the new structure.
            struct_new.add_bonds(**kwargs)

            # Return the new structure.
            return struct_new


    ###################################################################


    def add_atom(self,
                 atom_path,
                 atom_attributes,
                 in_place = False):
        """Add a single atom to the structure, given the 'path' to a
        point in the structure where the atom should be added.

        The 'path' comprises the identifier of the chain, segment, and
        residue the atom should be added to.

        Parameters
        ----------
        atom_path : ``tuple``
            The path to the point in the structure where the atom
            should be added.

        atom_attributes : ``dict``
            A dictionary containing the attributes for the atom.
            
            The dictionary must contain at least the following
            keys:

            - 'label_atom_id', associated with the atom's name.
            - 'label_alt_id', associated with the atom's alternate
              location.
            - 'label_entity_id', associated with the atom's
              entity ID.
            - 'type_symbol', associated with the atom's type.
            - 'Cartn_x', associated with the X coordinate of the
              atom's position.
            - 'Cartn_y', associated with the Y coordinate of the
              atom's position.
            - 'Cartn_z', associated with the Z coordinate of the
              atom's position.
            - 'pdbx_formal_charge', associated with the atom's
              charge.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new structure will be returned.

        Returns
        -------
        :class:`pdbcraft.structure.Structure` if ``in_place`` is ``False``;
        otherwise, ``None``.

        Examples
        --------
        .. code-block:: python

            # Add bonds between pairs of atoms.

            # Import 'pdbcraft'.
            >>> import pdbcraft

            # Load the example '19G.pdb' - it is a single-model
            # structure obtained using X-ray crystallography.
            >>> struct = pdbcraft.load_example("19G.pdb")

            # Get the number of atoms in residue 19G, chain 'B',
            # segment ''.
            >>> len(list(struct.atoms(model = 1,
                                      chain = "B",
                                      segment = "",
                                      residue = (1, "", "19G"))))
            43

            # Define the attribute of the atom to be added.
            >>> atom_attributes = \
                    {"label_atom_id" : "HAC",
                     "label_alt_id" : "",
                     "label_entity_id" : "1",
                     "type_symbol" : "H",
                     "Cartn_x" : 34.300,
                     "Cartn_y" : -8.500,
                     "Cartn_z" : 17.900,
                     "pdbx_formal_charge" : "0",
                     "occupancy" : 1.0,
                     "B_iso_or_equiv" : 0.0}

            # Add the missing atom to the structure.
            >>> atom_path = ("B", "", (1, "", "19G"))
            >>> struct.add_atom(atom_path = atom_path,
                                atom_attributes = atom_attributes,
                                in_place = True)

            # Get the number of atoms in residue 19G, chain 'B',
            # segment ''.
            >>> len(list(struct.atoms(model = 1,
                                      chain = "B",
                                      segment = "",
                                      residue = (1, "", "19G"))))

            # Define the bond(s) that need to be added.
            >>> bonds = \
                    [(("B", "", (1, "", "19G"), "OAC"),
                      ("B", "", (1, "", "19G"), "HAC"),
                      "covale",
                      "sing")]

            # Add the bond(s) to the structure.
            >>> struct.add_bonds(bonds = bonds,
                                 in_place = True)

        See also
        --------
        add_bonds : Add bonds between pairs of atoms in the structure.
        
        """

        # Set the keyword arguments to pass to the internal method.
        kwargs = \
            {"atom_path" : atom_path,
             "atom_attributes" : atom_attributes}

        # If the structure needs to be modified in place
        if in_place:
            
            # Add the atom to the structure.
            self._add_atom(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure.
            struct_new = self.get_copy()

            # Add the atom to the new structure.
            struct_new._add_atom(**kwargs)

            # Return the new structure.
            return struct_new


    ###################################################################


    def get_items(self,
                  level,
                  action = "get_unique",
                  squeeze = None,
                  models = None,
                  models_attributes = None,
                  chains = None,
                  chains_attributes = None,
                  segments = None,
                  segments_attributes = None,
                  residues = None,
                  residues_attributes = None,
                  atoms = None,
                  atoms_attributes = None):
        """Get the identifiers of specific items.

        Parameters
        ----------
        level : ``str``, {``"models"``, ``"chains"``, \
            ``"segments"``, ``"residues"``, ``"atoms"``}
            At which level of the hierarchy we should retrieve items.

        squeeze : ``str``, {``"models"``, ``"chains"``, \
            ``"segments"``, ``"residues"``, ``"atoms"``}, optional
            At which level of the hierarchy items should be 'squeezed',
            meaning that the 'get' action will involve all values found
            at deeper levels and consider them all together.

            For instance, if ``action = "get_unique"``, 
            ``level = "residues"``, and ``squeeze = "chains"``, the
            result will be a dictionary mapping each model to the
            chains it contains, and each chain will be mapped to a list
            containing all unique residues found in the chain. This
            will ignore the segments they residues belong to, meaning
            that if residues with the same ID are found in different
            segments, they will be reported only once.

            If ``squeeze`` is not specified, it will automatically be
            set to the level above ``level``.

            The hierarchy goes: models > chains > segments > residues >
            atoms.

        action : ``str``, {``"get_unique"``, ``"get_min"``, \
            ``"get_max"``}, ``"get_unique"`` 
            The 'get' action to be performed:
            - ``"get_unique"`` retrieves unique identifiers.
            - ``"get_min"`` retrieves the minimum among the
            identifiers.
            - ``"get_max"`` retrieved the maximum among the
            identifiers.

        models : ``list``, optional
            The unique models to be retrieved if ``level = "models"``.

            Otherwise, the models from which the selected chains,
            segments, residues, or atoms will be retrieved. If not
            provided, all models will be considered.

            Each model is defined by its number.

        models_attributes : ``dict``, optional
            A dictionary mapping names of models' attributes to lists
            of their accepted values.

            Only models whose attributes' values are among the accepted
            ones will be considered.

        chains : ``list``, optional
            The unique models to be retrieved if ``level = "chains"``.

            Otherwise, the chains from which the selected segments,
            residues, or atoms will be retrieved. If not provided, all
            chains in the selected models will be considered.

            Each chain is defined by its identifier.

        chains_attributes : ``dict``, optional
            A dictionary mapping names of chains' attributes to lists
            of their accepted values.

            Only chains whose attributes' values are among the accepted
            ones will be considered.

        segments : ``list``, optional
            The unique segments to be retrieved if 
            ``level = "segments"``.

            Otherwise, the segments from which the selected residues or
            atoms will be retrieved. If not provided, all segments in
            the selected models and chains will be considered.

            Each segments is defined by its identifier.

        segments_attributes : ``dict``, optional
            A dictionary mapping names of segments' attributes to lists
            of their accepted values.

            Only segments whose attributes' values are among the
            accepted ones will be considered.

        residues : ``list``, optional
            The unique residues to be retrieved if 
            ``level = "residues"``.

            Otherwise, the residues from which the selected atoms will
            be retrieved. If not provided, all residues will be
            considered.

            Each residue is defined by either its full identifier (a
            tuple containing the residues' sequence number, insertion
            code, and name), its sequence number, or its name.

        residues_attributes : ``dict``, optional
            A dictionary mapping names of residues' attributes to lists
            of their accepted values.

            Only residues whose attributes' values are among the 
            accepted ones will be considered in the selected segments,
            chains, and models.

        atoms : ``list``, optional
            The unique atoms to be retrieved if ``level = "atoms"``.
            
            Each atom is defined by its serial number.

        atoms_attributes : ``dict``, optional
            A dictionary mapping names of atoms' attributes to lists
            of their accepted values.

            Only atoms whose attributes' values are among the accepted
            ones will be considered for retrieval in the selected
            models, chains, segments, and residues.

        Examples
        --------
        .. code-block:: python

            # Get the identifiers of specific items.

            # Import 'pdbcraft'.
            >>> import pdbcraft

            # Load the example '6X8O_AB.pdb' - it is a single-model
            # structure obtained using X-ray crystallography.
            >>> struct = pdbcraft.load_example("6X8O_AB.pdb")

            # Get all residues in a specific chain, considering
            # all segments.
            #
            # The result is a nested dictionary mapping each
            # model to the chain of interest, which is then
            # in turn mapped to the segments included in
            # the chain, which are in turn mapped to the
            # residues in them.
            >>> struct.get_items(level = "residues",
                                 chains = ["A"])
            {1: {'A': {'': [(141, '', 'ASP'), (142, '', 'MET'),
            (143, '', 'ARG'), (144, '', 'PRO'), (145, '', 'GLU'),
            ..., ]}}}

            # If we want to get rid of the segments' level
            # and list all the residues in a chain regardless of
            # the segment they belong to, we use the 'squeeze'
            # keyword with the 'chains' option.
            >>> struct.get_items(level = "residues",
                                 chains = ["A"],
                                 squeeze = "chains")
            {1: {'A': [(141, '', 'ASP'), (142, '', 'MET'),
            (143, '', 'ARG'), (144, '', 'PRO'), (145, '', 'GLU'),
            ..., ]}}

            # We can also get all residues having a specific
            # attribute, such as residues marked for being
            # heteroresidues.
            #
            # In this case, we use the 'residues_attributes'
            # keyword with a dictionary containing the
            # attribute(s) of choice and the accepted values
            # for such attributes.
            #
            # The 'is_het' attribute is automatically defined
            # for each residue when parsing the structure from
            # the PDB file (using the HETATM records).
            >>> res_attrs = {"is_het" : [True]}
            >>> struct.get_items(level = "residues",
                                 chains = ["A"],
                                 squeeze = "chains",
                                 residues_attributes = res_attrs))
            {1: {'A': [(167, '', 'NH2'), (201, '', 'HOH'),
            (202, '', 'HOH'), (203, '', 'HOH'), ..., ]}}

        See also
        --------
        get_attribute : Get the unique values, minimum value, or \
        maximum value of a specific attribute.
        """

        # Set the keyword arguments to pass to the internal method.
        kwargs = \
            {"action" : action,
             "level" : level,
             "squeeze" : squeeze,
             "elements_type" : "items",
             "selected" : \
                self._get_selected(\
                        models = models,
                        models_attributes = models_attributes,
                        chains = chains,
                        chains_attributes = chains_attributes,
                        segments = segments,
                        segments_attributes = segments_attributes,
                        residues = residues,
                        residues_attributes = residues_attributes,
                        atoms = atoms,
                        atoms_attributes = atoms_attributes)}

        # Return the result.
        return self._get(**kwargs)


    def get_attribute(self,
                      attribute,
                      level,
                      action = "get_unique",
                      squeeze = None,
                      models = None,
                      models_attributes = None,
                      chains = None,
                      chains_attributes = None,
                      segments = None,
                      segments_attributes = None,
                      residues = None,
                      residues_attributes = None,
                      atoms = None,
                      atoms_attributes = None):
        """Get the unique values, minimum value, or maximum value
        of a specific attribute.

        Parameters
        ----------
        attribute : ``str``
            The name of the attribute.

        level : ``str``, {``"models"``, ``"chains"``, \
            ``"segments"``, ``"residues"``, ``"atoms"``}
            At which level of the hierarchy the attribute belongs to.

        action : ``str``, {``"get_unique"``, ``"get_min"``, \
            ``"get_max"``}, ``"get_unique"``
            The 'get' action to be performed:
            - ``"get_unique"`` retrieves all unique values for the
            attribute.
            - ``"get_min"`` retrieves the minimum value of the
            attribute.
            - ``"get_max"`` retrieved the maximum value of the
            attribute.

        squeeze : ``str``, {``"models"``, ``"chains"``, \
            ``"segments"``, ``"residues"``, ``"atoms"``}, optional
            At which level of the hierarchy items should be 'squeezed',
            meaning that the 'get' action will involve all values found
            at deeper levels and consider them all together.

            For instance, if ``action = "get_unique"``, 
            ``level = "residues"``, and ``squeeze = "chains"``, the
            result will be a dictionary mapping each model to the
            chains it contains, and each chain will be mapped to a list
            containing all unique residues found in the chain. This
            will ignore the segments they residues belong to, meaning
            that if residues with the same ID are found in different
            segments, they will be reported only once.

            If ``squeeze`` is not specified, it will automatically be
            set to the level above ``level``.

            The hierarchy goes: models > chains > segments > residues >
            atoms.

        models : ``list``, optional
            If ``level = "models"``, the models whose ``attribute``
            values need to be retrieved.

            Otherwise, the models from which chains, segments,
            residues, and atoms will be selected. If not provided,
            all models will be considered.

            Each model is defined by its number.

        models_attributes : ``dict``, optional
            A dictionary mapping names of models' attributes to lists
            of their accepted values.

            Only items in models whose attributes' values are among the
            accepted ones will be considered.

            This option is ignored if ``level = "models"``.

        chains : ``list``, optional
            If ``level = "chains"``, the chains whose ``attribute``
            values need to be retrieved.

            Otherwise, the chains from which segments, residues, or
            atoms will be selected. If not provided, all chains in the
            selected models will be considered.

            Each chain is defined by its identifier.

            This option is ignored if ``level = "models"``.

        chains_attributes : ``dict``, optional
            A dictionary mapping names of chains' attributes to lists
            of their accepted values.

            Only items in chains whose attributes' values are among the
            accepted ones will be considered.

            This option is ignored if ``level = "chains"`` or above.

        segments : ``list``, optional
            If ``level = "segments"``, the segments whose ``attribute``
            values need to be retrieved.

            Otherwise, the segments from which residues or atoms will
            be selected. If not provided, all segments in the selected
            models and chains will be considered.

            Each segment is defined by its identifier.

            This option is ignored if ``level = "chains"`` or above.

        segments_attributes : ``dict``, optional
            A dictionary mapping names of segments' attributes to lists
            of their accepted values.

            Only items in segments whose attributes' values are among
            the accepted ones will be considered.

            This option is ignored if ``level = "segments"`` or above.

        residues : ``list``, optional
            If ``level = "residues"``, the residues whose ``attribute``
            values need to be retrieved.  

            Otherwise, the residues from which the selected atoms will
            be selected. If not provided, all residues will be
            considered.

            Each residue is defined by either its full identifier (a
            tuple containing the residues' sequence number, insertion
            code, and name), its sequence number, or its name.

            This option is ignored if ``level = "segments"`` or above.

        residues_attributes : ``dict``, optional
            A dictionary mapping names of residues' attributes to lists
            of their accepted values.

            Only residues whose attributes' values are among the
            accepted ones will be considered in the selected segments,
            chains, and models.

            This option is ignored if ``level = "residues"`` or above.

        atoms : ``list``, optional
            If ``level = "atoms"``, the atoms whose ``attribute``
            values need to be retrieved.  
            
            Each atom is defined by its serial number.

        Examples
        --------
        .. code-block:: python

            # Get the values of a specific attribute.

            # Import 'pdbcraft'.
            >>> import pdbcraft

            # Load the example '6X8O_AB.pdb' - it is a single-model
            # structure obtained using X-ray crystallography.
            >>> struct = pdbcraft.load_example("6X8O_AB.pdb")

            # Get whether each segment in chain 'A' contains both
            # canonical residues and heteroresidues.
            #
            # To do this, we use the 'is_het' attribute, which is
            # automatically defined for each residue when parsing
            # the structure from the PDB file (using the HETATM
            # records).
            #
            # The result is a nested dictionary mapping each
            # model to the chain of interest, which is then
            # in turn mapped to the segments included in
            # the chain, which are in turn mapped to the unique
            # values of the 'is_het' attribute found for the
            # residues in each segment. If both 'True' and
            # 'False' are found, it means that the segment
            # contains both canonical residues and
            # heteroresidues. If only 'True' is found, it
            # means that the segment contains only
            # heteroresidues. If only 'False' is found, it means
            # that the segment contains only canonical residues.
            >>> struct.get_attribute(attribute = "is_het", 
                                     level = "residues",
                                     chains = ["A"])
            {1: {'A': {'': [False, True]}}}

            # Get whether chain 'A' contains both canonical
            # residues and heteroresidues. To skip the segments'
            # level, we use the 'squeeze' keyword with the
            # 'chains' option.
            >>> struct.get_attribute(attribute = "is_het",
                                     level = "residues",
                                     chains = ["A"],
                                     squeeze = "chains")
            {1: {'A': [False, True]}}

        See also
        --------
        get_items : Get the identifiers of specific items.
        """

        # Set the keyword arguments to pass to the internal method.
        kwargs = \
            {"action" : action,
             "level" : level,
             "squeeze" : squeeze,
             "elements_type" : "attributes",
             "attribute" : attribute,
             "selected" : \
                self._get_selected(\
                        models = models,
                        models_attributes = models_attributes,
                        chains = chains,
                        chains_attributes = chains_attributes,
                        segments = segments,
                        segments_attributes = segments_attributes,
                        residues = residues,
                        residues_attributes = residues_attributes,
                        atoms = atoms,
                        atoms_attributes = None)}

        # Return the result.
        return self._get(**kwargs)


    ###################################################################


    def keep_items(self,
                   models = None,
                   models_attributes = None,
                   chains = None,
                   chains_attributes = None,
                   segments = None,
                   segments_attributes = None,
                   residues = None,
                   residues_attributes = None,
                   atoms = None,
                   atoms_attributes = None,
                   in_place = False):
        """Keep only selected items in the structure.

        Parameters
        ----------
        models : ``list``, optional
            The models to be kept.

            Each model is defined by its number.

            If not provided, all models will be kept.

        models_attributes : ``dict``, optional
            A dictionary mapping names of models' attributes to lists
            of their accepted values.

            Only models whose attributes' values are among the accepted
            ones will be kept.

        chains : ``list``, optional
            The chains to be kept in the selected models.
            
            Each chain is defined by its chain identifier.

            If not provided, all chains will be kept.

        chains_attributes : ``dict``, optional
            A dictionary mapping names of chains' attributes to lists
            of their accepted values.

            Only chains whose attributes' values are among the accepted
            ones will be kept in the selected models.

        segments : ``list``, optional
            The segments to be kept in the selected models and chains.
            
            Each segment is defined by its segment identifier.

            If not provided, all segments will be kept.

        segments_attributes : ``dict``, optional
            A dictionary mapping names of segments' attributes to lists
            of their accepted values.

            Only segments whose attributes' values are among the
            accepted ones will be kept in the selected models and
            chains.

        residues : ``list``, optional
            The residues to be kept in the selected models, chains,
            and segments.

            Each residue is defined by either its full identifier (a
            tuple containing the residues' sequence number, insertion
            code, and name), its sequence number, or its name.

            The reserved keywords ``"protein"``, ``"dna"``, and
            ``"rna"`` can be used instead of the residues' identifiers,
            sequence numbers, or names, to keep only protein, DNA,
            and RNA residues, respectively. The reserved keyword
            ``"het"`` can also be used to keep only heteroresidues.

            If not provided, all residues will be kept.

        residues_attributes : ``dict``, optional
            A dictionary mapping names of residues' attributes to lists
            of their accepted values.

            Only residues whose attributes' values are among the
            accepted ones will be kept in the selected segments,
            chains, and models.

        atoms : ``list``, optional
            The atoms to be kept in the selected models, chains,
            segments, and residues.
            
            Each atom is defined by its serial number.
            
            If not provided, all atoms will be kept.

        atoms_attributes : ``dict``, optional
            A dictionary mapping names of atoms' attributes to lists
            of their accepted values.

            Only atoms whose attributes' values are among the accepted
            ones will be kept in the selected models, chains, segments,
            and residues.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new structure will be returned.

        Returns
        -------
        :class:`pdbcraft.structure.Structure` if ``in_place`` is ``False``;
        otherwise, ``None``.

        Examples
        --------
        .. code-block:: python

            # Keep only selected items in a structure.

            # Import 'pdbcraft'.
            >>> import pdbcraft
            
            # Load the example '3B9F.cif' - it is a single-model
            # structure obtained using X-ray crystallography.
            >>> struct = pdbcraft.load_example("3B9F.cif")

            # Get the structure's chains.
            >>> list(struct.chains(model = 1))
            ['L', 'H', 'I', 'A', 'B']

            # Keep only chains 'L' and 'H' and create a new structure
            # containing only these chains
            >>> chains = ["L", "H"]
            >>> struct_chainsLH = struct.keep_items(chains = chains)

            # Get the chains in the new structure.
            >>> list(struct_chainsLH.chains(model = 1))
            ['L', 'H']

            # Get the residues in chain 'L', segment ''.
            >>> list(struct_chainsLH.residues(model = 1,
                                              chain = "L",
                                              segment = ""))
            [(1, 'R', 'SER'), (1, 'Q', 'GLU'), (1, 'P', 'TYR'),
            (1, 'O', 'GLN'), (1, 'N', 'THR'), (1, 'M', 'PHE'),
            (1, 'L', 'PHE'), (1, 'K', 'ASN'), (1, 'J', 'PRO'),
            (1, 'I', 'ARG'), (1, 'H', 'THR'), (1, 'G', 'PHE'),
            (1, 'F', 'GLY'), (1, 'E', 'SER'), (1, 'D', 'GLY'),
            (1, 'C', 'GLU'), (1, 'B', 'ALA'), (1, 'A', 'ASP'),
            ..., ]

            # Keep only residues (1, 'A', 'SER') and (14, 'A', 'LYS')
            # in chain 'L', segment '' in the new structure. Since the
            # structure contains residues with the same sequence number
            # but different insertion codes, we must define the
            # residues using their full identifiers.
            >>> residues = [(1, "A", "ASP"), (14, "A", "LYS")]
            >>> struct_chainsLH.keep_items(chains = ["L"],
                                       residues = residues,
                                       in_place = True)

            # Get the residues now in chain 'L', segment ''.
            >>> list(struct_chainsLH.residues(model = 1,
                                              chain = "L",
                                              segment = "")))
            [(1, 'A', 'ASP'), (14, 'A', 'LYS')]
        
        See also
        --------
        remove_items : Remove items (models, chains, segments, \
        residues, or atoms) from a structure.
        rename_items : Rename items (models, chains, segments, \
        residues, or atoms) in a structure.
        renumber_items : Renumber items (models and residues) \
        in a structure.
        """

        # Set the keyword arguments to pass to the internal method.
        kwargs = \
            {"action" : "keep",
             "elements_type" : "items",
             "selected" : \
                self._get_selected(\
                        models = models,
                        models_attributes = models_attributes,
                        chains = chains,
                        chains_attributes = chains_attributes,
                        segments = segments,
                        segments_attributes = segments_attributes,
                        residues = residues,
                        residues_attributes = residues_attributes,
                        atoms = atoms,
                        atoms_attributes = atoms_attributes)}

        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected items.
            self._modify(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure.
            struct_new = self.get_copy()

            # Keep only the selected items in the new structure.
            struct_new._modify(**kwargs)

            # Return the new structure.
            return struct_new


    ###################################################################


    def remove_items(self,
                     models = None,
                     models_attributes = None,
                     chains = None,
                     chains_attributes = None,
                     segments = None,
                     segments_attributes = None,
                     residues = None,
                     residues_attributes = None,
                     atoms = None,
                     atoms_attributes = None,
                     extra_residues = None,
                     in_place = False):
        """Remove selected items.

        Parameters
        ----------
        models : ``list``, optional
            The models to be removed if no chains, segments, residues,
            atoms, or any of their attributes are specified.

            Otherwise, the models from which the selected chains,
            segments, residues, and atoms will be removed.

            Each model is defined by its number.

            If not provided, all models will be considered.

        models_attributes : ``dict``, optional
            A dictionary mapping names of models' attributes to lists
            of values the attributes can have.

            Only items in models whose attributes' values are among
            the ones provided will be considered for removal.

        chains : ``list``, optional
            The chains to be removed if no segments, residues, atoms,
            or any of their attributes are specified.

            Otherwise, the chains from which the selected segments,
            residues, and atoms will be removed.
            
            Each chain is defined by its chain identifier.

            If not provided, all chains in the selected models will be
            considered.

        chains_attributes : ``dict``, optional
            A dictionary mapping names of chains' attributes to lists
            of values the attributes can have.

            Only items in chains whose attributes' values are among the
            ones provided will be considered for removal.

        segments : ``list``, optional
            The segments to be removed if no residues, atoms, or any
            of their attributes are specified.

            Otherwise, the segments from which the selected residues
            and atoms will be removed.
            
            Each segment is defined by its segment identifier.

            If not provided, all segments in the selected chains and
            models will be considered.

        segments_attributes : ``dict``, optional
            A dictionary mapping names of segments' attributes to lists
            of values the attributes can have.

            Only items in segments whose attributes' values are among
            the ones provided will be considered for removal.

        residues : ``list``, optional
            The residues to be removed if no atoms or any of their
            attributes are specified.

            Otherwise, the residues from which the selected atoms will
            be removed.

            The reserved keywords ``"protein"``, ``"dna"``, and
            ``"rna"`` can be used instead of the residues' identifiers,
            sequence numbers, or names, to remove protein, DNA, and RNA
            residues, respectively. The reserved keyword ``"het"`` can
            also be used to remove heteroresidues.

            If not provided, all residues in the selected models,
            chains, and segments will be considered.

        residues_attributes : ``dict``, optional
            A dictionary mapping names of residues' attributes to lists
            of values the attributes can have.

            Only atoms in residues whose attributes' values are among
            the ones provided will be considered for removal.

        atoms : ``list``, optional
            The atoms to be removed in the selected models, chains,
            segments, and residues.
            
            Each atom is defined by its serial number.
            
            If not provided, all atoms in the selected models, chains,
            segments, and residues will be considered.

        atoms_attributes : ``dict``, optional
            A dictionary mapping names of atoms' attributes to lists
            of values they can have.

            Only atoms (in the selected models, chains, segments, and
            residues) whose attributes' values are among the ones
            provided will be considered for removal.

        extra_residues : ``dict``, optional
            A dictionary of names of residues that are not part of
            the canonital sets of protein, DNA, and RNA residues,
            but that should be considered as such. The dictionary
            accepts three keywords:

            - ``"protein"``, mapped to a list of names of extra
              protein residues.
            - ``"dna"``, mapped to a list of names of extra DNA
              residues.
            - ``"rna"``, mapped to a list of names of extra RNA
              residues.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new structure will be returned.

        Returns
        -------
        :class:`pdbcraft.structure.Structure` if ``in_place`` is ``False``;
        otherwise, ``None``.

        Examples
        --------
        .. code-block:: python

            # Remove selected items from a structure.

            # Import 'pdbcraft'.
            >>> import pdbcraft

            # Load the example '2M04.cif' - it is a multi-model
            # structure obtained using NMR spectroscopy.
            >>> struct = pdbcraft.load_example("2M04.cif")

            # Get the structure's chains in model 1.
            >>> list(struct.chains(model = 1))
            ['A', 'B']

            # Remove chain 'A' from the structure. The chain is
            # removed from all models.
            >>> chains = ["A"]
            >>> struct_chainB = struct.remove_items(chains = ["A"])

            # Get the chains in the new structure.
            >>> list(struct_chainB.chains(model = 1))
            ['B']

            # Get the residues in model 1, chain 'L', segment ''.
            >>> list(struct_chainB.residues(model = 1,
                                            chain = "B",
                                            segment = ""))
            [(68, '', 'GLU'), (69, '', 'GLU'), (70, '', 'GLN'),
            (71, '', 'TRP'), (72, '', 'ALA'), (73, '', 'ARG'),
            (74, '', 'GLU'), ..., ]

            # Remove residues (68, '', 'GLU') and (69, '', 'GLU')
            # from chain 'B' in the new structure. We can define
            # the residues using just their sequence number since
            # all residues in the structure have the same insertion
            # code. We can also avoid specifying the chain(s) and
            # segment(s) from which the residues need to be
            # removed since we have only one chain and one segment
            # in the structure. The residues will be removed from
            # all models.
            >>> residues = [68, 69]
            >>> struct_chainB.remove_items(residues = residues)

            # Get the residues now in model 1, chain 'L',
            # segment ''.
            >>> list(struct_chainB.residues(model = 1,
                                            chain = "B",
                                            segment = ""))
            [(70, '', 'GLN'), (71, '', 'TRP'), (72, '', 'ALA'),
            (73, '', 'ARG'), (74, '', 'GLU'), ..., ]

        See also
        --------
        keep_items : Keep only selected items (models, chains, \
        segments, residues, or atoms) in a structure.
        rename_items : Rename items (models, chains, segments, \
        residues, or atoms) in a structure.
        renumber_items : Renumber items (models and residues) \
        in a structure.
        """

        # Set the keyword arguments to pass to the internal method.
        kwargs = \
            {"action" : "remove",
             "elements_type" : "items",
             "selected" : \
                self._get_selected(\
                        models = models,
                        models_attributes = models_attributes,
                        chains = chains,
                        chains_attributes = chains_attributes,
                        segments = segments,
                        segments_attributes = segments_attributes,
                        residues = residues,
                        residues_attributes = residues_attributes,
                        atoms = atoms,
                        atoms_attributes = atoms_attributes)}

        # If the structure needs to be modified in place.
        if in_place:
            
            # Keep only the selected atoms.
            self._modify(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure.
            struct_new = self.get_copy()

            # Keep only the selected atoms in the new structure.
            struct_new._modify(**kwargs)

            # Return the new structure.
            return struct_new


    ###################################################################


    def rename_items(self,
                     level,
                     mapping,
                     models = None,
                     models_attributes = None,
                     chains = None,
                     chains_attributes = None,
                     segments = None,
                     segments_attributes = None,
                     residues = None,
                     residues_attributes = None,
                     atoms = None,
                     atoms_attributes = None,
                     in_place = False):
        """Rename selected items.

        Parameters
        ----------
        level : ``str``, {``"models"``, ``"chains"``, \
            ``"segments"``, ``"residues"``, ``"atoms"``}
            Where the items to be renamed are in the hierarchy.

        mapping : ``dict``
            A dictionary containing the mapping between the old
            names and the new ones.

        models : ``list``, optional
            The models to be considered when renaming.
            
            If no chains, segments, residues, or atoms are
            provided, the mapping will be used to rename 
            the models.

            Otherwise, only elements in the selected
            models will be considered for renaming.

            Each model is defined by its number.

            If not provided, all models will be considered.

        models_attributes : ``dict``, optional
            A dictionary mapping names of models' attributes
            to lists of their accepted values.

            If chains, segments, residues, or atoms are
            provided, only elements in models whose attributes'
            values match the accepted ones will be considered
            for renaming. 

        chains : ``list``, optional
            The chains to be considered when renaming.
            
            If no segments, residues, or atoms are provided,
            the mapping will be used to rename chains.

            Otherwise, only elements in the selected
            chains (and models) will be considered for
            renaming.
            
            Each chain is defined by its chain identifier.

            If not provided, all chains in the selected
            models will be considered.

        chains_attributes : ``dict``, optional
            A dictionary mapping names of chains' attributes
            to lists of their accepted values.

            If segments, residues, or atoms are provided, only
            elements in chains whose attributes' values match the
            accepted ones will be considered for renaming. 

        segments : ``list``, optional
            The segments to be considered when renaming.
            
            If no residues or atoms are provided, the mapping
            will be used to rename segments.

            Otherwise, only elements in the selected
            segments (and chains and models) will be
            considered for renaming.
            
            Each segment is defined by its segment identifier.

            If not provided, all segments in the selected models
            and chains will be considered.

        segments_attributes : ``dict``, optional
            A dictionary mapping names of segments' attributes
            to lists of their accepted values.

            If residues or atoms are provided, only elements
            in segments whose attributes' values match the
            accepted ones will be considered for renaming.

        residues : ``set``, optional
            The residues to be considered when renaming.
            
            If no atoms are provided, the mapping
            will be used to rename residues.

            Otherwise, only elements in the selected
            residues (and models, chains, and segments)
            will be considered for renaming.
            
            Each residues can be defined by its full identifier,
            its sequence number, or its name.

            If not provided, all residues in the selected models,
            chains, and segments will be considered.

        residues_attributes : ``dict``, optional
            A dictionary mapping names of residues' attributes
            to lists of their accepted values.

            If atoms are provided, only atoms in
            residues whose attributes' values match the
            accepted ones will be considered for renaming.

        atoms : ``list``, optional
            The atoms to be renamed.

            If atoms are provided, the mapping will
            be used to rename them.

            Each atom is defined by its serial number.

            If not provided, all atoms in the selected models,
            chains, segments, and residues will be considered.

        atoms_attributes : ``dict``, optional
            A dictionary mapping names of atoms' attributes
            to lists of their accepted values.

            If atoms are provided, only atoms in
            whose attributes' values match the accepted
            ones will be renamed.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new structure will be returned.

        Returns
        -------
        :class:`pdbcraft.structure.Structure` if ``in_place`` is ``False``;
        otherwise, ``None``.

        Examples
        --------
        .. code-block:: python

            # Rename selected items in a structure.

            # Import 'pdbcraft'.
            >>> import pdbcraft

            # Load the example '2M04.cif' - it is a multi-model
            # structure obtained using NMR spectroscopy.
            >>> struct = pdbcraft.load_example("2M04.cif")

            # Get the structure's chains in model 1.
            >>> list(struct.chains(model = 1))
            ['A', 'B']
            
            # Rename chains 'A' and 'B' to 'C' and 'D' in all
            # models.
            >>> struct.rename_items(level = "chains",
                                    mapping = {"A" : "C", "B" : "D"},
                                    in_place = True)

            # Get the structure's chains in model 1.
            >>> list(struct.chains(model = 1))
            ['C', 'D']
        
        See also
        --------
        keep_items : Keep only selected items (models, chains, \
        segments, residues, or atoms) in a structure.
        remove_items : Remove items (models, chains, segments, \
        residues, or atoms) from a structure.
        renumber_items : Renumber items (models and residues) \
        in a structure.
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "rename",
             "level" : level,
             "elements_type" : "items",
             "selected" : \
                self._get_selected(\
                        models = models,
                        models_attributes = models_attributes,
                        chains = chains,
                        chains_attributes = chains_attributes,
                        segments = segments,
                        segments_attributes = segments_attributes,
                        residues = residues,
                        residues_attributes = residues_attributes,
                        atoms = atoms,
                        atoms_attributes = atoms_attributes),
             "mapping" : mapping}


        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected atoms
            self._modify(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected atoms in the new structure
            struct_new._modify(**kwargs)

            # Return the new structure
            return struct_new


    ###################################################################


    def renumber_items(self,
                       level,
                       start,
                       models = None,
                       models_attributes = None,
                       chains = None,
                       chains_attributes = None,
                       segments = None,
                       segments_attributes = None,
                       residues = None,
                       residues_attributes = None,
                       in_place = False):
        """Renumber selected items.

        Parameters
        ----------
        start : ``int``
            The new starting number.

        level : ``str``, {``"models"``, ``"chains"``, \
            ``"segments"``, ``"residues", ``"atoms"``}
            The level of the hierarchy where the items to
            be renumbered are.

        models : ``list``, optional
            The models to be considered when renumbering.
            
            If no chains, segments, or residues are provided,
            the renumbering will be applied to the models.

            Otherwise, only elements in the selected
            models will be considered for renumbering.

            Each model is defined by its number.

            If not provided, all models will be considered.

        models_attributes : ``dict``, optional
            A dictionary mapping names of models' attributes
            to lists of their accepted values.

            If chains, segments, or residues are provided,
            only elements in models whose attributes' values
            match the accepted ones will be considered for
            renumbering.

        chains : ``list``, optional
            The chains to be considered when renumbering.
            
            If no segments or residues are provided,
            the renumbering will be applied to the chains.

            Otherwise, only elements in the selected
            chains (and models) will be considered for
            renumbering.
            
            Each chain is defined by its chain identifier.

            If not provided, all chains in the selected
            models will be considered.

        chains_attributes : ``dict``, optional
            A dictionary mapping names of chains' attributes
            to lists of their accepted values.

            If segments or residues are provided, only
            elements in chains whose attributes' values match the
            accepted ones will be considered for renumbering. 

        segments : ``list``, optional
            The segments to be considered when renumbering.
            
            If no residues are provided, the renumbering
            will be applied to the segments.

            Otherwise, only elements in the selected
            segments (and chains and models) will be
            considered for renumbering.
            
            Each segment is defined by its segment identifier.

            If not provided, all segments in the selected models
            and chains will be considered.

        segments_attributes : ``dict``, optional
            A dictionary mapping names of segments' attributes
            to lists of their accepted values.

            If residues are provided, only elements
            in segments whose attributes' values match the
            accepted ones will be considered for renumbering.

        residues : ``set``, optional
            The residues to be considered when renumbering.
            
            If they are provided, they will be the items
            being renumbered.
            
            Each residues can be defined by its full identifier,
            its sequence number, or its name.

            If not provided, all residues in the selected models,
            chains, and segments will be considered.

        residues_attributes : ``dict``, optional
            A dictionary mapping names of residues' attributes
            to lists of their accepted values.

            Only residues whose attributes' values match the
            accepted ones will be considered for renumbering.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new structure will be returned.

        Returns
        -------
        :class:`pdbcraft.structure.Structure` if ``in_place`` is ``False``;
        otherwise, ``None``.

        Examples
        --------
        .. code-block:: python

            # Renumber selected items in a structure.

            # Import 'pdbcraft'.
            >>> import pdbcraft

            # Load the example '2LPC.pdb' - it is a multi-model
            # structure obtained using NMR spectroscopy.
            >>> struct = pdbcraft.load_example("2LPC.pdb")

            # Get the structure's models.
            >>> list(struct.models())
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
            16, 17, 18, 19, 20, 21]

            # Renumber the structure's models starting
            # from 3.
            >>> struct.renumber_items(level = "models",
                                      start = 3,
                                      in_place = True)

            # Get the structure's models.
            >>> list(struct.models())
            [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
            18, 19, 20, 21, 22, 23]

        See also
        --------
        keep_items : Keep only selected items (models, chains, \
        segments, residues, or atoms) in a structure.
        remove_items : Remove items (models, chains, segments, \
        residues, or atoms) from a structure.
        rename_items : Rename items (models, chains, segments, \
        residues, or atoms) in a structure.
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "renumber",
             "level" : level,
             "elements_type" : "items",
             "start" : start,
             "selected" : \
                self._get_selected(\
                        models = models,
                        models_attributes = models_attributes,
                        chains = chains,
                        chains_attributes = chains_attributes,
                        segments = segments,
                        segments_attributes = segments_attributes,
                        residues = residues,
                        residues_attributes = residues_attributes)}

        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected atoms
            self._modify(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected atoms in the new structure
            struct_new._modify(**kwargs)

            # Return the new structure
            return struct_new


    ###################################################################


    def assign_attribute_value(self,
                               attribute,
                               new_value,
                               level,
                               models = None,
                               models_attributes = None,
                               chains = None,
                               chains_attributes = None,
                               segments = None,
                               segments_attributes = None,
                               residues = None,
                               residues_attributes = None,
                               atoms = None,
                               atoms_attributes = None,
                               in_place = False):
        """Assign a new value to a specific attribute.

        Parameters
        ----------
        attribute : ``str``
            The name of the attribute.

        new_value : ``str``, ``int``, ``float``, ``tuple``
            The new value of the attribute

        level : ``str``, {``"models"``, ``"chains"``, \
            ``"segments"``, ``"residues", ``"atoms"``}
            Where the attribute is in the hierarchy.

        models : ``list``, optional
            If ``level = "models"``, the models whose
            ``attribute``'s value will be re-assigned.

            Otherwise, only elements in the selected
            models will be considered when re-assigning
            the attribute's value.

            Each model is defined by its number.

            If not provided, all models will be considered.

        models_attributes : ``dict``, optional
            A dictionary mapping names of models' attributes
            to lists of their accepted values.

            If ``level = "models"``, this option will
            be ignored.

            Otherwise, only elements in models whose values
            match the accepted ones will be considered when
            re-assigning the attribute's value.

        chains : ``list``, optional
            If ``level = "chains"``, the chains whose
            ``attribute``'s value will be re-assigned.

            If ``level = "models"``, this option will
            be ignored.

            Otherwise, only elements in the selected
            chains and models will be considered when
            re-assigning the attribute's value.
            
            Each chain is defined by its chain identifier.

            If not provided, all chains in the selected
            models will be considered.

        chains_attributes : ``dict``, optional
            A dictionary mapping names of chains' attributes
            to lists of their accepted values.

            If ``level = "chains"`` or above, this option
            will be ignored.

            Otherwise, only elements in chains whose values
            match the accepted ones will be considered when
            re-assigning the attribute's value.

        segments : ``list``, optional
            If ``level = "segments"``, the segments whose
            ``attribute``'s value will be re-assigned.

            If ``level = "chains"`` or above, this option
            will be ignored.

            Otherwise, only elements in the selected
            segments, chains and models will be considered
            when re-assigning the attribute's value.
            
            Each segment is defined by its segment identifier.

            If not provided, all segments in the selected models
            and chains will be considered.

        segments_attributes : ``dict``, optional
            A dictionary mapping names of segments' attributes
            to lists of their accepted values.

            If ``level = "segments"`` or above, this option
            will be ignored.

            Otherwise, only elements in segments whose values
            match the accepted ones will be considered when
            re-assigning the attribute's value.

        residues : ``set``, optional
            If ``level = "residues"``, the residues whose
            ``attribute``'s value will be re-assigned.

            If ``level = "segments"`` or above, this option
            will be ignored.

            Otherwise, only atoms in the selected residues,
            segments, chains and models will be considered
            when re-assigning the attribute's value.
            
            Each residues can be defined by its full identifier,
            its sequence number, or its name.

            If not provided, all residues in the selected models,
            chains, and segments will be considered.

        residues_attributes : ``dict``, optional
            A dictionary mapping names of residues' attributes
            to lists of their accepted values.

            If ``level = "residues"`` or above, this option
            will be ignored.

            Otherwise, only atoms in residues whose values
            match the accepted ones will be considered when
            re-assigning the attribute's value.

        atoms : ``list``, optional
            If ``level = "atoms"``, the atoms whose
            ``attribute``'s value will be re-assigned.

            Otherwise, this option will be ignored.

            Each atom is defined by its serial number.

            If not provided, all atoms in the selected models,
            chains, segments, and residues will be considered.

        in_place : ``bool``, ``False``
            Whether to modify the structure in place. If ``False``,
            a new structure will be returned.

        Returns
        -------
        :class:`pdbcraft.structure.Structure` if ``in_place`` is ``False``;
        otherwise, ``None``.

        Examples
        --------
        .. code-block:: python

            # Assign a new value to an attribute.

            # Import 'pdbcraft'.
            >>> import pdbcraft

            # Load the example '19G.pdb' - it is a single-model
            # structure obtained using X-ray crystallography.
            # Here, the heteroatoms are incorrectly stored in
            # ATOM records.
            >>> struct = pdbcraft.load_example("19G_nohet.pdb")

            # Get whether the residues in chain 'B' are
            # flagged as heteroresidues by using the
            # 'is_het' attribute, which is automatically
            # defined for each residue when parsing the
            # structure from the PDB file (using the HETATM
            # records).
            >>> struct.get_attribute(attribute = "is_het", 
                                     level = "residues",
                                     squeeze = "chains",
                                     chains = ["B"])
            {1: {'B': [False]}}

            # Assign the 'True' value to the 'is_het' attribute
            # to all residues in chain 'B'.
            >>> struct.assign_attribute_value(attribute = "is_het",
                                              new_value = True,
                                              level = "residues",
                                              chains = ["B"],
                                              in_place = True)

            # Get whether the residues in chain 'B' are now
            # flagged as heteroresidues.
            >>> struct.get_attribute(attribute = "is_het", 
                                     level = "residues",
                                     squeeze = "chains",
                                     chains = ["B"])
            {1: {'B': [True]}}

        See also
        --------
        get_attribute : Get the unique values, minimum value, or \
        maximum value of a specific attribute.
        """

        # Set the keyword arguments to pass to the internal method
        kwargs = \
            {"action" : "assign",
             "elements_type" : "attributes",
             "attribute" : attribute,
             "new_value" : new_value,
             "level" : level,
             "selected" : \
                self._get_selected(\
                        models = models,
                        models_attributes = models_attributes,
                        chains = chains,
                        chains_attributes = chains_attributes,
                        segments = segments,
                        segments_attributes = segments_attributes,
                        residues = residues,
                        residues_attributes = residues_attributes,
                        atoms = atoms,
                        atoms_attributes = atoms_attributes)}


        # If the structure needs to be modified in place
        if in_place:
            
            # Keep only the selected atoms
            self._modify(**kwargs)

        # Otherwise
        else:

            # Create a copy of the current structure
            struct_new = self.get_copy()

            # Keep only the selected atoms in the new structure
            struct_new._modify(**kwargs)

            # Return the new structure
            return struct_new


    ###################################################################


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
            a new structure will be returned.

        Returns
        -------
        :class:`pdbcraft.structure.Structure` if ``in_place`` is ``False``;
        otherwise, ``None``.

        Examples
        --------
        .. code-block:: python

            # Assign a new value to an attribute.

            # Import 'pdbcraft'.
            >>> import pdbcraft

            # Load the example '4WI.pdb' - it is a single-model
            # structure obtained using X-ray crystallography.
            >>> struct = pdbcraft.load_example("4WI.pdb")

            # Get the serial number of atom 'H16' in residue '4WI'.
            atom_path = (1, "C", "", (1, "", "4WI"), "H16")
            >>> struct.atom_serial(path = atom_path)
            68

            # Define the sorting of the atoms in residue '4WI'.
            >>> atoms_names = \
                    [# Carbon atoms
                     "C8", "C3", "C1", "C13", "C21", "C6", "C7", "C10",
                     "C5", "C9", "C11", "C12", "C18", "C19", "C20",
                     "C16", "C17", "C23", "C4", "C2", "C14", "C22",
                     "C15",
                     # Nitrogen atoms
                     "N5", "N2", "N3", "N1", "N4",
                     # Oxygen atoms
                     "O1", "O2", "O3", "O4",
                     # Fluorine atoms
                     "F1", "F2", "F3",
                     # Hydrogen atoms
                     "H27", "H28", "H24", "H26", "H25", "H36", "H15",
                     "H17", "H16", "H18", "H20", "H19", "H11", "H9",
                     "H10", "H7", "H6", "H8", "H14", "H12", "H13",
                     "H29", "H30", "H31", "H32", "H1", "H3", "H4",
                     "H5", "H23"]

            # Sort the atoms in the structure.
            >>> struct.sort_atoms_residue_name(\
                    atoms_names = atoms_names,
                    residue_name = "4WI",
                    in_place = True)

            # Get the new serial number of atom 'H16' in residue '4WI'.
            atom_path = (1, "C", "", (1, "", "4WI"), "H16")
            >>> struct.atom_serial(path = atom_path)
            44
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
