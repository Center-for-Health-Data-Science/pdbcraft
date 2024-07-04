#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    pdbparsers.py
#
#    Parsers for PDB files.
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
__doc__ = "Parsers for PDB files."


#######################################################################


# Import from the standard library.
from collections import defaultdict
import logging as log
import os
import re
# Import from 'pdbcraft'.
from .structure import Structure
from . import _defaults


#######################################################################


# Get the module's logger.
logger = log.getLogger(__name__)


#######################################################################


class PDBParser:

    """
    A parser for PDB files.
    """


    ###################################################################


    # Set the character to be used for unnamed chains in the structure
    # built from the PDB file.
    UNKNOWN_CHAIN_ID = ""

    # Set the character to be used for unnamed segments in the
    # structure built from the PDB file.
    UNKNOWN_SEG_ID = ""


    ###################################################################


    def _get_atom_data(self,
                       file,
                       is_extended_pdb,
                       fix_residues_overflow,
                       only_alt_loc):
        """Parse the MODEL/ENDMDL/ATOM/HETATM records of a PDB file.

        Parameters
        ----------
        file : ``str``
            The PDB file to be parsed.

        is_extended_pdb: ``bool``, ``False``
            Whether the PDB file follows the "extended" PDB format,
            where the field usually reserved for a residue's
            insertion code is used as an additional fied to store
            the residue's sequence number.

        fix_residues_overflow : ``bool``, ``True``
            If the PDB file resets the residue numbering mid-chain
            (which can happen when a chain has more than 9999
            residues and the program writing the file could not
            hande it properly), automatically fix the issue by
            continuing the numbering past 9999.

            If ``True``, a warning will be raised and the numbering
            fixed.

            If ``False``, an exception will be raised.

        only_alt_loc : ``str``, ``"A"``
            For atoms having multiple alternate locations, keep only
            those with ``only_alt_loc`` alternate location.

        Returns
        -------
        atom_data : ``dict``
            A dictionary containing the atomic coordinates.

        serials2ids : ``dict``
            A dictionary mapping the atoms' serial numbers to the
            atoms' original identifiers.

        oldids2newids : ``dict``
            A dictionary mapping the atoms' old identifiers to the
            ones used in the structure.
        """

        # Initialize an empty dictionary to store the data parsed from
        # the atomic coordinates' section of the file.
        atom_data = {}

        #-------------------------------------------------------------#

        # Initialize an empty dictionary to store the mapping between
        # the atoms' serials and the atoms' unique identifiers.
        #
        # For each model, an atom's unique identifier contains:
        #
        # - The ID of the chain the atom belongs to.
        # - The ID of the segment the atom belongs to.
        # - The ID of the residue the atom belongs to.
        # - The atom's name.
        serials2ids = {}

        # Initialize an empty dictionary to store the mapping between
        # the original atom identifiers and the new ones (it is
        # especially useful to keep track of what was what when we
        # renumber residues because of 'overflow' problems)
        oldids2newids = {}

        #-------------------------------------------------------------#

        # Initialize an empty set to store the unique models found
        # in the file.
        unique_models = set()

        #-------------------------------------------------------------#

        # Initialize the current model number to 1 (it will be the
        # number assigned to the only model present in the structure,
        # if there is only one model).
        model = 1

        #-------------------------------------------------------------#

        # Initialize a variable to keep track of which model we were
        # on at the previous ATOM/HETATM record to None.
        prev_model = None

        # Initialize a variable to keep track of which chain ID we were
        # on at the previous ATOM/HETATM record to None.
        prev_chain_id = None

        # Initialize a variable to keep track of the sequence number of
        # the residue we were on at the previous ATOM/HETATM record to
        # None.
        prev_res_seq = None

        # Initialize a variable to keep track of the original sequence
        # number of the residue we were on at the previous ATOM/HETATM
        # record to None.
        prev_res_seq_orig = None

        # Initialize a set to store the chains whose residues had
        # sequence numbers that needed fixing.
        renumbered_chain_ids = set()

        #-------------------------------------------------------------#

        # Open the PDB file.
        with open(file, "r") as f:

            # For each line in the file (start enumerating from 1)
            for line_ix, line in enumerate(f, start = 1):

                #-----------------------------------------------------#

                # If the line is empty
                if re.match(r"^\s*$", line):

                    # Skip the line.
                    continue

                #-----------------------------------------------------#

                # Remove the trailing newline character from the line.
                line = line.rstrip("\n")

                #-----------------------------------------------------#

                # Get the record's header.
                header = line[0:6]

                # Log the current record for debug purporses.
                debugstr = \
                    f"Now parsing the '{header}' record on line " \
                    f"{line_ix}..."
                logger.debug(debugstr)

                #-----------------------------------------------------#

                # If the header contains a model's number
                if header == "MODEL ":

                    # Get the current model's number.
                    model = int(line[10:14])

                    # If there are duplicate models
                    if model in unique_models:

                        # Raise an error.
                        errstr = \
                            f"Multiple models numbered '{model}' " \
                            "were found. Models must have unique " \
                            "numbers, and must be numbered " \
                            "sequentially in the PDB file."
                        raise ValueError(errstr)

                    # Log the model's number for debug purporses
                    debugstr = f"Model: {model}."
                    logger.debug(debugstr)

                    # Add it to the list of unique models found in the
                    # file.
                    unique_models.add(model)

                #-----------------------------------------------------#

                # If the line contains an ATOM or HETATM record
                elif header in {"ATOM  ", "HETATM"}:

                    # If the record refers to a heteroatom
                    if header == "HETATM":

                        # Set the heteroresidue flag to True.
                        is_het = True

                    # If the record refers to a normal atom
                    elif header == "ATOM  ":

                        # Set the heteroresidue flag to False.
                        is_het = False

                    #-------------------------------------------------#

                    # Get the atom's name.
                    atom_name = line[12:16].strip()

                    #-------------------------------------------------#

                    # Get the identifier for the atom's alternate
                    # location.
                    alt_loc = line[16].rstrip()

                    # If an alternate location is specified, but it
                    # differs from the one we want
                    if alt_loc and alt_loc != only_alt_loc:

                        # Skip the atom.
                        continue

                    #-------------------------------------------------#

                    # Get the atom's serial.
                    serial = int(line[6:11])

                    #-------------------------------------------------#

                    # Get the atom's residue name.
                    res_name = line[17:20]

                    #-------------------------------------------------#

                    # Get the atom's chain identifier.
                    chain_id = line[21]

                    # If no chain identifier is specified
                    if not chain_id.strip():

                        # Assign an arbitrary chain identifier.
                        chain_id = self.UNKNOWN_CHAIN_ID

                    #-------------------------------------------------#

                    # If we are dealing with the 'extended' PDB format,
                    # where the field usually dedicated to the
                    # insertion code is used to store the extra digit
                    # for the residues' sequence numbers for structures
                    # with more than 9999 residues
                    if is_extended_pdb:

                        # Get the current residue's sequence number
                        # from the fields that usually include the
                        # sequence number plus the field usually
                        # reserved for the residue's insertion code.
                        res_seq = int(line[22:27].split()[0])

                        # The current residue's insertion code is
                        # going to be an empty string.
                        i_code = ""

                    # Otherwise
                    else:

                        # Get the sequence number for the current
                        # residue.
                        res_seq = int(line[22:26].split()[0])

                        # Get the insertion code for the current
                        # residue.
                        i_code = line[26].strip()

                    # Save the residue's original sequence number.
                    res_seq_orig = res_seq

                    # If we are in the same model and chain as before,
                    # but the current residue's sequence number is
                    # lower than the sequence number of the previous
                    # residue
                    if ((prev_model == model) \
                         and (chain_id == prev_chain_id) \
                         and (prev_res_seq is not None \
                              and res_seq < prev_res_seq)):

                        # It means that we are parsing a misformatted
                        # PDB file where numbering re-starts when
                        # we reach 10000 residues.

                        # If we are in a new residue
                        if prev_res_seq_orig != res_seq_orig:

                            # Initialize the error string to an empty
                            # string.
                            errstr = ""

                            # If we had not encountered the overflow
                            # yet for the current chain
                            if chain_id not in renumbered_chain_ids:
                                renumbered_chain_ids.add(chain_id)

                                # Warn the user about the overflow.
                                errstr += \
                                    "The residue numbering for " \
                                    f"chain '{chain_id}' re-starts " \
                                    f"from {res_seq_orig} after " \
                                    f"residue {prev_res_seq_orig}.\n"
                            
                            # If we should fix the overflow
                            if fix_residues_overflow:

                                # Complete the error string and
                                # just raise a warning.
                                errstr += \
                                    f"Residue {res_seq_orig} will " \
                                    f"be residue {prev_res_seq + 1} " \
                                    "in the structure."
                                logger.warning(errstr)

                            # Otherwise
                            else:

                                # Raise an error.
                                raise ValueError(errstr)

                            # Assign a new sequence number to the
                            # current residue.
                            res_seq = prev_res_seq + 1

                        # Otherwise
                        else:

                            # We are in the same residue as the one
                            # before.
                            res_seq = prev_res_seq

                    #-------------------------------------------------#

                    # Try to get the atom's coordinates.
                    try:

                        # Get the x coordinate.
                        x = float(line[30:38])

                        # Get the y coordinate.
                        y = float(line[38:46])

                        # Get the z coordinate.
                        z = float(line[46:54])

                    # If something went wrong
                    except Exception:

                        # Raise an error.
                        errstr = \
                            "Missing or invalid value(s) found in " \
                            f"the coordinates of atom {serial} " \
                            f"(line {line_ix})."
                        raise ValueError(errstr)

                    #-------------------------------------------------#

                    # Try to get the atom's occupancy.
                    try:
                        
                        occupancy = float(line[54:60])

                    # If something went wrong
                    except Exception:

                        # Raise an error
                        errstr = \
                            "Missing or invalid value found " \
                            f"for the occupancy of atom {serial} " \
                            f"(line {line_ix})."
                        raise ValueError(errstr)

                    #-------------------------------------------------#

                    # Try to get the atom's temperature factor
                    # (b-factor).
                    try:
                        
                        temp_factor = float(line[60:66])

                    # If something went wrong
                    except Exception:

                        # Raise an error
                        errstr = \
                            "Missing or invalid value found " \
                            f"for the temperature factor (b-factor)" \
                            f"of atom {serial} (line {line_ix})."
                        raise ValueError(errstr)

                    #-------------------------------------------------#

                    # Get the identifier of the segment the atom
                    # belongs to.
                    seg_id = line[72:76].lstrip().rstrip()

                    # If no segment identifier is specified
                    if not seg_id.strip():

                        # Assign an arbitrary segment identifier.
                        seg_id = self.UNKNOWN_SEG_ID

                    #-------------------------------------------------#

                    # Get the atom's element's name and convert it
                    # into all-uppercase characters.
                    element = line[76:78].strip().upper()

                    #-------------------------------------------------#

                    # Try to get the atom's charge.
                    try:

                        charge = line[78:80].strip()

                    # If something went wrong
                    except Exception:

                        # Raise an error.
                        errstr = \
                            "Missing or invalid value found " \
                            f"for the charge of atom {serial} " \
                            f"(line {line_ix})."
                        raise ValueError(errstr)

                    #-------------------------------------------------#

                    # If the current model's number is not in the
                    # dictionary
                    if model not in atom_data.keys():

                        # Add it.
                        atom_data[model] = \
                            {"_items" : {},
                             "_attributes" : {}}

                    #-------------------------------------------------#

                    # If the current chain's identifier is not in the
                    # dictionary
                    if chain_id not in atom_data[model][\
                        "_items"].keys():

                        # Add it.
                        atom_data[model]["_items"][chain_id] = \
                            {"_items" : {},
                             "_attributes" : {}}

                    #-------------------------------------------------#

                    # If the current segment's identifier is not in the
                    # dictionary
                    if seg_id not in atom_data[model]["_items"][\
                        chain_id]["_items"].keys():

                        # Add it.
                        atom_data[model]["_items"][\
                            chain_id]["_items"][seg_id] = \
                                {"_items" : {},
                                 "_attributes" : {}}

                    #-------------------------------------------------#

                    # If the current residue's identifier is not in the
                    # dictionary
                    if (res_seq, i_code, res_name) not in \
                        atom_data[model]["_items"][chain_id][\
                            "_items"][seg_id]["_items"].keys():

                        # Add it.
                        atom_data[model]["_items"][chain_id][\
                            "_items"][seg_id]["_items"][\
                            (res_seq, i_code, res_name)] = \
                                    {"_items" : {},
                                     "_attributes" : \
                                        {"is_het" : is_het}}

                    #-------------------------------------------------#
                    
                    # Add the data for the current atom to the
                    # dictionary.
                    atom_data[model]["_items"][chain_id]["_items"][\
                        seg_id]["_items"][(res_seq, i_code, \
                            res_name)]["_items"][serial] = \
                                {"_attributes" : \
                                     {"label_atom_id" : atom_name,
                                      "label_alt_id" : alt_loc,
                                      "Cartn_x" : x,
                                      "Cartn_y" : y,
                                      "Cartn_z" : z,
                                      "occupancy" : occupancy,
                                      "B_iso_or_equiv" : temp_factor,
                                      "type_symbol" : element,
                                      "pdbx_formal_charge"  : charge,
                                      "label_entity_id" : ""}}

                    # Log the data about the curren atom for debug
                    # purporses.
                    debugstr = \
                        f"Atom: {serial}. " \
                        f"Atom name: {atom_name}. " \
                        f"Alt. loc.: {alt_loc}. " \
                        f"Res. name: {res_name}. " \
                        f"Chain ID: {chain_id}. " \
                        f"Res. seq.: {res_seq}. " \
                        f"I-Code: {i_code}. " \
                        f"Coordinates: ({x}, {y}, {z}). " \
                        f"Occupancy: {occupancy}. " \
                        f"Temp. factor: {temp_factor}. " \
                        f"Element: {element}. " \
                        f"Charge: {charge}. " \
                        f"Heteroresidue? {is_het}."
                    logger.debug(debugstr)

                    #-------------------------------------------------#

                    # Update the variable storing the number of the
                    # current model.
                    prev_model = model
                    
                    # Update the variable storing the identifier of
                    # the current chain.
                    prev_chain_id = chain_id

                    # Update the variable storing the sequence number
                    # of the current residue.
                    prev_res_seq = res_seq

                    # Update the variable storing the original sequence
                    # number of the current residue.
                    prev_res_seq_orig = res_seq_orig

                    #-------------------------------------------------#

                    # Set the atom's original identifier, with the
                    # original residue number.
                    atom_id_old = \
                        (chain_id,
                         seg_id,
                         (res_seq_orig, i_code, res_name),
                         atom_name)

                    # Set the atom's new identifier, with the residue
                    # number used in the structure.
                    atom_id_new = \
                        (chain_id,
                         seg_id,
                         (res_seq, i_code, res_name),
                         atom_name)

                    #-------------------------------------------------#

                    # If the model has not been encountered yet
                    if model not in serials2ids:
                        
                        # Add an entry for the model in the dictionary
                        # mapping the atoms' serials to their original
                        # idenfitiers.
                        serials2ids[model] = \
                            {serial : atom_id_old}

                    # Otherwise
                    else:

                        # Update the model's entry in the dictionary
                        # mapping the atoms' serials to their original
                        # idenfitiers.
                        serials2ids[model][serial] = atom_id_old

                    #-------------------------------------------------#

                    # If the model has not been encountered yet
                    if model not in oldids2newids:
                        
                        # Add an entry for the model in the dictionary
                        # mapping the atoms' original idenfitiers to
                        # the ones used in the structure.
                        oldids2newids[model] = \
                            {atom_id_old : atom_id_new}

                    # Otherwise
                    else:

                        # Update the model's entry in the dictionary
                        # mapping the atoms' original idenfitiers to
                        # the ones used in the structure.
                        oldids2newids[model][atom_id_old] = \
                            atom_id_new

        #-------------------------------------------------------------#

        # If no atomic coordinates were found
        if not atom_data:

            # Raise an error.
            errstr = \
                "Empty PDB file: no atomic coordinates found " \
                f"in '{file}'."
            raise ValueError(errstr)

        #-------------------------------------------------------------#

        # Return the parsed data.
        return atom_data, serials2ids, oldids2newids


    def _get_conect_data(self,
                         file,
                         serials2ids):
        """Parse the CONECT records of a PDB file.

        Parameters
        ----------
        file : ``str``
            The PDB file to be parsed.

        serials2ids : ``dict``
            A dictionary mapping the atoms' serial numbers to
            the atoms' original identifiers.
        
        Returns
        -------
        conect_data : ``dict``
            A dictionary containing the data parsed from the
            CONECT records.
        """

        # Initialize an empty dictionary to store the data parsed
        # from the CONECT records.
        conect_data = defaultdict(lambda: defaultdict(dict))

        #-------------------------------------------------------------#

        # Open the PDB file.
        with open(file, "r") as f:

            # Log the fact that we are getting the connectivity data.
            debugstr = \
                "Now getting the connectivity data..."
            logger.debug(debugstr)

            # For each line in the file (start enumerating from 1)
            for line_ix, line in enumerate(f, start = 1):

                #-----------------------------------------------------#

                # If the line is empty
                if re.match(r"^\s*$", line):

                    # Skip the line.
                    continue

                #-----------------------------------------------------#

                # Remove the trailing newline character from the line.
                line = line.rstrip("\n")

                #-----------------------------------------------------#

                # Get the record's header.
                header = line[0:6]

                # Log the current record for debug purporses.
                debugstr = \
                    f"Now parsing the '{header}' record on line " \
                    f"{line_ix}..."
                logger.debug(debugstr)

                #-----------------------------------------------------#

                # If the header indicates a CONECT record
                if header == "CONECT":

                    # Get the record.
                    record = line.lstrip("CONECT")

                    # Split the record into chunks containing the
                    # serial numbers of the bonded atoms.
                    atoms = \
                        [int(record[i:i+5]) \
                         for i in range(0, len(record), 5) \
                         if record[i:i+5].strip() != ""]

                    # Get the atom the record refers to and the list
                    # of atoms bonded to it.
                    atom_1_serial, atoms_bonded_serials = \
                        atoms[0], atoms[1:]

                    # Get the first atom's identifier.
                    atom_1 = serials2ids[1][atom_1_serial]

                    # Initialize a new CONECT record.
                    new_record = {}
                    
                    # For each bonded atom
                    for atom_2_serial in atoms_bonded_serials:

                        # Get the second atom's identifier.
                        atom_2 = serials2ids[1][atom_2_serial]

                        # If the atom is already in the record
                        # (= multiple bonds are specified by
                        # repeating the bonded atoms as many times
                        # as the order of the bond)
                        if atom_2 in new_record:

                            # Take the bond order.
                            bond_order = \
                                new_record[atom_2]["pdbx_value_order"]
                            
                            # If the bond was a single bond
                            if bond_order == "sing":

                                # It is now a double bond.
                                new_bond_order = "doub"

                            # If the bond was a double bond
                            elif bond_order == "doub":

                                # If is now a triple bond.
                                new_bond_order = "trip"

                            # If the bond was a triple bond
                            elif bond_order == "trip":

                                # It is now a quadruple bond.
                                new_bond_order = "quad"

                            # Update the bond order.
                            new_record[atom_2][\
                                "pdbx_value_order"] = new_bond_order

                        # If the bonded atom is not in the record yet
                        else:

                            # Create a record for the bonded atom.
                            new_record[atom_2] = \
                                {"conn_type_id" : "covale",
                                 "pdbx_value_order" : "sing"}

                    #-------------------------------------------------#

                    # Add the new record to the dictionary storing
                    # the connectivity data.
                    conect_data[1][atom_1] = new_record

        #-------------------------------------------------------------#

        # If there are connectivity data
        if conect_data:

            # Repeat the connectivity data for as many models as there
            # are in the structure and return the connectivity data.
            return \
                {mod : \
                    {atom1 : \
                        {atom2 : data for atom2, data \
                         in conect_data[mod][atom1].items()} \
                     for atom1 in conect_data[mod]} \
                 for mod in serials2ids}

        # Otherwise
        else:

            # Return an empty dictionary.
            return {}


    def _get_ssbond_data(self,
                         file):
        """Parse the SSBOND records in a PDB file.

        Parameters
        ----------
        file : ``str``
            The PDB file to be parsed.

        Returns
        -------
        ssbond_data : ``dict``
            A dictionary containing the data parsed from the SSBOND
            records.
        """

        # Initialize a dictionary to store the disulfide bonds parsed
        # from the SSBOND records.
        ssbond_data = {}

        #-------------------------------------------------------------#

        # Open the PDB file.
        with open(file, "r") as f:

            # Log the fact that we are getting the disulfide bonds'
            # data.
            debugstr = \
                "Now getting the disulfide bonds' data from " \
                "the SSBOND records..."
            logger.debug(debugstr)

            # For each line in the file (start enumerating from 1)
            for line_ix, line in enumerate(f, start = 1):

                #-----------------------------------------------------#

                # If the line is empty
                if re.match(r"^\s*$", line):

                    # Skip the line.
                    continue

                #-----------------------------------------------------#

                # Remove the trailing newline character from the line.
                line = line.rstrip("\n")

                #-----------------------------------------------------#

                # Get the record's header.
                header = line[0:6]

                # Log the current record for debug purporses.
                debugstr = \
                    f"Now parsing the '{header}' record on line " \
                    f"{line_ix}..."
                logger.debug(debugstr)

                #-----------------------------------------------------#

                # If the header indicates the section storing disulfide
                # bonds
                if header == "SSBOND":

                    # Get the disulfide bond's number.
                    ssbond_num = int(line[7:10])

                    # Get the name of the first cysteine.
                    cys1_name = line[11:14].strip()

                    # Get the chain identifier of the first cysteine.
                    cys1_chain_id = line[15]

                    # Get the sequence number of the first cysteine.
                    cys1_seq = int(line[17:21])

                    # Get the insertion code of the first cysteine.
                    cys1_i_code = line[21].strip()

                    # Get the name of the second cysteine.
                    cys2_name = line[25:28].strip()

                    # Get the chain identifier of the second cysteine.
                    cys2_chain_id = line[29]

                    # Get the sequence number of the second cysteine.
                    cys2_seq = int(line[31:35])

                    # Get the insertion code of the second cysteine.
                    cys2_i_code = line[35].strip()

                    # Get the symmetry operator of the first cysteine.
                    cys1_symop = line[59:65].strip()

                    # Get the symmetry operator of the second cysteine.
                    cys2_symop = line[66:72].strip()

                    # Get the disulfide bond's length.
                    bond_length = float(line[73:78])

                    #-------------------------------------------------#

                    # Set the identifier for the first SG atom's (to
                    # add bond later).
                    atom1_id = \
                        (cys1_chain_id, "",
                         (cys1_seq, cys1_i_code, cys1_name), 
                         "SG")

                    # Set the identifier for the second SG atom's (to
                    # add bond later).
                    atom2_id = \
                        (cys2_chain_id, "",
                         (cys2_seq, cys2_i_code, cys2_name), 
                         "SG")

                    # Set the disulfide bond's identifier (to add the
                    # bond later).
                    bond_id = (atom1_id, atom2_id)

                    #-------------------------------------------------#

                    # Add the bond to the dictionary.
                    ssbond_data[bond_id] = \
                            {"conn_type_id" : "disulf",
                             "ptnr1_symmetry" : cys1_symop,
                             "ptnr2_symmetry" : cys2_symop,
                             "pdbx_dist_value" : bond_length,
                             "pdbx_value_order" : "sing"}

        #-------------------------------------------------------------#

        # Return the dictionary storing the data parsed from the
        # SSBOND records.
        return ssbond_data


    def _get_link_data(self,
                       file,
                       only_alt_loc):
        """Parse the LINK records in a PDB file.

        Parameters
        ----------
        file : ``str``
            The PDB file to be parsed.

        only_alt_loc : ``str``, ``"A"``
            For atoms having multiple alternate locations, keep only
            those with ``only_alt_loc`` alternate location.

        Returns
        -------
        link_data : ``dict``
            A dictionary containing the data parsed from the LINK
            records.
        """

        # Initialize a dictionary to store the data parsed from the
        # LINK records.
        link_data = {}

        #-------------------------------------------------------------#

        # Open the PDB file.
        with open(file, "r") as f:

            # Log the fact that we are getting the LINK data.
            debugstr = \
                "Now getting the 'link' bonds' data from the LINK " \
                "records..."
            logger.debug(debugstr)

            # For each line in the file (start enumerating from 1)
            for line_ix, line in enumerate(f, start = 1):

                #-----------------------------------------------------#

                # If the line is empty
                if re.match(r"^\s*$", line):

                    # Skip the line.
                    continue

                #-----------------------------------------------------#

                # Remove the trailing newline character from the line.
                line = line.rstrip("\n")

                #-----------------------------------------------------#

                # Get the record's header.
                header = line[0:6].strip()

                # Log the current record for debug purporses.
                debugstr = \
                    f"Now parsing the '{header}' record on line " \
                    f"{line_ix}..."
                logger.debug(debugstr)

                #-----------------------------------------------------#

                # If the header indicates the section storing LINK
                # records.
                if header == "LINK":

                    # Get the name of the atom in the first residue.
                    res1_atom_name = line[12:16].strip()

                    # Get the alternate location of the atom in the
                    # first residue.
                    res1_atom_alt_loc = line[16].strip()

                    # Get the name of the first residue.
                    res1_name = line[17:20].strip()

                    # Get the chain identifier of the first residue.
                    res1_chain_id = line[21]

                    # Get the sequence number of the first residue.
                    res1_seq = int(line[22:26])

                    # Get the insertion code of the first residue.
                    res1_i_code = line[26].strip()

                    # Get the name of the atom in the second residue.
                    res2_atom_name = line[42:46].strip()

                    # Get the alternate location of the atom in the
                    # second residue.
                    res2_atom_alt_loc = line[46].strip()

                    # Get the name of the second residue.
                    res2_name = line[47:50].strip()

                    # Get the chain identifier of the second residue.
                    res2_chain_id = line[51]

                    # Get the sequence number of the second residue.
                    res2_seq = int(line[52:56])

                    # Get the insertion code of the second residue.
                    res2_i_code = line[56].strip()

                    # Get the symmetry operator of the first residue.
                    res1_symop = line[59:65].strip()

                    # Get the symmetry operator of the second residue.
                    res2_symop = line[66:72].strip()

                    # Get the bond's length.
                    bond_length = float(line[73:78])

                    #-------------------------------------------------#

                    # If the first atom has an alternate location, but
                    # it is not the one we want
                    if res1_atom_alt_loc \
                    and res1_atom_alt_loc != only_alt_loc:

                        # Go to the next record.
                        continue

                    # If the second atom has an alternate location, but
                    # it is not the one we want
                    if res2_atom_alt_loc \
                    and res2_atom_alt_loc != only_alt_loc:

                        # Go to the next record.
                        continue

                    #-------------------------------------------------#

                    # Set the identifier of the first atom (to add the
                    # bond later).
                    atom1_id = \
                        (res1_chain_id, "",
                         (res1_seq, res1_i_code, res1_name), 
                         res1_atom_name)

                    # Set the identifier of the second atom (to add
                    # the bond later).
                    atom2_id = \
                        (res2_chain_id, "",
                         (res2_seq, res2_i_code, res2_name), 
                         res2_atom_name)

                    # Set the bond's identifier (to add the bond
                    # later).
                    bond_id = (atom1_id, atom2_id)

                    #-------------------------------------------------#

                    # Add the bond to the dictionary.
                    link_data[bond_id] = \
                        {"conn_type_id" : "covale",
                         "ptnr1_symmetry" : res1_symop,
                         "ptnr2_symmetry" : res2_symop,
                         "pdbx_dist_value" : bond_length,
                         "pdbx_value_order" : "sing"}

        #-------------------------------------------------------------#

        # Return the dictionary storing the data parsed from the
        # LINK records.
        return link_data


    def _get_all_conect(self,
                        conect_data,
                        ssbond_data,
                        link_data,
                        oldids2newids):
        """Get all the connectivity data.

        Parameters
        ----------
        conect_data : ``dict``
            A dictionary containing the data parsed from the
            CONECT records.

        ssbond_data : ``dict``
            A dictionary containing the data parsed from the
            SSBOND records.

        link_data : ``dict``
            A dictionary containing the data parsed from the
            LINK records.

        oldids2newids : ``dict``
            A dictionary mapping the atoms' original identifiers to
            the ones used in the structure.

        Returns
        -------
        all_conect : ``dict``
            A dictionary containing all the connectivity data.
        """

        # Initialize an empty dictionary to store the final
        # connectivity data.
        all_conect = {}

        #-------------------------------------------------------------#

        # For each model
        for mod in conect_data:

            # Add it to the dictionary.
            all_conect[mod] = {}

            # Initialize a dictionary to store the current unique
            # index for disulfide and covalent bonds.
            bond_index = {"disulf" : 1, "covale" : 1}

            #---------------------------------------------------------#

            # For each disulfide bond in the disulfide bonds' data
            for disulf in ssbond_data:

                # Get the original identifiers of the two atoms.
                old_disulf_atom_1, old_disulf_atom_2 = disulf

                # Get the new identifier of the first atom.
                disulf_atom_1 = oldids2newids[mod][old_disulf_atom_1]

                # Get the new identifier for the second atom.
                disulf_atom_2 = oldids2newids[mod][old_disulf_atom_2]

                # Get the bond's data.
                disulf_bond_data = ssbond_data[disulf]

                # Update the bond's data with the bond's identifier.
                disulf_bond_data["id"] = \
                    f"disulf{bond_index['disulf']}"

                #-----------------------------------------------------#

                # If the first atom is already in the connectivity
                # data
                if disulf_atom_1 in all_conect[mod]:

                    # Update the corresponding record.
                    all_conect[mod][disulf_atom_1][disulf_atom_2] = \
                        disulf_bond_data

                # If the first atom is not in the connectivity data
                else:

                    # Create a new record for the atom.
                    all_conect[mod][disulf_atom_1] = \
                        {disulf_atom_2 : disulf_bond_data}

                #-----------------------------------------------------#

                # If the second atom is already in the connectivity
                # data
                if disulf_atom_2 in all_conect[mod]:

                    # Update the corresponding record.
                    all_conect[mod][disulf_atom_2][disulf_atom_1] = \
                        disulf_bond_data

                # If the second atom is not in the connectivity data
                else:

                    # Create a new record for the atom.
                    all_conect[mod][disulf_atom_2] = \
                        {disulf_atom_1 : disulf_bond_data}

                #-----------------------------------------------------#

                # Update the bond's index.
                bond_index["disulf"] += 1

            #---------------------------------------------------------#

            # For each LINK bond in the LINK data
            for link in link_data:

                # Get the original identifiers of the two atoms.
                old_link_atom_1, old_link_atom_2 = link

                # Get the new identifier of the first atom.
                link_atom_1 = oldids2newids[mod][old_link_atom_1]

                # Get the new identifier of the second atom.
                link_atom_2 = oldids2newids[mod][old_link_atom_2]

                # Get the bond's data.
                link_bond_data = link_data[link]

                # Update the bond's data with the bond's identifier.
                link_bond_data["id"] = \
                    f"covale{bond_index['covale']}"

                #-----------------------------------------------------#

                # If the first atom is already in the connectivity data
                if link_atom_1 in all_conect[mod]:

                    # If the second atom is already among the bonded
                    # atoms of the first atom
                    if link_atom_2 in all_conect[mod][link_atom_1]:
                        
                        # Ignore the bond (= we already saved it).
                        continue

                    # Otherwise
                    else:

                        # Update the corresponding record.
                        all_conect[mod][link_atom_1][link_atom_2] = \
                            link_bond_data

                # If the first atom is not in the connectivity data
                else:

                    # Create a new record for the atom.
                    all_conect[mod][link_atom_1] = \
                        {link_atom_2 : link_bond_data}

                #-----------------------------------------------------#

                # If the second atom is already in the connectivity
                # data
                if link_atom_2 in all_conect[mod]:

                    # If the first atom is already among the bonded
                    # atoms of the second atom
                    if link_atom_1 in all_conect[mod][link_atom_2]:
                        
                        # Ignore the bond (= we already saved it).
                        continue

                    # Otherwise
                    else:

                        # Update the corresponding record.
                        all_conect[mod][link_atom_2][link_atom_1] = \
                            link_bond_data

                # If the first atom is not in the connectivity data
                else:

                    # Create a new record for the atom.
                    all_conect[mod][link_atom_2] = \
                        {link_atom_1 : link_bond_data}

                #-----------------------------------------------------#

                # Update the bond's index.
                bond_index["covale"] += 1
            
            #---------------------------------------------------------#

            # For each 'first' atom in the dictionary of covalent bonds
            for old_coval_atom_1 in conect_data[mod]:

                # For each 'second' atom bonded to the first one
                for old_coval_atom_2 in \
                    conect_data[mod][old_coval_atom_1]:

                    # Get the new identifier of the first atom.
                    coval_atom_1 = oldids2newids[mod][old_coval_atom_1]

                    # Get the new identifier of the second atom.
                    coval_atom_2 = oldids2newids[mod][old_coval_atom_2]

                    # Get the bond's data.
                    coval_bond_data = \
                        conect_data[mod][old_coval_atom_1][\
                            old_coval_atom_2]

                    # Update the bond's data with the bond's
                    # identifier.
                    coval_bond_data["id"] = \
                        f"covale{bond_index['covale']}"

                    #-------------------------------------------------#

                    # If the first atom is already in the connectivity
                    # data
                    if coval_atom_1 in all_conect[mod]:

                        # If the second atom is already among the
                        # bonded atoms of the first atom
                        if coval_atom_2 in all_conect[mod][\
                            coval_atom_1]:
                            
                            # Ignore the bond (= we already saved it).
                            continue

                        # Otherwise
                        else:

                            # Add the bond to the data.
                            all_conect[mod][coval_atom_1][\
                                coval_atom_2] = coval_bond_data
                    
                    # If the first atom is not in the connectivity data
                    else:

                        # Save the current bond for the first atom.
                        all_conect[mod][coval_atom_1] = \
                           {coval_atom_2 : coval_bond_data}

                    #-------------------------------------------------#

                    # If the first second is already in the
                    # connectivity data
                    if coval_atom_2 in all_conect[mod]:

                        # If the first atom is already among the bonded
                        # atoms of the second atom
                        if coval_atom_1 in all_conect[mod][\
                            coval_atom_2]:
                            
                            # Ignore the bond (= we already saved it).
                            continue

                        # Otherwise
                        else:

                            # Add the bond to the data.
                            all_conect[mod][coval_atom_2][\
                                coval_atom_1] = coval_bond_data
                    
                    # If the second atom is not in the connectivity
                    # data
                    else:

                        # Save the current bond for the second atom.
                        all_conect[mod][coval_atom_2] = \
                           {coval_atom_1 : coval_bond_data}

                    #-------------------------------------------------#

                    # Update the bond's index.
                    bond_index["covale"] += 1

        #-------------------------------------------------------------#

        # Return the connectivity data.
        return all_conect


    ###################################################################


    def parse(self,
              file,
              only_alt_loc = "A",
              parse_conect_records = True,
              parse_ssbond_records = True,
              parse_link_records = True,
              is_extended_pdb = False,
              fix_residues_overflow = True):
        """Parse a PDB file.

        Parameters
        ----------
        file : ``str``
            The PDB file to be parsed.

        only_alt_loc : ``str``, ``"A"``
            For atoms having multiple alternate locations, keep only
            those with ``only_alt_loc`` alternate location.

        parse_conect_data : ``bool``, ``True``
            Whether to parse CONECT records.

        parse_ssbond_data : ``bool``, ``True``
            Whether to parse SSBOND records.

        parse_link_data : ``bool``, ``True``
            Whether to parse LINK records.

        is_extended_pdb: ``bool``, ``False``
            Whether the PDB file follows the "extended" PDB format,
            where the field usually reserved for a residue's
            insertion code is used as an additional fied to store
            the residue's sequence number.

        fix_residues_overflow : ``bool``, ``True``
            If the PDB file resets the residue numbering mid-chain
            (which can happen when a chain has more than 9999
            residues and the program writing the file could not
            hande it properly), automatically fix the issue by
            continuing the numbering past 9999.

            If ``True``, a warning will be raised and the numbering
            fixed.

            If ``False``, an exception will be raised.

        Returns
        -------
        struct : ``pdbcraft.structure.Structure``
            The parsed structure.
        """

        # Get the atomic coordinates.
        atom_data, serials2ids, oldids2newids = \
            self._get_atom_data(\
                file = file,
                is_extended_pdb = is_extended_pdb,
                fix_residues_overflow = fix_residues_overflow,
                only_alt_loc = only_alt_loc)

        # Initialize the connectivity data to empty dictionaries.
        conect_data = {}
        ssbond_data = {}
        link_data = {}

        #-------------------------------------------------------------#

        # If we need to parse the CONECT records
        if parse_conect_records:
            
            # Get the connectivity data from the CONECT records.
            conect_data = \
                self._get_conect_data(file = file,
                                      serials2ids = serials2ids)

        #-------------------------------------------------------------#

        # If we need to parse the SSBOND records
        if parse_ssbond_records:
            
            # Get the disulfide bonds' data from the SSBOND records.
            ssbond_data = \
                self._get_ssbond_data(file = file)

        #-------------------------------------------------------------#

        # If we need to parse the LINK records
        if parse_link_records:

            # Get the LINK data from the LINK records.
            link_data = \
                self._get_link_data(file = file,
                                    only_alt_loc = only_alt_loc)

        #-------------------------------------------------------------#

        # Get all connectivity data.
        all_conect = \
            self._get_all_conect(conect_data = conect_data,
                                 ssbond_data = ssbond_data,
                                 link_data = link_data,
                                 oldids2newids = oldids2newids)

        #-------------------------------------------------------------#

        # Get the name of the file (without the path leading to it,
        # if any).
        file_name, _ = os.path.splitext(os.path.basename(file))

        #-------------------------------------------------------------#

        # Create the structure.
        struct = \
            Structure(atom_data = atom_data,
                      conect_data = all_conect,
                      name = file_name)

        # Renumber the atoms since atoms' serial numbers may not be
        # continuous if there were atoms with alternate locations.
        struct._update_atom_and_conect()

        # Inform the user about the structure's creation.
        infostr = \
            f"The structure '{file_name}' was created from '{file}'."
        logger.info(infostr)

        #-------------------------------------------------------------#

        # Return the structure.
        return struct
