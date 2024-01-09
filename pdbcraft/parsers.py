#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    parsers.py
#
#    Parsers for files containing atomic structures.
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
import abc
import logging as log
import os
import re
# pdbcraft
from .structure import Structure


# Get the module's logger
logger = log.getLogger(__name__)


class StructureParser(metaclass = abc.ABCMeta):

    """
    Abstract class implementing a generic parser for structures'
    files.
    """

    @classmethod
    def __subclasshook__(cls,
                         subclass):
        """Every class with a 'parse' method will be a concrete
        implementation of ``StructureParser``.
        """

        return (hasattr(subclass, "parse") \
                and callable(subclass.parse))


class PDBParser:

    """
    A parser for PDB files.
    """


    #------------------------ Public methods -------------------------#


    def parse(self,
              file):
        """Parse a PDB file.

        Parameters
        ----------
        file : ``str``
            The PDB file to be parsed.

        Returns
        -------
        struct : ``pdbcraft.structure.Structure``
            The parsed structure.
        """


        #---------------------- Initialization -----------------------#


        # Initialize an empty dictionary to store the data parsed
        # from the atomic coordinates' section of the file
        atom_data = {}

        # Initialize an empty dictionary to store the data parsed
        # from the CONECT section of the file
        conect_data = {}

        # Initialize an empty set to store the unique models found
        # in the file
        unique_models = set()

        # Initialize the current model number to 1 (it will be the
        # number assigned to the only model present in the structure,
        # if there is only one model)
        model = 1

        # Initialize the ID to be used for unnamed chains
        unknown_chain_id = ""

        # Initialize the ID to be used for unknown segments
        unknown_seg_id = ""


        #-------------------------- Parsing --------------------------#


        # Open the PDB file
        with open(file, "r") as f:

            # For each line in the file (start enumerating from 1)
            for line_ix, line in enumerate(f, start = 1):

                # If the line is empty:
                if re.match(r"^\s*$", line):

                    # Skip the line
                    continue

                # Remove the trailing newline character from the line
                line = line.rstrip("\n")

                # Get the record's header
                header = line[0:6]

                # Log the current record for debug purporses
                debugstr = \
                    f"Now parsing '{header}' record on line " \
                    f"{line_ix}..."
                logger.debug(debugstr)

                #------------------- CONECT record -------------------#

                # If the header indicates a CONECT record
                if header == "CONECT":

                    # Get the record
                    record = line.lstrip("CONECT")

                    # Split the record into chunks containing the
                    # serial numbers of the bonded atoms
                    atoms = \
                        [int(record[i:i+5]) \
                         for i in range(0, len(record), 5) \
                         if record[i:i+5].strip() != ""]

                    # Get the atom the record refers to and the
                    # atoms bonded to the first atom
                    atom_record, atoms_bonded = atoms[0], atoms[1:]

                    # If there is no CONECT record for the
                    # current model yet
                    if not model in conect_data:

                        # Create the section for the current model
                        # and add the record
                        conect_data[model] = \
                            {atom_record : atoms_bonded}

                    # Otherwise
                    else:

                        # Just add the record
                        conect_data[model][atom_record] = atoms_bonded

                    # Log the atoms for debug purporses
                    debugstr = \
                        f"First atom: {atom_record}. " \
                        "Bonded atoms: " \
                        f"{', '.join([f'{a}' for a in atoms_bonded])}."
                    logger.debug(debugstr)

                #------------------- MODEL record --------------------#

                # If the header contains a model number
                elif header == "MODEL ":

                    # Get the current model number
                    model = int(line[10:14])

                    # If there are duplicate models
                    if model in unique_models:

                        # Raise an error
                        errstr = \
                            "Multiple models numbered " \
                            f"'{model}' were found. Models " \
                            "must have unique numbers, and " \
                            "be numbered sequentially in the " \
                            "PDB file."
                        raise ValueError(errstr)

                    # Log the model for debug purporses
                    debugstr = f"Model: {model}."
                    logger.debug(debugstr)

                    # Add it to the list of unique models
                    unique_models.add(model)

                #---------------- ATOM/HETATM record -----------------#

                # If the line contains an ATOM or HETATM record
                elif header in {"ATOM  ", "HETATM"}:

                    #-------------------- Header ---------------------#

                    # If the record refers to a heteroatom
                    if header == "HETATM":

                        # Set the heteroatom flag to True
                        is_hetatm = True

                    # If the record refers to a normal atom
                    elif header == "ATOM  ":

                        # Set the heteroatom flag to False
                        is_hetatm = False

                    #-------------------- Serial ---------------------#
                    
                    # Get the atom's serial number
                    serial = int(line[6:11])

                    #------------------- Atom name -------------------#

                    # Get the atom's name
                    atom_name = line[12:16].strip()

                    #-------------- Alternate location ---------------#

                    # Get the atom's alternate location identifier
                    alt_loc = line[16]

                    #----------------- Residue name ------------------#

                    # Get the atom's residue name
                    res_name = line[17:20]

                    #------------------- Chain ID --------------------#

                    # Get the atom's chain identifier
                    chain_id = line[21]

                    # If no chain identifier is specified
                    if not chain_id.strip():

                        # Assign an arbitrary chain identifier
                        # based on the fact that the chain identifier
                        # is unknown
                        chain_id = unknown_chain_id

                    #---------------- Residue number -----------------#

                    # Get the atom's residue's position in the sequence
                    res_seq = int(line[22:26].split()[0])

                    #--------------------- iCODE ---------------------#

                    # Get the atom's residue's insertion code
                    i_code = line[26]

                    #----------------- Atom position -----------------#

                    # Try to get the atom's position
                    try:

                        # Get the x coordinate of the atom's position
                        x = float(line[30:38])

                        # Get the y coordinate of the atom's position
                        y = float(line[38:46])

                        # Get the z coordinate of the atom's position
                        z = float(line[46:54])

                    # If something went wrong
                    except Exception:

                        # Raise an error
                        errstr = \
                            f"Missing or invalid value(s) found in " \
                            f"atom {serial} coordinates " \
                            f"(line {line_ix})."
                        raise ValueError(errstr)

                    #---------------- Atom occupancy -----------------#

                    # Try to get the atom's occupancy
                    try:
                        
                        occupancy = float(line[54:60])

                    # If something went wrong
                    except Exception:

                        # Raise an error
                        errstr = \
                            f"Missing or invalid value found " \
                            f"for atom {serial}'s occupancy " \
                            f"(line {line_ix})."
                        raise ValueError(errstr)

                    #------------ Atom temperature factor ------------#

                    # Try to get the atom's temperature factor
                    # (b-factor)
                    try:
                        
                        temp_factor = float(line[60:66])

                    # If something went wrong
                    except Exception:

                        # Raise an error
                        errstr = \
                            f"Missing or invalid value found " \
                            f"for atom {serial}'s temperature " \
                            f"factor (b-factor) (line {line_ix})."
                        raise ValueError(errstr)

                    #------------------ Segment ID -------------------#

                    # Get the identifier of the segment the atom
                    # belongs to
                    seg_id = line[72:76]

                    # If no segment identifier is specified
                    if not seg_id.strip():

                        # Assign an arbitrary segment identifier
                        # based on the fact that the segment
                        # identifier is unknown
                        seg_id = unknown_seg_id

                    #----------------- Atom element ------------------# 

                    # Get the atom's element's name and convert it
                    # into all uppercase characters
                    element = line[76:78].strip(" ").upper()

                    #------------------ Atom charge ------------------#

                    # Try to get the atom's charge
                    try:

                        charge = line[78:80].strip()

                    # If something went wrong
                    except Exception:

                        # Raise an error
                        errstr = \
                            f"Missing or invalid value found " \
                            f"for atom {serial}'s charge " \
                            f"(line {line_ix})."
                        raise ValueError(errstr)

                    #-------------------- Update ---------------------#

                    # If the model is not yet in the dictionary
                    if model not in atom_data.keys():

                        # Add it
                        atom_data[model] = {}

                    # If the chain identifier is not yet in the
                    # dictionary
                    if chain_id not in atom_data[model].keys():

                        # Add it
                        atom_data[model][chain_id] = {}

                    # If the segment identifier is not yet in the
                    # dictionary
                    if seg_id not in atom_data[model][chain_id].keys():

                        # Add it
                        atom_data[model][chain_id][seg_id] = {}

                    # If the residue number is not yet in the
                    # dictionary
                    if (res_seq, i_code) not in \
                        atom_data[model][chain_id][seg_id].keys():

                        # Add it
                        atom_data[model][chain_id][seg_id][\
                            (res_seq, i_code)] = \
                                {"res_name" : res_name, "atoms" : {}}
                    
                    # Add the data for the current atom
                    atom_data[model][chain_id][seg_id][\
                        (res_seq, i_code)][\
                        "atoms"][serial] = \
                            {"atom_name" : atom_name,
                             "alt_loc" : alt_loc,
                             "x" : x,
                             "y" : y,
                             "z" : z,
                             "occupancy" : occupancy,
                             "temp_factor" : temp_factor,
                             "element" : element,
                             "charge"  : charge,
                             "is_hetatm" : is_hetatm}

                    # Log the atoms' data for debug purporses
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
                        f"Heteroatom? {is_hetatm}."
                    logger.debug(debugstr)

        # Set the atom data to be used to create the structure
        # to None if the dictionary is empty
        atom_data_final = atom_data if atom_data else None

        # Set the CONECT data to be used to create the structure
        # to None if the dictionary is empty
        conect_data_final = conect_data if conect_data else None

        # Get the name of the file (without the path leading
        # to it, if any)
        file_name, _ = os.path.splitext(os.path.basename(file))

        # Create the structure
        struct = \
            Structure(atom_data = atom_data_final,
                      conect_data = conect_data_final,
                      name = file_name)

        # Inform the user about the structure's creation
        infostr = \
            f"Structure '{file_name}' created " \
            f"from '{file}'."
        logger.info(infostr)

        # Return the structure
        return struct
