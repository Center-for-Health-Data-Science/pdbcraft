#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    pdbwriters.py
#
#    Writers for PDB files.
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
__doc__ = "Writers for PDB files."


#######################################################################


# Import from the standard library.
from collections import defaultdict
import logging as log
# Import from 'pdbcraft'.
from . import _defaults


#######################################################################


# Get the module's logger.
logger = log.getLogger(__name__)


#######################################################################


class PDBWriter:

    """
    A writer for PDB files.
    """

    ###################################################################


    def _write_atom_records(self,
                            struct,
                            file_handle,
                            write_models_records):
        """Write ATOM/HETATM/MODEL/ENDMDL/TER records.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written.

        file_handle : file handle
            The handle of the file where to write the structure.

        write_models_records : ``bool``
            Write the MODEL and ENDMDL records for each model found
            in the structure.

            You can set it to ``False`` only if you have one model
            in the structure. 
        """

        # Get the format strings for different types of records.
        fmt_atom = _defaults.PDB_SECTIONS["ATOM"]
        fmt_hetatm = _defaults.PDB_SECTIONS["HETATM"]
        fmt_model = _defaults.PDB_SECTIONS["MODEL"]
        fmt_endmdl = _defaults.PDB_SECTIONS["ENDMDL"]
        fmt_ter = _defaults.PDB_SECTIONS["TER"]

        #-------------------------------------------------------------#

        # For each model and associated data
        for i, (mod, mod_d) in enumerate(struct.atom_data.items()):

            # If the models' records need to be written
            if write_models_records:

                # Write out the corresponding MODEL record.
                file_handle.write(fmt_model.format(mod))

            # Otherwise
            else:

                # If there are multiple models in the structure
                if len(struct.atom_data) > 1:

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

                    # Write out the records anyway.
                    file_handle.write(fmt_model.format(mod))

            #---------------------------------------------------------#

            # For each chain and associated data
            for ch, ch_d in mod_d["_items"].items():

                # For each segment and associated data
                for seg, seg_d in ch_d["_items"].items():

                    # For each residue and associated data
                    for res, res_d in seg_d["_items"].items():

                        # If we need to write a HETATM record
                        if res_d["_attributes"]["is_het"]:
                            
                            # Select the corresponding format string.
                            fmt_str = fmt_hetatm

                        # If we need to write an ATOM record
                        else:

                            # Select the corresponding format string.
                            fmt_str = fmt_atom

                        # For each atom and associated data
                        for atom, atom_d in \
                            res_d["_items"].items():

                            # Get the atom's attributes.
                            atom_attrs = atom_d["_attributes"]

                            # Write out the ATOM/HETATM record.
                            file_handle.write(\
                              fmt_str.format(\
                                atom,
                                atom_attrs["label_atom_id"],
                                atom_attrs["label_alt_id"],
                                res[2],
                                ch,
                                res[0],
                                res[1],
                                atom_attrs["Cartn_x"],
                                atom_attrs["Cartn_y"],
                                atom_attrs["Cartn_z"],
                                atom_attrs["occupancy"],
                                atom_attrs["B_iso_or_equiv"],
                                atom_attrs["type_symbol"],
                                atom_attrs["pdbx_formal_charge"],
                                atom_attrs["label_entity_id"]))

                            # Save the current atom serial to write
                            # the TER record for the current chain
                            # (plus one since the TER record has a
                            # serial one unit higher than the last
                            # atom of the chain).
                            ter_atom_serial = atom + 1

                        # Save the current residue's sequence number to
                        # write the TER record for the current chain.
                        ter_res_seq = res[0]

                        # Save the current residue's name to write the
                        # TER record for the current chain.
                        ter_res_name = res[2]

                # Save the current chain identifier to write the TER
                # record for the current chain.
                ter_chain_id = ch

                # Write out the TER record for the current chain.
                file_handle.write(fmt_ter.format(ter_atom_serial,
                                                 ter_res_name,
                                                 ter_chain_id,
                                                 ter_res_seq))

            #---------------------------------------------------------#

            # If the models' records need to be written
            if write_models_records:
                
                # Write out the ENDMDL record for the current model.
                file_handle.write(fmt_endmdl)

            # Otherwise
            else:

                # If there are multiple models in the structure
                if len(struct.atom_data) > 1:

                    # Write them out anyway (do not warn the user here
                    # since they were already warned before).
                    file_handle.write(fmt_model.format(mod))


    def _write_conect_records(self,
                              struct,
                              file_handle):
        """Write CONECT records.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written.

        file_handle : file handle
            The handle of the file where to write the structure.
        """

        # If there are no connectivity data
        if not struct.conect_data:

            # Warn the user.
            warnstr = \
                "No connectivity data found in the structure."
            logger.warning(warnstr)

            # Return.
            return

        #-------------------------------------------------------------#

        # Define a helper function to sort the records before writing
        # them.
        def get_sorted_records(conect_data):

            # Initialize an empty dictionary to store the records,
            # sorted according to the serial number of the atoms
            # they refer to.
            sorted_records = defaultdict(lambda: defaultdict(dict))

            # For each model
            for mod in conect_data:

                # For each connectivity record in the model
                for atom1, rec in conect_data[mod].items():

                    # Get the serial number for the first atom.
                    atom1_serial = \
                        struct.atom_serial(path = (mod,) + atom1)

                    # Initialize an empty list to store the bonded
                    # atoms, sorted in ascending order of serial
                    # number.
                    sorted_bonded_atoms = []

                    # For each atom in the connectivity record
                    for atom2, data in rec.items():

                        # Get the serial number for the second atom.
                        atom2_serial = \
                            struct.atom_serial(path = (mod,) + atom2)

                        # If the bond is a hydrogen bond
                        if data["conn_type_id"] == "hydrog":

                            # Ignore it.
                            continue

                        # Get how many time we have to write the atom
                        # from the bond order.
                        times_to_write = \
                            _defaults.STRUCT_BOND_ORDERS[\
                                data["pdbx_value_order"]]

                        # For each time we have to write the atom
                        for time in range(times_to_write):

                            # Add the bonded atom to the list.
                            sorted_bonded_atoms.append(atom2_serial)

                    # Sort the list and save it as a record for the
                    # current atom.
                    sorted_records[mod][atom1_serial] = \
                        sorted(sorted_bonded_atoms)

            # Sort the atoms in the dictionary and return it.
            return {mod : \
                        {atom1 : data for atom1, data \
                         in sorted(sorted_records[mod].items())} \
                    for mod in sorted_records}

        #-------------------------------------------------------------#

        # Get the records, sorted.
        sorted_records = \
            get_sorted_records(conect_data = struct.conect_data)

        #-------------------------------------------------------------#

        # For each model in the sorted records
        for mod in sorted_records:

            # For each connectivity record in the model
            for atom1_serial, rec in sorted_records[mod].items():

                # Write out the header.
                file_handle.write("CONECT")

                # Write out the atom to which the record refers to.
                file_handle.write("{:5d}".format(atom1_serial))

                # For each atom in the connectivity record
                for atom2_serial in rec:

                    # Write out the atom.
                    file_handle.write(\
                        "{:5d}".format(atom2_serial))

                # Write out a newline character at the end of the
                # record.
                file_handle.write("\n")


    def _write_ssbond_records(self,
                              struct,
                              file_handle):
        """Write SSBOND records.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written.

        file_handle : file handle
            The handle of the file where to write the structure.
        """

        # If there are no connectivity data
        if not struct.conect_data:

            # Warn the user.
            warnstr = \
                "No connectivity data found in the structure."
            logger.warning(warnstr)

            # Return.
            return

        #-------------------------------------------------------------#

        # Get the format string for the SSBOND records.
        fmt_ssbond = _defaults.PDB_SECTIONS["SSBOND"]

        #-------------------------------------------------------------#

        # For each model
        for mod in struct.conect_data:

            # Set an index to number the disulfide bonds.
            bond_index = 1

            # Initialize an empty set to store the bonds already
            # written out.
            written = set()

            # For each connectivity record in the model
            for atom1, rec \
                in struct.conect_data[mod].items():

                # Get:
                # - The dentifier of the chain the first atom belongs
                #   to.
                # - The sequence number of the residue the first atom
                #   belongs to.
                # - The insertion code of the residue the first atom
                #   belongs to.
                # - The name of the residue the first atom belongs to.
                atom1_chain_id, atom1_res_seq, \
                    atom1_i_code, atom1_res_name = \
                        atom1[0], atom1[2][0], \
                        atom1[2][1], atom1[2][2]

                # For each bonded atom
                for atom2, data in rec.items():

                    # If it is not a disulfide bond or we have
                    # already written it out
                    if data["conn_type_id"] != "disulf" \
                    or (atom1, atom2) in written:

                        # Ignore it.
                        continue

                    # Get:
                    # - The dentifier of the chain the second atom
                    #   belongs to.
                    # - The sequence number of the residue the second
                    #   atom belongs to.
                    # - The insertion code of the residue the second
                    #   atom belongs to.
                    # - The name of the residue the second atom belongs
                    #   to.
                    atom2_chain_id, atom2_res_seq, \
                        atom2_i_code, atom2_res_name = \
                            atom2[0], atom2[2][0], \
                            atom2[2][1], atom2[2][2]  

                    # Get the symmetry operator for the first atom.
                    atom1_symop = data["ptnr1_symmetry"]

                    # Get the symmetry operator for the second atom.
                    atom2_symop = data["ptnr2_symmetry"]

                    # Get the bond's length.
                    bond_length = data["pdbx_dist_value"]

                    # Write out the record.
                    file_handle.write(\
                        fmt_ssbond.format(bond_index,
                                          atom1_res_name,
                                          atom1_chain_id,
                                          atom1_res_seq,
                                          atom1_i_code,
                                          atom2_res_name,
                                          atom2_chain_id,
                                          atom2_res_seq,
                                          atom2_i_code,
                                          atom1_symop,
                                          atom2_symop,
                                          bond_length))

                    # Write out a newline character at the end of the
                    # record.
                    file_handle.write("\n")

                    # Update the set of 'written out' bonds.
                    written.add((atom1, atom2))
                    written.add((atom2, atom1))

                    # Update the bond's index.
                    bond_index += 1


    ###################################################################


    def write(self,
              struct,
              file,
              write_conect_records = True,
              write_ssbond_records = True,
              write_models_records = True):
        """Write a structure to a PDB file.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure to be written.

        file : ``str``
            The file where to write the structure.

        write_conect_records : ``bool``, ``True``
            Write the CONECT data from the connectivity data
            associated with the structure, if any are present.

        write_ssbond_records : ``bool``, ``True``
            Write the SSBOND data from the connectivity data
            associated with the structure, if any are present.

        write_models_records : ``bool``, ``True``
            Write the MODEL and ENDMDL records for each model found
            in the structure.

            You can set it to ``False`` only if you have one model
            in the structure. 
        """

        # Open the output PDB file.
        with open(file, "w") as out:

            # If we need to write out SSBOND records
            if write_ssbond_records:

                # Write out the records.
                self._write_ssbond_records(\
                    struct = struct,
                    file_handle = out)

            #---------------------------------------------------------#

            # Write out the atomic coordinates.
            self._write_atom_records(\
                struct = struct,
                file_handle = out,
                write_models_records = write_models_records)

            #---------------------------------------------------------#

            # If we need to write out CONECT records
            if write_conect_records:

                # Write out the records.
                self._write_conect_records(\
                    struct = struct,
                    file_handle = out)

            #---------------------------------------------------------#

            # Write out the end of the PDB file.
            out.write("END   \n")
