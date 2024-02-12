#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    fastawriters.py
#
#    Writers for FASTA files.
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
import logging as log
import os
# pdbcraft
from . import _defaults


# Get the module's logger
logger = log.getLogger(__name__)


class FASTAWriter:

    """
    A writer for FASTA files.
    """


    #------------------------ Public methods -------------------------#


    def write(self,
              struct,
              file,
              split_models = False,
              split_chains = False,
              disc_chains_mode = "join_with_gaps",
              res_i_code = " ",
              gap_char = "-",
              wrap_at = None,
              resnames_3to1 = None):
        """Write a FASTA file containing the sequence corresponding
        to a ``Structure``.

        Parameters
        ----------
        struct : ``pdbcraft.structure.Structure``
            The structure whose sequence needs to be written.

        file : ``str``
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
        The file will have the name provided in the ``file``
        option.

        If ``split_models = False`` and ``split_chains = True``,
        one FASTA file per chain will be written, containing
        the FASTA sequence of that chain in all models.
        Assuming the ``file`` option has the form 
        ``{file_name}.{extension}``, each file will be named
        ``{file_name}_{chain}.{extension}``.

        If ``split_models = True`` and ``split_chains = False``,
        one FASTA file per model will be written, containing
        the FASTA sequences of the chains in that model.
        Assuming the ``file`` option has the form 
        ``{file_name}.{extension}``, each file will be named
        ``{file_name}_{model}.{extension}``.

        If ``split_models = True`` and ``split_chains = True``,
        one FASTA file per chain in each model will be written,
        containing the FASTA sequence of that chain in that
        model.
        Assuming the ``file`` option has the form 
        ``{file_name}.{extension}``, each file will be named
        ``{file_name}_{model}_{chain}.{extension}``.
        """

        # Create an empty set to store the files that have
        # been created
        files_created = set()

        #-------------------------------------------------------------#

        # Get the file's name (including the path leading to
        # it, if any) and the file's extension, separately
        fasta_name, fasta_ext = os.path.splitext(file)
        
        #-------------------------------------------------------------#

        # If the structure has no name
        if struct.name is None:

            # Set the prefix for the FASTA headers to an empty string
            prefix = ""

        # Otherwise
        else:

            # The prefix will contain the name of the structure
            prefix = f"{struct.name}_"

        #-------------------------------------------------------------#

        # If no mapping between the residues' three-letter
        # names and their one-letter names was passed
        if resnames_3to1 is None:

            # Set it to the default one
            resnames_3to1 = _defaults.FASTA_RESNAMES_3TO1

        #-------------------------------------------------------------#

        # If we do not need to split neither models nor chains
        if not split_models and not split_chains:

            # Set the name of the only FASTA file that will be
            # written
            file_name = file

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

        #-------------------------------------------------------------#

        # For each model and associated data in the structure
        for mod, mod_data in struct.atom_data.items():

            #---------------------------------------------------------#

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

                    # Add the name of the file to the set of files
                    # created
                    files_created.add(file_name)

                # Open the file in the selected mode
                out = open(file_name, file_mode)

            #---------------------------------------------------------#
            
            # For each chain and associated data in the current model
            for ch, ch_data in mod_data["_items"].items():

                # Inizialite a list to store the FASTA sequences
                # (and corresponding headers) for the current chain
                # in the current model
                fasta_seqs = [[f">{prefix}{mod}_{ch}_", ""]]

                #-----------------------------------------------------#

                # For each chain segment and associated data
                for seg, seg_data in ch_data["_items"].items():

                    # Initialize the variable storing the
                    # previously visited residue's sequence
                    # number to None
                    prev_res_num = None

                    #-------------------------------------------------#

                    # For each residue in the segment
                    for i, res in \
                        enumerate(seg_data["_items"].items()):

                        # Get the residue's sequence number,
                        # insertion code, and name
                        res_num, i_code, res_name_3 = res

                        # If the insertion code is not the one
                        # of the residues to be considered
                        if i_code != res_i_code:

                            # Skip the current residue
                            continue

                        #---------------------------------------------#

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

                        #---------------------------------------------#

                        # If we are at the first residue of the
                        # segment
                        if i == 0:

                            # Update the current chain's sequence
                            fasta_seqs[-1][1] += res_name_1

                            # Update the current chain's header
                            fasta_seqs[-1][0] += f"{res_num}"
                            
                            # Set the previous residue to the
                            # current residue
                            prev_res_num = res_num
                            
                            # Go to the next residue
                            continue

                        #---------------------------------------------#

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

                            #-----------------------------------------#

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

                            #-----------------------------------------#

                            # If we should split discontinuous portions
                            # of the chain
                            elif disc_chains_mode == "split":

                                # Update the current chain's sequence
                                fasta_seqs.append(\
                                    [f">{struct.name}_{res_num}",
                                     res_name_1])

                                # Update the previous chain's header
                                # by adding the end of the previous
                                # chain
                                fasta_seqs[-2][0] += \
                                    f"-{prev_res_num}"

                        #---------------------------------------------#

                        # If we are at a continuous point in the chain
                        else:

                            # Update the current chain's sequence
                            fasta_seqs[-1][1] += res_name_1

                        # If we are at the last residue of the
                        # segment
                        if i == len(seg_data) - 1:

                            # Update the current chain's header
                            # by adding the end of the current
                            # chain
                            fasta_seqs[-1][0] += f"_{res_num}"

                        # Update the 'previous residue number'
                        # variable
                        prev_res_num = res_num

                #-----------------------------------------------------#

                # If we need to split chains and not models
                if split_models and split_chains:

                    # Set the file's name
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

                #-----------------------------------------------------#

                # If we need to split chains but not models
                elif not split_models and split_chains:

                    # Set the file's name
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

                #-----------------------------------------------------#

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
