#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    util.py
#
#    Miscellanea utilities.
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
__doc__ = "Miscellanea utilities."


#######################################################################


# Import from the standard library.
import logging as log
import os
import tempfile
import urllib.request
# Import from 'pdbcraft'.
from . import (
    _defaults,
    mmcifparsers,
    mmcifwriters,
    pdbparsers,
    pdbwriters,
    structure,
    )


#######################################################################


# Get the module's logger.
logger = log.getLogger(__name__)


#######################################################################


# Set a dictionary mapping the parsers' names to the classes they
# represent and the options to be used with their 'parse()' method.
NAME2PARSER = \
    {"PDBParser" : \
        {"class" : pdbparsers.PDBParser,
         "options" : \
            {"only_alt_loc" : "A",
             "parse_conect_records" : True,
             "parse_ssbond_records" : True,
             "parse_link_records" : True,
             "is_extended_pdb" : False,
             "fix_residues_overflow" : True}},
     "MMCIFParser" : \
        {"class" : mmcifparsers.MMCIFParser,
         "options" : \
            {"only_alt_loc" : "A",
             "use_auth_chain_ids" : True,
             "use_auth_res_seqs" : True,
             "use_auth_res_names" : True,
             "use_auth_atom_names" : True,
             "strict" : False}}}

# Set a dictionary mapping the file types to the default parser that
# will be used for them and the parser's options.
FILETYPE2PARSER = \
    {"pdb" : \
        {"name" : "PDBParser",
         "options" : NAME2PARSER["PDBParser"]["options"]},
     "mmcif" : \
        {"name" : "MMCIFParser",
         "options" : NAME2PARSER["MMCIFParser"]["options"]}}

#---------------------------------------------------------------------#

# Set a dictionary mapping the writers' names to the classes they
# represent and the options to be used with their 'write()' method.
NAME2WRITER = \
    {"PDBWriter" : \
        {"class" : pdbwriters.PDBWriter,
         "options" : \
            {"write_conect_records" : True,
             "write_ssbond_records" : True,
             "write_models_records" : True}},
     "MMCIFWriter" : \
        {"class" : mmcifwriters.MMCIFWriter,
         "options" : \
            {"write_conect_data" : True}}}

# Set a dictionary mapping the file types to the default writer that
# will be used for them and the writer's options.
FILETYPE2WRITER = \
    {"pdb" : \
        {"name" : "PDBWriter",
         "options" : NAME2WRITER["PDBWriter"]["options"]},
     "mmcif" : \
        {"name" : "MMCIFWriter",
         "options" : NAME2WRITER["MMCIFWriter"]["options"]}}

#---------------------------------------------------------------------#

# Set the mapping between the source and where to retrieve the
# the structures from.
SOURCE2DB = \
    {# The Protein Data Bank
     "pdb" : \
        {# The generic URL to retrieve PDB or PDBx/MMCIF files
         # from the Protein Data Bank
         "url": "https://files.rcsb.org/view/{:s}.{:s}",
         # The name of the database 
         "db_name" : "Protein Data Bank"},
     # The AlphaFold Protein Structure Database
     "alphafold" : \
        {# The generic URL to retrieve PDB or PDBx/MMCIF files from
         # the AlphaFold Protein Structure Database
         "url" : "https://www.alphafold.ebi.ac.uk/files/{:s}.{:s}",
         # The name of the database
         "db_name" : "AlphaFold Protein Structure Database"},
    }

# Set the mapping between the file type and the corresponding file
# extension.
FILETYPE2FILEEXT = \
    {# PDB files
     "pdb" : ".pdb",
     # PDBx/MMCIf files
     "mmcif" : ".cif"}

# Set the mapping between a file's extension and the file type.
FILEEXT2FILETYPE = {v : k for k, v in FILETYPE2FILEEXT.items()}


#######################################################################


def load_structure(struct_file,
                   parser_name = None,
                   parser_options = None):
    """Load a structure from a PDB or PDBx/mmCIF file.

    Parameters
    ----------
    struct_file : ``str``
        The PDB or PDBx/mmCIF file.

    parser_name : ``str``, optional, \
        {``"PDBParser"``, ``"MMCIFParser"``}
        The name of the parser to be used to read in the structure.

        By default:

        - If the file is a PDB file, the
          :class:`pdbcraft.pdbparsers.PDBParser` parser will be used.
        - If the file is a mmCIF file, the
          :class:`pdbcraft.mmcifparsers.MMCIFParser` parser will be
          used.

    parser_options : ``dict``, optional
        The options to be passed to the ``parse()`` method
        of the PDB or PDBx/MMCIF parser defined by ``parser_name``.

        This option is ignored if no ``parser_name`` is passed.

    Returns
    -------
    struct : :class:`pdbcraft.structure.Structure`
        The structure.
    """

    # Get the file's extension.
    file_name, file_ext = os.path.splitext(struct_file)

    #-----------------------------------------------------------------#

    # If the user passed the name of a parser
    if parser_name is not None:

        # If the name of the parser is invalid
        if parser_name not in NAME2PARSER:

            # Set a string representing the valid names.
            valid_names_str = \
                ", ".join((f"'{n}'" for n in NAME2PARSER))

            # Raise an error
            errstr = \
                f"Invalid 'parser_name': '{parser_name}'. Valid " \
                f"parser names are: {valid_names_str}."
            raise ValueError(errstr)

        # If no options for the parser were provided
        if parser_options is None:

            # Set it to the default options.
            parser_options = NAME2PARSER[parser_name]["options"]

    # Otherwise
    else:

        # Fall back to the default parser for this specific type of
        # file.
        parser_name = \
            FILETYPE2PARSER[FILEEXT2FILETYPE[file_ext]]["name"]

        # Fall back to the default parser's options for this specific
        # type of file.
        parser_options = \
            FILETYPE2PARSER[FILEEXT2FILETYPE[file_ext]]["options"]

    #-----------------------------------------------------------------#

    # Create the parser.
    parser = NAME2PARSER[parser_name]["class"]()

    # Create a string contaning the parser's options.
    parser_options_str = \
        ", ".join(\
            [f"{opt} = '{val}'" if isinstance(val, str) \
             else f"{opt} = {val}" \
             for opt, val in parser_options.items()])

    # Inform the user about the parser that will be used and the
    # options for the 'parse()' method.
    infostr = \
        f"Reading the structure using the '{parser_name}' parser. " \
        f"Options used for the 'parse()' method: {parser_options_str}."
    logger.info(infostr)

    # Read in the structure.
    struct = parser.parse(file = struct_file,
                          **parser_options)

    #-----------------------------------------------------------------#

    # Return the structure.
    return struct


def fetch_structure(struct_id,
                    source = "pdb",
                    file_type = "mmcif",
                    save_file = False,
                    file = None,
                    parser_name = None,
                    parser_options = None):
    """Given the ID of a structure in either the Protein Data Bank
    or the AlphaFold Protein Structure Database, return a
    :class:`pdbcraft.structure.Structure` instance containing the
    structure.

    Parameters
    ----------
    struct_id : ``str``
        The structure's ID.

    source : ``str``, {``"pdb"``, ``"alphafold"``}, ``"pdb"``
        Where to fetch the structure from.

        If ``"pdb"``, the structure will be retrieved from the Protein
        Data Bank.

        If ``"alphafold"``, the structure will be retrieved from the
        AlphaFold Protein Structure Database.

    file_type : ``str``, {``"pdb"``, ``"mmcif"``}, ``"mmcif"``
        Whether to generate the :class:`pdbcraft.structure.Structure`
        instance from the structure's PDB (``"pdb"``) or PCBx/MMCIF
        (``"mmcif"``) file.

    save_file : ``bool``, ``"False"``
        Whether to save the file the structure was created from.

    file : ``str``, optional
        The file where to save the structure, if ``save_file`` is
        ``True``.

        If not provided, the file will be saved in the current
        working directory and named ``"{struct_id}.pdb"`` or
        ``"{struct_id}.cif"`` depending on the ``file_type``
        provided.

    parser_name : ``str``, optional, \
        {``"PDBParser"``, ``"MMCIFParser"``}
        The name of the parser to be used to read in the structure.

        By default:
        
        - If ``file_type`` is ``"pdb"``, the
          :class:`pdbcraft.pdbparsers.PDBParser` parser will be used.
        - If ``file_type`` is ``"mmcif"``, the
          :class:`pdbcraft.mmcifparsers.MMCIFParser` parser will be
          used.

    parser_options : ``dict``, optional
        The options to be passed to the ``parse()`` method
        of the PDB or PDBx/MMCIF parser defined by ``parser_name``.

        This option is ignored if no ``parser_name`` is passed.

    Returns
    -------
    struct : :class:`pdbcraft.structure.Structure`
        The structure.
    """

    # If the source is invalid
    if source not in SOURCE2DB:

        # Set a string representing the valid sources.
        valid_sources_str = \
            ", ".join((f"'{s}'" for s in SOURCE2DB))

        # Raise an error.
        errstr = \
            f"Invalid 'source': '{source}'. Valid sources are: " \
            f"{valid_sources_str}."
        raise ValueError(errstr)

    #-----------------------------------------------------------------#

    # If the file type is invalid
    if file_type not in FILETYPE2FILEEXT:

        # Set a string representing the valid file types.
        valid_types_str = \
            ", ".join((f"'{t}'" for t in FILETYPE2FILEEXT))

        # Raise an error
        errstr = \
            f"Invalid 'file_type': '{file_type}'. Valid file " \
            f"types are: {valid_types_str}."
        raise ValueError(errstr)

    #-----------------------------------------------------------------#

    # If the user passed the name of a parser
    if parser_name is not None:

        # If the name of the parser is invalid
        if parser_name not in NAME2PARSER:

            # Set a string representing the valid names.
            valid_names_str = \
                ", ".join((f"'{n}'" for n in NAME2PARSER))

            # Raise an error
            errstr = \
                f"Invalid 'parser_name': '{parser_name}'. Valid " \
                f"parser names are: {valid_names_str}."
            raise ValueError(errstr)

        # If no options for the parser were provided
        if parser_options is None:

            # Set it to the default options.
            parser_options = NAME2PARSER[parser_name]["options"]

    # Otherwise
    else:

        # Fall back to the default parser for this specific type of
        # file.
        parser_name = FILETYPE2PARSER[file_type]["name"]

        # Fall back to the default parser's options for this specific
        # type of file.
        parser_options = FILETYPE2PARSER[file_type]["options"]

    #-----------------------------------------------------------------#

    # If the source is 'pdb'
    if source == "pdb":

        # Capitalize the structure's ID, in case it isn't already.
        struct_id = struct_id.upper()

    #-----------------------------------------------------------------#

    # Open the URL.
    response = \
        urllib.request.urlopen(\
            source2url[source]["url"].format(\
                struct_id,
                FILETYPE2FILEEXT[file_type]))

    # If something went wrong
    if response.status != 200:

        # Raise an error.
        errstr = \
            "It was not possible to retrieve the structure with " \
            f"ID '{struct_id}' from the " \
            f"{SOURCE2DB['source']['db_name']}. Error code when " \
            f"retrieving the structure: {response.status}."
        raise Exception(errstr)

    # Read the data from the URL.
    data = response.read()

    #-----------------------------------------------------------------#

    # If we need to save the file
    if save_file:

        # If the user did not provide a custom file path
        if file is None:

            # Set it.
            file = \
                os.path.join(\
                    os.getcwd(),
                    f"{struct_id}{FILETYPE2FILEEXT[file_type]}")

        # Open the file.
        f = open(file, "r")

        # Write out the data.
        f.write(data)

        # Inform the user that the file was saved.
        infostr = \
            f"The structure with ID '{struct_id}' was successfully " \
            f"saved in '{file}'."
        logger.info(infostr)

    # Otherwise
    else:

        # Open a named temporary file.
        f = tempfile.NamedTemporaryFile(delete = False)

        # Write out the data.
        f.write(data)

        # Set the file's name to the name of the temporary file.
        file = f.name

    #-----------------------------------------------------------------#

    # Create the parser.
    parser = NAME2PARSER[parser_name]["class"]()

    # Create a string contaning the parser's options.
    parser_options_str = \
        ", ".join(\
            [f"{opt} = '{val}'" if isinstance(val, str) \
             else f"{opt} = {val}" \
             for opt, val in parser_options.items()])

    # Inform the user about the parser that will be used and the
    # options for the 'parse()' method.
    infostr = \
        f"Reading the structure using the '{parser_name}' parser. " \
        f"Options used for the 'parse()' method: {parser_options_str}."
    logger.info(infostr)

    # Read in the structure.
    struct = parser.parse(file = file,
                          **parser_options)

    #-----------------------------------------------------------------#

    # Close the file.
    f.close()

    #-----------------------------------------------------------------#

    # Return the structure.
    return struct


def write_structure(struct,
                    file,
                    writer_name = None,
                    writer_options = None):
    """Write a structure as a PDB or PDBx/mmCIF file.

    Parameters
    ----------
    struct : :class:`pdbcraft.structure.Structure`
        The structure.

    file : ``str``
        The file where to write the structure. The file's type will
        be guessed from the file's extension.

    writer_name : ``str``, optional, \
        {``"PDBWriter"``, ``"MMCIFWriter"``}
        The name of the writer to be used to write the structure.

        By default:

        - If the file is a PDB file, the
          :class:`pdbcraft.pdbwriters.PDBWriter` writer will be used.
        - If the file is a mmCIF file, the
          :class:`pdbcraft.mmcifwriters.MMCIFWriter` writer will be
          used.

    writer_options : ``dict``, optional
        The options to be passed to the ``write()`` method
        of the PDB or PDBx/MMCIF writer defined by ``writer_name``.

        This option is ignored if no ``writer_name`` is passed.
    """

    # Get the file's extension.
    file_name, file_ext = os.path.splitext(file)

    #-----------------------------------------------------------------#

    # If the user passed the name of a writer
    if writer_name is not None:

        # If the name of the writer is invalid
        if writer_name not in NAME2WRITER:

            # Set a string representing the valid names.
            valid_names_str = \
                ", ".join((f"'{n}'" for n in NAME2WRITER))

            # Raise an error
            errstr = \
                f"Invalid 'writer_name': '{writer_name}'. Valid " \
                f"writer names are: {valid_names_str}."
            raise ValueError(errstr)

        # If no options for the writer were provided
        if writer_options is None:

            # Set it to the default options.
            writer_options = NAME2WRITER[writer_name]["options"]

    # Otherwise
    else:

        # Fall back to the default writer for this specific type of
        # file.
        writer_name = \
            FILETYPE2WRITER[FILEEXT2FILETYPE[file_ext]]["name"]

        # Fall back to the default writer's options for this specific
        # type of file.
        writer_options = \
            FILETYPE2WRITER[FILEEXT2FILETYPE[file_ext]]["options"]

    #-----------------------------------------------------------------#

    # Create the writer.
    writer = NAME2WRITER[writer_name]["class"]()

    # Create a string contaning the writer's options.
    writer_options_str = \
        ", ".join(\
            [f"{opt} = '{val}'" if isinstance(val, str) \
             else f"{opt} = {val}" \
             for opt, val in writer_options.items()])

    # Inform the user about the writer that will be used and the
    # options for the 'write()' method.
    infostr = \
        f"Writing the structure using the '{writer_name}' writer. " \
        f"Options used for the 'write()' method: {writer_options_str}."
    logger.info(infostr)

    # Write the structure.
    writer.write(struct = struct,
                 file = file,
                 **writer_options)


def list_examples():
    """List all structures available as examples.

    Returns
    -------
    examples : ``list``
        A list containing the names of the files of the structures
        available as examples.
    """
    
    # Return the names of the files in the directory containing the
    # structures available as examples.
    return [f for f in os.listdir(_defaults.EXAMPLES_DIR)]


def load_example(example_id):
    """Load one of the structures available as examples.

    Parameters
    ----------
    example_id : ``str``
        The name of the file of a structure available as an example.

        You can see all available structures by calling the
        :func:`pdbcraft.util.list_examples` function.

    Returns
    -------
    struct : :class:`pdbcraft.structure.Structure`
        The structure.
    """

    # Get the list of available structures.
    examples = list_examples()

    #-----------------------------------------------------------------#

    # If the provided ID is not valid
    if example_id not in examples:

        # Set a string representing the available examples.
        examples_str = ", ".join([f"'{ex}'" for ex in examples])

        # Raise an error.
        errstr = \
            f"Example '{example_id}' is not available. " \
            f"Available examples are: {examples_str}."
        raise ValueError(errstr)

    #-----------------------------------------------------------------#

    # If the structure is a PDB structure
    if os.path.splitext(example_id)[1] == ".pdb":

        # We are going to load it using the standard PDB parser.
        parser = pdbparsers.PDBParser()

    #-----------------------------------------------------------------#

    # If the structure is a mmCIF parser
    elif os.path.splitext(example_id)[1] == ".cif":

        # We are going to load it using the standard mmCIF parser.
        parser = mmcifparsers.MMCIFParser()

    #-----------------------------------------------------------------#

    # Load the structure.
    struct = \
        parser.parse(\
            file = os.path.join(_defaults.EXAMPLES_DIR, example_id))

    #-----------------------------------------------------------------#

    # Return the structure.
    return struct
