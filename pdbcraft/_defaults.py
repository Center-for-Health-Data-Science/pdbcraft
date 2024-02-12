#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    _defaults.py
#
#    Private default/constant values.
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


#-------------------- Structure-related defaults ---------------------#


# Which residues are considered protein, DNA, or RNA residues
STRUCT_RESNAMES = \
  {# Protein residues
   "protein" : \
      {"ALA", "CYS", "ASP", "GLU",
       "PHE", "GLY", "HIS", "ILE",
       "LYS", "LEU", "MET", "ASN",
       "PRO", "GLN", "ARG", "SER",
       "THR", "VAL", "TRP", "TYR",
       "HIE", "HID", "HIP"},
   # DNA residues
   "dna" : \
      {"A", "C", "G", "U"},
   # RNA residues
   "rna" : \
      {"DA", "DC", "DG", "DT"}}

# Which keywords are supported to indicate a bond's order in
# the connectivity data. These keywords constitute the controlled
# vocabulary for the '_struct_conn.pdbx_value_order' item of
# mmCIF files
STRUCT_BOND_ORDERS = \
  {"sing" : 1, "doub" : 2, "trip" : 3, "quad" : 4}

# Which keywords are supported to indicate a bond's type in
# the connectivity data. Thess keywords constitute the controlled
# vocabulary for the '_struct_conn.conn_type_id' item of
# mmCIF files
STRUCT_BOND_TYPES = \
    {"covale", "disulf", "hydrog", "metalc"}


#---------------------- FASTA-related defaults -----------------------#


# A mapping between residues' three-letter names and residues'
# one-letter names (canonical residues + histidine
# protonation states)
FASTA_RESNAMES_3TO1 = \
    {"ALA" : "A", "CYS" : "C", "ASP" : "D", "GLU" : "E",
     "PHE" : "F", "GLY" : "G", "HIS" : "H", "ILE" : "I",
     "LYS" : "K", "LEU" : "L", "MET" : "M", "ASN" : "N",
     "PRO" : "P", "GLN" : "Q", "ARG" : "R", "SER" : "S",
     "THR" : "T", "VAL" : "V", "TRP" : "W", "TYR" : "Y",
     "HIE" : "H", "HID" : "H", "HIP" : "H"}


#----------------------- PDB-related defaults ------------------------#


# The format strings to use for each record type in PDB files
PDB_SECTIONS = \
    {"ATOM" : \
        "ATOM  {:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}   " \
        "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          " \
        "{:>2s}{:2s}\n",
     "HETATM" : \
        "HETATM{:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}   " \
        "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          " \
        "{:>2s}{:2s}\n",
     "MODEL" : \
        "MODEL     {:>4}\n",
     "ENDMDL" : \
        "ENDMDL\n",
     "TER" : \
        "TER   {:5d}      {:3s} {:1s}{:4d}\n",
     "SSBOND" : \
        "SSBOND {:>3d} {:<3s} {:>1s}{:>4d}{:>1s}   " \
        "{:<3s} {:>1s}{:>4d}{:>1s}                      " \
        "{:>5s} {:>5s} {:>5.2f}"}


#---------------------- mmCIF-related defaults -----------------------#

# The directive indicating the beginning of a data block
# in mmCIF files
MMCIF_DATA_DIRECTIVE = "data_"

# The directive indicating the beginning of a loop section
# in mmCIF files
MMCIF_LOOP_DIRECTIVE = "loop_"

# The characters used as quotes in mmCIF files
MMCIF_QUOTE_CHARACTERS = ["'", '"']

# The characters used to delimit white spaces in mmCIF files
MMCIF_WHITESPACE_CHARACTERS = [" ", "\t"]

# The characters used to specify unassigned values in mmCIF files
MMCIF_UNASSIGNED_CHARACTERS = [".", "?"]

# The character used at the beginning of lines to indicate
# comments in mmCIF files
MMCIF_COMMENT_CHARACTER = "#"

# The character used at the beginning of lines to indicate
# a multi-line quoted string
MMCIF_MULTILINE_CHARACTER = ";"

# The 'MMCIF_CATEGORIES' dictionary is a nested dictionary with
# two layers replicating the structure of mmCIF files' "data
# categories" and "data items". The keys of the outer dictionary
# correspond to mmCIF files' data categories, while the keys
# of each inner dictionary are names of data items of the
# corresponding data category.
#
# Each key in the inner dictionary is associated with a tuple
# with several values:
#
# - The first element defines the type of value stored by the data
#   item.
#
# - The second element defines, for data items whose value is a
#   floating point number, the number of decimal digits to use when
#   writing out the value. 'None' is used when the values of the
#   data item are strings or integers.
#
# - The third element defines either a controlled vocabulary to
#   be used to check the values of the data item, or the numerical
#   boundaries for the values, if they are numbers. This element can
#   be either a set of strings constituting the controlled vocabulary,
#   a tuple pointing to another data category and data item where
#   the controlled vocabulary is defined, or a list whose first
#   item is the lower boundary for the values, and the second item
#   is the higher boundary. 'None' is used when neither a controlled
#   vocabulary nor boundaries are defined.
#
# - The fourth element defines whether the item must be present
#   in the mmCIF file.
#
# - The fifth element defines whether there is a default value for
#   the item when its value is undefined.
MMCIF_CATEGORIES = \
    {\

     #----------------------------------------------------------------#

     # The fields reported are those flagged as "used in current
     # PDB entries" in at least 50% of entries in:
     #
     # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/
     # Categories/atom_type.html
     "_atom_type" : {\

         # The code used to identify the atom species (singular
         # or plural) representing this atom type. Normally this
         # code is the element symbol
         "symbol" : \
            (str, None, None, True, None),
      },

     #----------------------------------------------------------------#

     # The fields reported are those flagged as "used in current
     # PDB entries" in:
     #
     # https://mmcif.wwpdb.org/dictionaries/
     # mmcif_pdbx_v50.dic/Categories/atom_site.html
     #
     # We excluded the '_pdbx_refine_tls_group.id' data item
     "_atom_site" : {\
         
         # The group of atoms to which the atom site belongs.
         # This data item is provided for compatibility with the
         # original Protein Data Bank format, and only for that
         # purpose
         "group_PDB" : \
            (str, None, {"ATOM", "HETATM"}, True, None),
         
         # The value of _atom_site.id must uniquely identify a
         # record in the '_atom_site' list
         "id" : \
            (int, None, None, True, None),
         
         # A pointer to 'symbol' in the '_atom_type' category
         "type_symbol" : \
            (str, None, ("_atom_type", "symbol"), False, None),
         
         # A pointer to 'atom_id' in the '_chem_comp_atom' category
         "label_atom_id" : \
            (str, None, ("_chem_comp_atom", "atom_id"), True, None),
         
         # A pointer to 'id' in the '_atom_sites_alt' category
         "label_alt_id" : \
            (str, None, None, True, None),
         
         # A pointer to 'id' in the '_chem_comp' category
         "label_comp_id" : \
            (str, None, ("_chem_comp", "id"), True, None),
         
         # A pointer to 'id' in the '_struct_asym' category
         "label_asym_id" : \
            (str, None, None, True, None),
         
         # A pointer to 'id' in the '_entity' category
         "label_entity_id" : \
            (int, None, ("_entity", "id"), True, None),
         
         # A pointer to 'num' in the '_entity_poly_seq' category
         "label_seq_id" : \
            (int, None, ("_entity_poly_seq", "num"), True, None),
         
         # The PDB insertion code
         "pdbx_PDB_ins_code" : \
            (str, None, None, False, None),
         
         # The x atom-site coordinate in angstroms specified
         # according to a set of orthogonal Cartesian axes related
         # to the cell axes as specified by the description given
         # in '_atom_sites.Cartn_transform_axes'
         "Cartn_x" : \
            (float, 3, None, True, None),
         
         # The y atom-site coordinate in angstroms specified
         # according to a set of orthogonal Cartesian axes related
         # to the cell axes as specified by the description given
         # in '_atom_sites.Cartn_transform_axes'
         "Cartn_y" : \
            (float, 3, None, True, None),
         
         # The z atom-site coordinate in angstroms specified
         # according to a set of orthogonal Cartesian axes related
         # to the cell axes as specified by the description given
         # in '_atom_sites.Cartn_transform_axes'
         "Cartn_z": \
            (float, 3, None, True, None),
         
         # The fraction of the atom type present at this site.
         # The sum of the occupancies of all the atom types at
         # this site may not significantly exceed 1.0 unless it
         # is a dummy site
         "occupancy" : \
            (float, 2, None, False, None),
         
         # The isotropic atomic displacement parameter,
         # or equivalent isotropic atomic displacement parameter,
         # B_eq, calculated from the anisotropic displacement
         # parameters
         "B_iso_or_equiv" : \
            (float, 2, None, False, None),
         
         # The net integer charge assigned to this atom. This is
         # the formal charge assignment normally found in
         # chemical diagrams
         "pdbx_formal_charge" : \
            (float, 2, [-8, 8], False, None),
         
         # An alternative identifier for '_atom_site.label_seq_id'
         # that may be provided by an author in order to match the
         # identification used in the publication that describes
         # the structure
         "auth_seq_id" : \
            (int, None, None, False, None),
         
         # An alternative identifier for '_atom_site.label_comp_id'
         # that may be provided by an author in order to match the
         # identification used in the publication that describes
         # the structure.
         "auth_comp_id" : \
            (str, None, None, False, None),
         
         # An alternative identifier for '_atom_site.label_asym_id'
         # that may be provided by an author in order to match the
         # identification used in the publication that describes
         # the structure
         "auth_asym_id" : \
            (str, None, None, False, None),
         
         # An alternative identifier for '_atom_site.label_atom_id'
         # that may be provided by an author in order to match the
         # identification used in the publication that describes
         # the structure
         "auth_atom_id" : \
            (str, None, None, False, None),
         
         # The PDB model number
         "pdbx_PDB_model_num" : \
            (int, None, [0.0, float("inf")], False, None),
         
         # The standard uncertainty (estimated standard deviation)
         # of '_atom_site.Cartn_x'
         "Cartn_x_esd" : \
            (float, 2, None, False, None),
         
         # The standard uncertainty (estimated standard deviation)
         # of '_atom_site.Cartn_y'
         "Cartn_y_esd" : \
            (float, 2, None, False, None),
         
         # The standard uncertainty (estimated standard deviation)
         # of '_atom_site.Cartn_z'
         "Cartn_z_esd" : \
            (float, 2, None, False, None),
         
         # The standard uncertainty (estimated standard deviation)
         # of '_atom_site.occupancy'.
         "occupancy_esd" : \
            (float, 2, None, False, None),
         
         # The standard uncertainty (estimated standard deviation)
         # of '_atom_site.B_iso_or_equiv'
         "B_iso_or_equiv_esd" : \
            (float, 2, None, False, None),
         
         # A standard code to signal whether the site coordinates
         # have been determined from the intensities or calculated
         # from the geometry of surrounding sites, or have been
         # assigned dummy values.
         "calc_flag" : \
            (str, None, {"c", "calc", "d", "dum"}, False, None),
         },

     #----------------------------------------------------------------#

     # The fields reported are those flagged as "used in current
     # PDB entries" in:
     #
     # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/
     # Categories/atom_sites_alt.html
     "_atom_sites_alt" : {\

         # This value must uniquely identify a record in the
         # '_atom_sites_alt' list
         "id" : \
            (str, None, None, True, None),
         },

     #----------------------------------------------------------------#

     # The fields reported are those flagged as "used in current
     # PDB entries" in:
     #
     # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/
     # Categories/chem_comp.html
     "_chem_comp" : {\

         # This value must uniquely identify each item in the
         # '_chem_comp' list
         "id" : \
            (str, None, None, True, None),
         
         # The formula for the chemical component
         "formula" : \
            (str, None, None, False, None),
         
         # Formula mass in daltons of the chemical component
         "formula_weight" : \
            (float, None, None, [1.0, float("inf"), False, None]),
         
         # 'yes' indicates that this is a 'standard' monomer,
         # 'no' indicates that it is 'nonstandard'
         "mon_nstd_flag" : \
            (str, None, None, {"n", "no", "y", "yes"}, False, None),
         
         # The full name of the component
         "name" : \
            (str, None, None, False, None),
         
         # Synonym list for the component
         "pdbx_synonyms" : \
            (str, None, None, False, None),
         
         # For standard polymer components, the type of the monomer
         "type" : \
            (str, None, {"D-beta-peptide, C-gamma linking",
                         "D-gamma-peptide, C-delta linking",
                         "D-peptide COOH carboxy terminus",
                         "D-peptide NH3 amino terminus",
                         "D-peptide linking",
                         "D-saccharide",
                         "D-saccharide, alpha linking",
                         "D-saccharide, beta linking",
                         "DNA OH 3 prime terminus",
                         "DNA OH 5 prime terminus",
                         "DNA linking",
                         "L-DNA linking",
                         "L-beta-peptide, C-gamma linking",
                         "L-gamma-peptide, C-delta linking",
                         "L-peptide COOH carboxy terminus",
                         "L-peptide NH3 amino terminus",
                         "L-peptide linking",
                         "L-saccharide",
                         "L-saccharide, alpha linking",
                         "L-saccharide, beta linking",
                         "RNA OH 3 prime terminus",
                         "RNA OH 5 prime terminus",
                         "RNA linking",
                         "non-polymer",
                         "other",
                         "peptide linking",
                         "peptide-like",
                         "saccharide"},
               False, None),
         },

     #----------------------------------------------------------------#

     # The fields reported are those flagged as "used in current
     # PDB entries" in:
     #
     # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/
     # Categories/chem_comp_atom.html
     "_chem_comp_atom" : {\

         # This value must uniquely identify each atom in
         # each monomer in the '_chem_comp_atom' list
         "atom_id" : \
            (str, None, None, True, None),
         
         # A pointer to 'id' in the '_chem_comp' category
         "comp_id" : \
            (str, None, ("_chem_comp", "id"), True, None),
         
         # The code used to identify the atom species representing
         # this atom type. Normally, this code is the element
         # symbol
         "type_symbol" : \
            (str, None, None, False, None),
         },

     #----------------------------------------------------------------#

     # The fields reported are those flagged as "used in current
     # PDB entries" in:
     #
     # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/
     # Categories/entity.html
     "_entity" : {\

         # This value must uniquely identify a record in the
         # '_entity' category
         "id" : \
            (str, None, None, True, None),
         
         # A description of special aspects of the entity
         "details" : \
            (str, None, None, False, None),
         
         # The formula mass in daltons of the entity
         "formula_weight" : \
            (float, None, [1.0, float("inf")], False, None),
         
         # The enzyme Commission (EC) number(s)
         "pdbx_ec" : \
            (str, None, None, False, None),
         
         # The entity fragment description(s)
         "pdbx_fragment" : \
            (str, None, None, False, None),
         
         # Details about any entity mutation(s)
         "pdbx_mutation" : \
            (str, None, None, False, None),
         
         # A placeholder for the number of molecules of the
         # entity in the entry
         "pdbx_number_of_molecules" : \
            (int, None, None, False, None),
         
         # The method by which the sample for the entity was
         # produced
         "src_method" : \
            (str, None, {"man", "nat", "syn"}, False, None),
         
         # This value defines the type of the entity
         "type" : \
            (str, None, {"branched", "macrolide", "non-polymer",
                         "polymer", "water"},
             False, None),
         },

     #----------------------------------------------------------------#

     # The fields reported are all those in "used in current
     # PDB entries" in:
     #
     # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/
     # Categories/entity_poly_seq.html
     "_entity_poly_seq" : {\

         # This value must uniquely and sequentially identify
         # a record in the '_entity_poly_seq' list
         "num" : \
            (int, None, None, True, None),

         # A pointer to 'id' in the '_entity' category
         "entity_id" : \
            (str, None, ("_entity", "id"), True, None),
         
         # A pointer to 'id' in the '_chem_comp' category
         "mon_id" : \
            (str, None, ("_chem_comp", "id"), True, None),
         
         # A flag to indicate whether this monomer in the
         # polymer is heterogeneous in sequence.
         "hetero" : \
            (str, None, {"n", "no", "y", "yes"}, False, None),
         },

     #----------------------------------------------------------------#
     
     # The fields reported are those flagged as "used in current
     # PDB entries" in:
     #
     # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/
     # Categories/struct_conn.html
     "_struct_conn" : {\

         # The value of '_struct_conn.id' must uniquely identify
         # a record in the '_struct_conn' category
         "id" : \
            (str, None, None, True, None),
         
         # A pointer to 'id' in the '_struct_conn_type' category
         "conn_type_id" : \
            (str, None, ("_struct_conn_type", "id"), True, None),
         
         # This data item identifies if the linkage has displaced
         # leaving atoms on both, one or none of the connected
         # atoms forming the linkage. Leaving atoms are defined
         # within their chemical defintions of each
         # connected component
         "pdbx_leaving_atom_flag" : \
            (str, None, {"both", "none", "one"}, False, None),
         
         # A pointer to 'auth_asym_id' in the '_atom_site' category
         "ptnr1_label_asym_id" : \
            (str, None, ("_atom_site", "label_asym_id"), True, None),
         
         # A pointer to 'label_comp_id' in the '_atom_site' category
         "ptnr1_label_comp_id" : \
            (str, None, ("_atom_site", "label_comp_id"), True, None),
         
         # A pointer to 'label_seq_id' in the '_atom_site' category
         "ptnr1_label_seq_id" : \
            (int, None, ("_atom_site", "label_seq_id"), True, None),
         
         # A pointer to 'atom_id' in the '_chem_comp_atom' category
         "ptnr1_label_atom_id" : \
            (str, None, ("_chem_comp_atom", "atom_id"), True, None),
         
         # A pointer to 'label_alt_id' in the '_atom_site' category
         "pdbx_ptnr1_label_alt_id" : \
            (str, None, ("_atom_site", "label_alt_id"), True, None),
         
         # A pointer to 'pdbx_PDB_ins_code' in the '_atom_site'
         # category
         "pdbx_ptnr1_PDB_ins_code" : \
            (str, None, ("_atom_site", "pdbx_PDB_ins_code"), 
             True, None),
         
         # This value describes the symmetry operation that
         # should be applied to the atom set specified by
         # '_struct_conn.ptnr1_label*' to generate the
         # first partner in the structure connection
         "ptnr1_symmetry" : \
            (str, None, None, False, None),
         
         # A pointer to 'auth_asym_id' in the '_atom_site' category
         "ptnr2_label_asym_id" : \
            (str, None, ("_atom_site", "label_asym_id"), True, None),
         
         # A pointer to 'label_comp_id' in the '_atom_site' category
         "ptnr2_label_comp_id" : \
            (str, None, ("_atom_site", "label_comp_id"), True, None),
         
         # A pointer to 'label_seq_id' in the '_atom_site' category
         "ptnr2_label_seq_id" : \
            (int, None, ("_atom_site", "label_seq_id"), True, None),
         
         # A pointer to 'atom_id' in the '_chem_comp_atom' category
         "ptnr2_label_atom_id" : \
            (str, None, ("_chem_comp_atom", "atom_id"), True, None),
         
         # A pointer to 'label_alt_id' in the '_atom_site' category
         "pdbx_ptnr2_label_alt_id" : \
            (str, None, ("_atom_site", "label_alt_id"), True, None),
         
         # A pointer to 'pdbx_PDB_ins_code' in the '_atom_site'
         # category
         "pdbx_ptnr2_PDB_ins_code" : \
            (str, None, ("_atom_site", "pdbx_PDB_ins_code"), 
             True, None),
         
         # A pointer to 'auth_asym_id' in the '_atom_site' category
         # for partner 1
         "ptnr1_auth_asym_id" : \
            (str, None, None, False, None),
         
         # A pointer to 'auth_comp_id' in the '_atom_site' category
         # for partner 1
         "ptnr1_auth_comp_id" : \
            (str, None, None, False, None),
         
         # A pointer to 'auth_seq_id' in the '_atom_site' category
         # for partner 1
         "ptnr1_auth_seq_id" : \
            (int, None, None, False, None),
         
         # A pointer to 'auth_asym_id' in the '_atom_site' category
         # for partner 2
         "ptnr2_auth_asym_id" : \
            (str, None, None, False, None),
         
         # A pointer to 'auth_comp_id' in the '_atom_site' category
         # for partner 2
         "ptnr2_auth_comp_id" : \
            (str, None, None, False, None),
         
         # A pointer to 'auth_seq_id' in the '_atom_site' category
         # for partner 2
         "ptnr2_auth_seq_id" : \
            (int, None, None, False, None),
         
         # This value describes the symmetry operation that
         # should be applied to the atom set specified by
         # '_struct_conn.ptnr2_label*' to generate the
         # second partner in the structure connection
         "ptnr2_symmetry" : \
            (str, None, None, False, None),
         
         # A description of special aspects of the connection 
         "details" : \
            (str, None, None, False, None),
         
         # The distance value for this contact
         "pdbx_dist_value" : \
            (float, 3, None, False, None),
         
         # The chemical or structural role of the interaction
         "pdbx_role" : \
            (str, None, {"C-Mannosylation", "N-Glycosylation",
                         "O-Glycosylation", "S-Glycosylation"},
             False, None),
         
         # The chemical bond order associated with the specified
         # atoms in this contact
         "pdbx_value_order" : \
            (str, None, {"doub", "quad", "sing", "trip"}, 
             True, "sing"),
         },

     #----------------------------------------------------------------#

     # The fields reported are those in:
     #
     # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/
     # Categories/struct_conn_type.html
     "_struct_conn_type" : {\
         
         # The chemical or structural type of the interaction
         "id" : \
            (str, None, {"covale", "covale_base", 
                         "covale_phosphate", "covate_sugar",
                         "disulf", "hydrog", "metalc",
                         "mismat", "modres", "saltbr"},
             True, None),
         
         # The criteria used to define the interaction
         "criteria" : \
            (str, None, None, False, None),
         
         # A reference that specifies the criteria used to define
         # the interaction
         "reference" : \
            (str, None, None, False, None),
         },
      }