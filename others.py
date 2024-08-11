#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
import shlex

MMCIF_SECTIONS_EXTENDED = \
    {\

     #----------------------------------------------------------------#

     # The fields reported are those flagged as "used in current
     # PDB entries" in at least 50% of entries in:
     # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/
     # Categories/atom_type.html
     "_atom_type" : {\

         # The Cromer-Mann scattering-factor coefficient a1 used to
         # calculate the scattering factors for this atom type
         "scat_Cromer_Mann_a1" : \
            (float, 6, None),
         # The Cromer-Mann scattering-factor coefficient a2 used to
         # calculate the scattering factors for this atom type
         "scat_Cromer_Mann_a2" : \
            (float, 6, None),
         # The Cromer-Mann scattering-factor coefficient a3 used to
         # calculate the scattering factors for this atom type
         "scat_Cromer_Mann_a3" : \
            (float, 6, None),
         # The Cromer-Mann scattering-factor coefficient a4 used to
         # calculate the scattering factors for this atom type
         "scat_Cromer_Mann_a4" : \
            (float, 6, None),
         # The Cromer-Mann scattering-factor coefficient b1 used to
         # calculate the scattering factors for this atom type
         "scat_Cromer_Mann_b1" : \
            (float, 6, None),
         # The Cromer-Mann scattering-factor coefficient b2 used to
         # calculate the scattering factors for this atom type
         "scat_Cromer_Mann_b2" : \
            (float, 6, None),
         # The Cromer-Mann scattering-factor coefficient b3 used to
         # calculate the scattering factors for this atom type
         "scat_Cromer_Mann_b3" : \
            (float, 6, None),
         # The Cromer-Mann scattering-factor coefficient b4 used to
         # calculate the scattering factors for this atom type
         "scat_Cromer_Mann_b4" : \
            (float, 6, None),
         # The Cromer-Mann scattering-factor coefficient c used to
         # calculate the scattering factors for this atom type
         "scat_Cromer_Mann_c" : \
            (float, 6, None),
         # Ahe reference to the source of the scattering factors
         # or scattering lengths used for this atom type
         "scat_source" : \
            (str, None, None),
      },
    }


def load_mmcif_dict(file):

    data = {}

    with open(file, "r") as f:

        parse = False
        is_multiline = False
        is_loop = False

        curr_multiline = ""

        curr_category = None
        curr_item = None

        for line_ix, line in enumerate(f, 1):

            line = line.lstrip()
            print(line)

            

            if line.startswith("save__"):
                item = line[5:].rstrip()
                curr_item = item
                #print(item)

            elif line.startswith("_item.category_id"):

                print(shlex.split(line), line_ix)
                category = "_" + shlex.split(line)[1]
                if category not in data:
                    data[category] = {}
                if item not in data[category]:
                    data[category][item] = {}
                
                #print(category)
                curr_category = category

            elif line.startswith("_item_type.code"):
                data_type = shlex.split(line)[1]
                if data_type == "float":
                    data_type = float
                elif data_type == "int":
                    data_type = int
                elif data_type == "str":
                    data_type = str
                
                data[curr_category][\
                    curr_item]["type"] = data_type

        return data

    def _parse_line_old(self,
                    line,
                    line_ix):
        """Parse a single line of the mmCIF file and yield the
        line's content as a list of single tokens.

        Parameters
        ----------
        line : ``str``
            The line.

        line_ix : ``int``
            The line's index.

        Yields
        ------
        token : ``str``
            The tokens that make up the line, one by one.
        """

        # Set a flag to indicate whether the current character parsed
        # is in a distinc token
        is_in_token = False

        # Set a variable that will hold, in case quotes are
        # encountered, what the open character is
        open_quote_char = None

        # Initialize the starting index (in the line) of the current
        # token to 0
        token_start = 0

        # For each character in the line, and its associated index
        for char_index, char in enumerate(line):


            #------------ The character is a white space -------------#


            # If the character is a white space
            if char in self.WHITESPACE_CHARACTERS:

                # If the character is inside a token and not inside
                # a quote
                if is_in_token and open_quote_char is None:

                    # Set the flag for whether the current character
                    # is inside a token to False (since it is a
                    # white space outside quotes)
                    is_in_token = False

                    # Add the previous token (which started at the
                    # saved starting index and ends at the current
                    # index) to the line's content
                    yield line[token_start:char_index]


            #--------------- The character is a quote ----------------#


            # If the character is a quote
            elif char in self.QUOTE_CHARACTERS:

                # Set a condition for whether we have reached
                # a closing quote matching the one that we have
                # open
                is_quote_closed = (char == open_quote_char)

                # Set a condition for whether the current
                # character is the last of the line
                is_last_character = \
                    (char_index + 1) == len(line)

                # Set a condition for whether the next character
                # is a white space (if we are in the middle of
                # the line)
                is_next_character_space = \
                    (char_index + 1) < len(line) and \
                    line[char_index + 1] in self.WHITESPACE_CHARACTERS

                # If we have no open quote and we are not inside
                # a token (= we are at the start of a quoted string)
                if not is_in_token and open_quote_char is None:

                    # Mark the quote as open
                    open_quote_char = char

                    # Set the flag for whether the character is
                    # inside a token to True
                    is_in_token = True

                    # Set the starting index of the new token
                    # to the one after the current one
                    token_start = char_index + 1

                # If the first condition is met, and either of the
                # following two are (= we reached the end of a
                # quote-enclosed token)
                elif is_quote_closed and \
                     (is_last_character or is_next_character_space):

                    # Mark the quote as closed
                    open_quote_char = None

                    # Set the flag for whether the character is
                    # inside a token to False
                    is_in_token = False

                    # Add the quotes string to the line's content
                    yield line[token_start:char_index]


            #----------- The character starts a new token ------------#


            # If the current character marks the start of a new
            # token
            elif not is_in_token:

                # Set the flag for whether the character is
                # inside a token to True
                is_in_token = True

                # Set the start index of a new token to the
                # current index
                token_start = char_index

        # If the end of the line is reached while still inside
        # a token (= the line ends with a non-white space and non-
        # quote character, because we deal with those cases above)
        if is_in_token:

            # Add the token to the line's content
            yield line[token_start:]

        # If the line ends while a quote is still open
        if open_quote_char:

            # Raise an error
            errstr = \
                f"Line {line_ix+1} ended while a quote " \
                "was open."
            raise ValueError(errstr)

load_mmcif_dict("/Users/lcb518/Documents/work/pdbcraft/pdbcraft/data/mmcif_pdbx_v50.dic")


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


    def _get_unique(self,
                    items_type,
                    attribute = None):
        """Get the unique values for items/attributes.

        Parameters
        ----------
        items_type : ``str``
            The type of items to get the last of (model, chain,
            segment, residue, atom, or atom attribute).

        attribute : ``str``
            If ``items_type`` is ``"atom_attribute"``, the name
            of the attribute whose unique values need to be found.
        """


        #-------------------- Define the recursion -------------------#


        def recursive_step(struct,
                           target_depth,
                           attribute,
                           current_depth = 1,
                           unique_values = set()):

            # If the current depth is the target depth
            if current_depth == target_depth:

                # If we need to find the values for a specific
                # attribute
                if attribute is not None:

                    # Add the values of the attribute at the
                    # current depth to the set of unique values
                    unique_values.update(struct[attribute])

                # Otherwise
                else:

                    # Add the items at the current depth to the
                    # set of unique values
                    unique_values.update(struct.keys())

            # Otherwise
            else:

                # If the current structure is a dictionary
                if isinstance(struct, dict):

                    # For each sub-structure
                    for sub_struct in struct.values():

                        # Recurse
                        recursive_step(\
                            struct = sub_struct,
                            target_depth = target_depth,
                            attribute = attribute,
                            current_depth = current_depth + 1,
                            unique_values = unique_values)

            # Return the unique values gathered while recursing,
            # sorted
            return sorted(list(unique_values))


        #------------------- Get the unique values -------------------#


        # Get the depth corresponding to where the items we want
        # to take the last of are stored in the structure
        target_depth = self._ITEM2DEPTH[items_type]

        # Return the result of the recursion, sorted
        return recursive_step(struct = self.atom_data,
                              target_depth = target_depth,
                              attribute = attribute)