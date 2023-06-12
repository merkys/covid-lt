class Atom:

    def __init__(self, parent, line):
        self.parent = parent
        self.line = line
        self.parse_line() 

    def element(self, element=None):
        old_element = self.parent.parent.parent.content[self.line][76:78].strip()
        if element:
            element = element.strip()
            if len(element) == 1:
                element = ' ' + element
            elif len(element) > 2:
                element = element[0:2]
            self.parent.parent.parent.content[self.line] = self.parent.parent.parent.content[self.line][:76] + element + self.parent.parent.parent.content[self.line][78:]
        return(old_element)
    
    def parse_line(self):
        """
        Parse an ATOM/HETATM line from a PDB file and return a dictionary of the fields.

        Parameters
        ----------
        line : str
            The ATOM/HETATM line to parse.

        Returns
        -------
        dict
            A dictionary containing the parsed fields, with the keys being the field names.
        """
        line = self.parent.parent.parent.content[self.line]
        fields = {
            'record_name': line[0:6].strip(),
            'atom_number': int(line[6:11]),
            'atom_name': line[12:16].strip(),
            'alt_loc': line[16],
            'res_name': line[17:20].strip(),
            'chain_id': line[21],
            'res_seq': int(line[22:26]),
            'i_code': line[26],
            'x': float(line[30:38]),
            'y': float(line[38:46]),
            'z': float(line[46:54]),
            'occupancy': float(line[54:60]),
            'temp_factor': float(line[60:66]),
            'element': line[76:78].strip(),
            'charge': line[78:80].strip()
        }
        self.fields = fields
    
    def coords(self):
        return(
                [
                self.fields['x'],
                self.fields['y'],
                self.fields['z'],
                ]
        )

    def name(self):
        return self.fields["atom_name"]

    @property
    def number(self):
        return self.fields["atom_number"]

    @property
    def line_from_fields(self):
        """
        Generate an ATOM/HETATM line for a PDB file from a dictionary of the fields.

        Parameters
        ----------
        fields : dict
            A dictionary containing the fields for the ATOM/HETATM line, with the keys being the field names.

        Returns
        -------
        str
            The generated ATOM/HETATM line.
        """
        fields = self.fields
        format_string = "{:<6}{:>5}  {:<3}{:1}{:<3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:>2}"
        args = (
            fields['record_name'],
            fields['atom_number'],
            fields['atom_name'],
            fields['alt_loc'],
            fields['res_name'],
            fields['chain_id'],
            fields['res_seq'],
            fields['i_code'],
            fields['x'],
            fields['y'],
            fields['z'],
            fields['occupancy'],
            fields['temp_factor'],
            fields['element'].strip(),
            fields['charge'].strip()
        )
        return format_string.format(*args)


    
    # Parent links

    @property
    def residue(self):
        return self.parent

    @property
    def chain(self):
        return self.parent.parent

    @property
    def file(self):
        return self.parent.parent.parent

    def within(self, distance):
        cKDTree, atoms = self.parent.parent.parent._get_cKDTree()
        return [atoms[x] for x in cKDTree.query_ball_point( self.coords(), distance )]
