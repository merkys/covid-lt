from Bio.Data.IUPACData import protein_letters_3to1_extended
from pdbio.residue import Residue
from warnings import warn

class Chain:

    anarci_chain_types = { 'H': 'H', 'K': 'L', 'L': 'L' }

    def __init__(self, parent, name):
        self.parent = parent
        self.name = name

    def __iter__(self):
        iter_residue_line = -1
        while iter_residue_line < len(self.parent.content)-1:
            iter_residue_line += 1
            if not self.parent.content[iter_residue_line].startswith('ATOM  ') or self.parent.content[iter_residue_line][21] != self.name:
                continue
            start, end = iter_residue_line, iter_residue_line
            atom = self.parent.content[start]
            this_chain, this_number, this_icode = atom[21], int(atom[22:26]), atom[26] 
            while iter_residue_line + 1 < len(self.parent.content):
                iter_residue_line += 1
                if not self.parent.content[iter_residue_line].startswith('ATOM  '):
                    continue
                atom = self.parent.content[iter_residue_line]
                chain, number, icode = atom[21], int(atom[22:26]), atom[26]
                if [this_chain, this_number, this_icode] == [chain, number, icode]:
                    end = iter_residue_line
                else:
                    break
            yield Residue(self, start, end)

    def __len__(self):
        return len([residue for residue in self])

    def antibody_numbering(self):
        from anarci import run_anarci
        _, numbered, details, _ = run_anarci([(self.name, self.sequence())], scheme='chothia', allow=set(self.anarci_chain_types.keys()))
        numbered = numbered[0]
        details = details[0]
        if numbered is None:
            return None
        if len(numbered) > 1:
            warn('more than one H or L fragment was found in chain {}, using the first'.format(self.name))
        return numbered[0]

    def antibody_type(self):
        from anarci import run_anarci
        _, numbered, details, _ = run_anarci([(self.name, self.sequence())], scheme='chothia', allow=set(self.anarci_chain_types.keys()))
        numbered = numbered[0]
        details = details[0]
        if numbered is None:
            return None
        if len(numbered) > 1:
            warn('more than one H or L fragment was found in chain {}, using the first'.format(self.name))
        return self.anarci_chain_types[details[0]['chain_type']]

    def is_contiguous(self):
        last = None
        for residue in self:
            if last is not None and residue.number() - last != 1:
                return False
            last = residue.number()
        return True

    def rename(self, name):
        self.parent.rename_chains({self.name: name})

    def renumber(self, func):
        def _renumber_this_chain(chain, number, icode):
            if chain != self.name:
                return None, None, None
            else:
                return None, *func(number, icode)
        self.parent.renumber(_renumber_this_chain)

    def sequence(self):
        sequence = self.sequence_seqres()
        if sequence:
            return sequence
        else:
            return self.sequence_atom()

    def sequence_atom(self):
        sequence = None
        for residue in self:
            if sequence is None:
                sequence = ''
            sequence = sequence + protein_letters_3to1_extended[residue.resname().capitalize()]
        return sequence

    def sequence_seqres(self, return_3L = False):
        """ 
        Returns single letter coded sequences as string  or single leter 
        sequence as string and a list of three letter representations
        """
            
        sequence = None
        sequence3L = []
        out=""
        
        for line in self.parent.get('SEQRES'):
            if line[11] != self.name:
                continue
            if sequence is None:
                sequence = ''
            for i in range(0,13):
                residue = line[19+i*4:22+i*4]
                if residue != '   ' and residue != '':
                    resU = residue.capitalize()
                    sequence3L.append(resU)
                    if (resU in protein_letters_3to1_extended.keys()):
                        sequence = sequence + protein_letters_3to1_extended[resU]
                    else:
                        sequence = sequence + "X"
                            
        if return_3L:
            return (sequence,sequence3L)
        else:
            return sequence

    def within(self, distance, result_class='Chain', prefilter=None):
        cKDTree, atoms = self.parent._get_cKDTree(prefilter=prefilter)

        coordinates = []
        for residue in self:
            for atom in residue:
                if not prefilter or prefilter(atom):
                    coordinates.append( atom.coords() )

        merged = set()
        for hits in cKDTree.query_ball_point( coordinates, distance ):
            merged.update(set(hits))
        result = [atoms[x] for x in merged]
        result = filter(lambda x: x.parent.parent.name != self.name, result)

        if result_class == 'Atom':
            return list(result)
        result = set(atom.parent for atom in result)
        if result_class == 'Residue':
            return list(result)
        result = set(residue.parent for residue in result)
        return list(result)
