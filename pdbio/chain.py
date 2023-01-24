from Bio.Data.IUPACData import protein_letters_3to1
from pdbio.residue import Residue
from warnings import warn

class Chain:

    anarci_chain_types = { 'H': 'H', 'K': 'L', 'L': 'L' }

    def __init__(self, parent, name):
        self.parent = parent
        self.name = name

    def __iter__(self):
        self.iter_residue_line = 0
        return self

    def __next__(self):
        while self.iter_residue_line < len(self.parent.content) and (not self.parent.content[self.iter_residue_line].startswith('ATOM  ') or self.parent.content[self.iter_residue_line][21] != self.name):
            self.iter_residue_line += 1
        if self.iter_residue_line + 1 >= len(self.parent.content): # File end has been reached
            raise StopIteration()
        start, end = self.iter_residue_line, self.iter_residue_line
        atom = self.parent.content[start]
        this_chain, this_number, this_icode = atom[21], int(atom[22:26]), atom[26]
        while self.iter_residue_line + 1 < len(self.parent.content):
            self.iter_residue_line += 1
            if not self.parent.content[self.iter_residue_line].startswith('ATOM  '):
                continue
            atom = self.parent.content[self.iter_residue_line]
            chain, number, icode = atom[21], int(atom[22:26]), atom[26]
            if [this_chain, this_number, this_icode] == [chain, number, icode]:
                end = self.iter_residue_line
            else:
                break
        return Residue(self, start, end)

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
            sequence = sequence + protein_letters_3to1[residue.resname().capitalize()]
        return sequence

    def sequence_seqres(self):
        sequence = None
        for line in self.parent.get('SEQRES'):
            if line[11] != self.name:
                continue
            if sequence is None:
                sequence = ''
            for i in range(0,13):
                residue = line[19+i*4:22+i*4]
                if residue != '   ':
                    sequence = sequence + protein_letters_3to1[residue.capitalize()]
        return sequence
