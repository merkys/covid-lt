from pdbio.residue import Residue

class Chain:

    def __init__(self, parent, name):
        self.parent = parent
        self.name = name

    def __iter__(self):
        self.iter_residue_line = 0
        return self

    def __next__(self):
        while self.iter_residue_line < len(self.parent.content) and not self.parent.content[self.iter_residue_line].startswith('ATOM  ') and self.parent.content[self.iter_residue_line][21] != self.name:
            self.iter_residue_line += 1
        if self.iter_residue_line == len(self.parent.content): # File end has been reached
            return None
        start, end = self.iter_residue_line, self.iter_residue_line
        atom = self.parent.content[start]
        this_chain, this_number, this_icode = atom[21], int(atom[22:26]), atom[26]
        while self.iter_residue_line < len(self.parent.content):
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

    def is_contiguous(self):
        last = None
        for residue in self:
            if last is not None:
                if residue.number() - last != 1:
                    return False
            else:
                last = residue.number()
        return True
