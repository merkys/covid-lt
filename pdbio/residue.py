from pdbio.atom import Atom

class Residue:

    def __init__(self, parent, start, end):
        self.parent = parent
        self.start = start
        self.end = end

    def __iter__(self):
        self.iter_atom_line = self.start
        return self

    def __next__(self):
        if self.iter_atom_line > self.end:
            raise StopIteration()
        while not self.parent.parent.content[self.iter_atom_line].startswith('ATOM  ') and self.iter_atom_line < self.end:
            self.iter_atom_line += 1
        if not self.parent.parent.content[self.iter_atom_line].startswith('ATOM  '):
            raise StopIteration()
        self.iter_atom_line += 1
        return Atom(self, self.iter_atom_line-1)

    def icode(self):
        return self.parent.parent.content[self.start][26]

    def number(self):
        return int(self.parent.parent.content[self.start][22:26])

    def resname(self):
        return self.parent.parent.content[self.start][17:20]

    def delete(self):
        # FIXME: Other fields should as well be updated to reflect the removal of this residue
        # FIXME: There may be other fields between start and end
        self.parent.parent.content = self.parent.parent.content[:self.start] + self.parent.parent.content[self.end+1]
