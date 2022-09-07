class Residue:

    def __init__(self, parent, start, end):
        self.parent = parent
        self.start = start
        self.end = end

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
