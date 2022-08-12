class Residue:

    def __init__(self, parent, start, end):
        self.parent = parent
        self.start = start
        self.end = end

    def number(self):
        return int(self.parent.parent[start][22:26])

    def icode(self):
        return self.parent.parent[start][26]
