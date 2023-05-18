class Atom:

    def __init__(self, parent, line):
        self.parent = parent
        self.line = line

    def element(self, element=None):
        old_element = self.parent.parent.parent.content[self.line][76:78].strip()
        if element:
            element = element.strip()
            if len(element) == 1:
                element = ' ' + element
            elif len(element) > 2:
                element = element[0:2]
            self.parent.parent.parent.content[self.line] = self.parent.parent.parent.content[self.line][:76] + element + self.parent.parent.parent.content[self.line][78:]
        return old_element

    def coords(self):
        return [float(self.parent.parent.parent.content[self.line][30:38]),
                float(self.parent.parent.parent.content[self.line][38:46]),
                float(self.parent.parent.parent.content[self.line][46:54])]

    def name(self, name=None):
        old_name = self.parent.parent.parent.content[self.line][12:16].strip()
        if name:
            if len(name) > 4:
                raise ValueError("atom name cannot be longer than 4 characters")
            if len(name) < 4:
                name = ' ' + name
            while len(name) < 4:
                name = name + ' '
            self.parent.parent.parent.content[self.line] = self.parent.parent.parent.content[self.line][:12] + element + self.parent.parent.parent.content[self.line][16:]
        return old_name

    @property
    def number(self):
        return int(self.parent.parent.parent.content[self.line][6:11])

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
