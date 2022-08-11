from Bio.Data.SCOPData import protein_letters_3to1

class PDBFile:

    def __init__(self, fileobject):
        if isinstance(fileobject, str):
            self.file = open(fileobject, 'r')
        elif isinstance(fileobject, file):
            self.file = fileobject
        else:
            ValueError('Cannot process {}' % type(fileobject))
        self.content = self.file.readlines()

    def get(self, keyword):
        while len(keyword) < 6:
            keyword = keyword + ' '
        return filter(lambda x: x.startswith(keyword), self.content)

    def chains(self):
        return set([x[21] for x in self.get('ATOM')])

    def sequence(self, chain):
        sequence = None
        for line in self.get('SEQRES'):
            if line[11] != chain:
                continue
            if sequence is None:
                sequence = ''
            for i in range(0,13):
                residue = line[19+i*4:22+i*4]
                if residue != '   ':
                    sequence = sequence + protein_letters_3to1[residue]
        return sequence
