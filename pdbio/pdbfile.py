from Bio.Data.SCOPData import protein_letters_3to1
from io import TextIOWrapper

class PDBFile:

    def __init__(self, fileobject):
        if isinstance(fileobject, str):
            self.file = open(fileobject, 'r')
        elif isinstance(fileobject, TextIOWrapper):
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

    def renumber(self, func):
        for i, atom in enumerate(self.content):
            if not atom.startswith('ATOM  '):
                continue
            chain, number, icode = func(atom[21], int(atom[22:26]), atom[26])
            if chain is not None:
                self.content[i] = self.content[i][0:21] + chain + self.content[i][22:]
            if number is not None:
                number = str(number)
                while len(number) < 4:
                    number = ' ' + number
                self.content[i] = self.content[i][0:22] + number + self.content[i][26:]
            if icode is not None:
                self.content[i] = self.content[i][0:26] + icode + self.content[i][27:]

    def __str__(self):
        return ''.join(self.content)
