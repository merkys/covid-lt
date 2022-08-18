from io import TextIOWrapper
from pdbio.chain import Chain
import sys

class PDBFile:

    def __init__(self, fileobject):
        if isinstance(fileobject, str):
            self.file = open(fileobject, 'r')
        elif isinstance(fileobject, TextIOWrapper):
            self.file = fileobject
        else:
            ValueError('Cannot process {}' % type(fileobject))
        self.content = self.file.readlines()

    def __iter__(self):
        self.iterchain = -1
        return self

    def __next__(self):
        self.iterchain += 1
        if self.iterchain < len(self.chains()):
            return self.chain(sorted(self.chains())[self.iterchain])
        else:
            raise StopIteration()

    def get(self, keyword):
        while len(keyword) < 6:
            keyword = keyword + ' '
        return filter(lambda x: x.startswith(keyword), self.content)

    def chain(self, chain):
        return Chain(self, chain)

    def chains(self):
        return set([x[21] for x in self.get('ATOM')])

    def sequence(self, chain):
        return self.chain(chain).sequence()

    def renumber(self, func):
        renumberings = {
            'ATOM  ': [ [ 21, 22, 25, 26 ] ],
            'ANISOU': [ [ 21, 22, 25, 26 ] ],
            'CISPEP': [ [ 15, 17, 20, 21 ], [ 29, 31, 34, 35 ] ],
            'HELIX ': [ [ 19, 21, 24, 25 ], [ 31, 33, 36, 37 ] ],
            'LINK  ': [ [ 21, 22, 25, 26 ], [ 51, 52, 55, 56 ] ],
            'SHEET ': [ [ 21, 22, 25, 26 ], [ 32, 33, 36, 37 ], [ 49, 50, 53, 54 ], [ 64, 65, 68, 69 ] ],
            'SITE  ': [ [ 22, 23, 26, 27 ], [ 33, 34, 37, 38 ], [ 44, 45, 48, 49 ], [ 55, 56, 59, 60 ] ],
            'SSBOND': [ [ 15, 17, 20, 21 ], [ 29, 31, 34, 35 ] ],
        }

        for i, atom in enumerate(self.content):
            if atom[0:6] not in renumberings:
                continue
            for renumbering in renumberings[atom[0:6]]:
                number = None
                try:
                    number = int(atom[renumbering[1]:renumbering[2]+1])
                except ValueError: # blank residue found, possible in HELIX and other fields
                    continue
                chain, number, icode = func(atom[renumbering[0]], number, atom[renumbering[3]])
                if chain is not None:
                    self.content[i] = self.content[i][0:renumbering[0]] + chain + self.content[i][renumbering[0]+1:]
                if number is not None:
                    number = str(number)
                    while len(number) < 4:
                        number = ' ' + number
                    self.content[i] = self.content[i][0:renumbering[1]] + number + self.content[i][renumbering[2]+1:]
                if icode is not None:
                    self.content[i] = self.content[i][0:renumbering[3]] + icode + self.content[i][renumbering[3]+1:]

    def __str__(self):
        return ''.join(self.content)
