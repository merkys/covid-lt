class PDBFile:

    def __init__(self, fileobject):
        if isinstance(fileobject, str):
            self.file = open(fileobject, 'r')
        elif isinstance(fileobject, file):
            self.file = fileobject
        else:
            ValueError("Cannot process {}" % type(fileobject))
        self.content = self.file.readlines()

    def get(self, keyword):
        while len(keyword) < 6:
            keyword = keyword + ' '
        return filter(lambda x: x.startswith(keyword), self.content)
