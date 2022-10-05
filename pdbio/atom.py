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
        return(old_element)

    def name(self):
        return self.parent.parent.parent.content[self.line][12:16].strip()
