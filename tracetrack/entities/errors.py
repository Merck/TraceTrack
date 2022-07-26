class AlignmentError(Exception):
    def __init__(self, message):
        self.message = message


class ReferenceError(Exception):
    def __init__(self, message):
        self.message = message


class ParsingError(Exception):
    def __init__(self, message):
        self.message = message


class UnresolvedConflictError(Exception):
    def __init__(self, message):
        self.message = message


