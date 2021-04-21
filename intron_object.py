class intron:
    def __init__(self, scaffold_id, donor, acceptor, strand):
        self.scaffold_id = scaffold_id
        self.donor = donor - 1
        self.acceptor = acceptor - 1
        self.strand = strand