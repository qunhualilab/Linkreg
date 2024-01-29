class gene_class(object):
    def __init__(self, cCRE=[], ID=None, cell_num=None, expression=None, tss=None, gene_body=None, strand=None, cCRE_loc=None, gene_symbol=None, chromosome=None):
        # cCRE format: [cCRE1_c1, cCRE1_c2, ...]
        self.gene_symbol = gene_symbol
        self.cell_num = cell_num
        self.cCRE = cCRE
        self.tss = gene_body[0] if strand == '+' else gene_body[1]
        self.gene_body = gene_body
        self.cCRE_loc = cCRE_loc
        self.cCRE_num = len(cCRE)//cell_num if cCRE != [] else 0
        self.expression = expression
        self.chromosome = chromosome
        self.strand = strand
        self.ID = ID
    
    @property
    def cCRE(self):
        return self._cCRE

    @cCRE.setter
    def cCRE(self, cCRE):
        self._cCRE = cCRE
        self.cCRE_num = len(cCRE)//self.cell_num
