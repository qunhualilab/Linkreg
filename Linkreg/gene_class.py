import numpy as np
import copy



class gene_class(object):
    def __init__(self, cCRE=None, cell_num=None, expression=None, tss=None, gene_loc=None, cCRE_loc=None, cCRE_dist=None, gene_name=None, chromosome=None, proximal_trsd=None, strain=None):
        # cCRE format: [cCRE1_c1, cCRE1_c2, ...]
        self.gene_name = gene_name
        self.cCRE = cCRE
        self.ideas_num = len(cCRE[0])
        self.cell_num = cell_num
        self.tss = tss
        self.gene_loc = gene_loc
        self.cCRE_loc = cCRE_loc
        self.cCRE_num = len(self.cCRE)//self.cell_num
        if cCRE_dist is None:
            self.cCRE_dist = np.ones(self.cCRE_num)
            self.proximal_trsd = 0
        else:
            self.cCRE_dist = cCRE_dist
            self.proximal_trsd = proximal_trsd
            self.original_cCRE_dist = copy.deepcopy(cCRE_dist)
        self.expression = expression
        self.cCRE_status = np.ones(self.cCRE_num)
        self.inclusion_status = np.ones(self.cCRE_num) # length stays the same
        self.chromosome = chromosome
        self.strain = strain
    
    def pre_select_cCRE(self, cCRE_status):
        self.inclusion_status = cCRE_status
        self.cCRE_num = np.sum(cCRE_status)
        self.cCRE_status = np.ones(np.sum(cCRE_status).astype(int))
        self.cCRE = self.cCRE[np.repeat(cCRE_status, self.cell_num).astype(bool)]
        self.cCRE_dist = self.cCRE_dist[cCRE_status.astype(bool)]
        idx = np.where(self.cCRE_dist<=self.proximal_trsd)[0]
        if len(idx) == 0:
            self.proximal = np.array([np.zeros(self.ideas_num) for i in range(self.cell_num)])
        else:
            idx = np.array([idx[i]*self.cell_num for i in range(len(idx))])
            self.proximal = np.array([np.mean(self.cCRE[idx+i], axis=0) for i in range(self.cell_num)])
    
    def get_distal(self):
        idx = np.where(self.cCRE_status==1)[0]
        if len(idx) == 0:
            self.distal = np.array([np.zeros(self.ideas_num) for i in range(self.cell_num)])
        else:
            idx = np.array([idx[i]*self.cell_num for i in range(len(idx))])
            self.distal = np.array([np.mean(self.cCRE[idx+i], axis=0) for i in range(self.cell_num)])
        return self.distal
        
    def get_effective_cCRE_num(self):
        return np.sum(self.cCRE_status)
    
    def refine_cCRE(self, cCRE_status, iter_num):
        idx = np.where(cCRE_status!=self.cCRE_status)[0]
        if 100//iter_num+1 < len(idx):
            idx = np.random.choice(idx, 100//iter_num+1, replace=False)
            self.cCRE_status[idx] = 1-self.cCRE_status[idx]
        else:
            self.cCRE_status = cCRE_status
    
    def add_eRP_scores(self, eRP_scores):
        # eRP_scores format: [[cCRE1, cCRE2, ...], ...]
        self.final_cCRE = {'eRP_scores': eRP_scores}
        self.final_cCRE['cCRE_idx'] = np.where(self.inclusion_status==1)[0][np.where(self.cCRE_status==1)[0]]
        if self.cCRE_loc is not None:
            self.final_cCRE['cCRE_loc'] = self.cCRE_loc[self.inclusion_status==1][self.cCRE_status==1]
        if self.gene_name is not None:
            self.final_cCRE['gene_name'] = self.gene_name
        return self.final_cCRE
    
    def get_effective_cCRE(self, num):
        self.effective_cCRE = {}
        effectives = []
        for i in range(self.cell_num):
            effectives.append(self.final_cCRE['cCRE_idx'][np.argsort(self.final_cCRE['eRP_scores'][i])[-num:]])
        self.effective_cCRE['effective_cCRE_idx'] = effectives
        if self.cCRE_loc is not None:
            self.effective_cCRE['cCRE_loc'] = [self.cCRE_loc[effectives[i]] for i in range(self.cell_num)]
        if self.gene_name is not None:
            self.effective_cCRE['gene_name'] = self.gene_name
        return self.effective_cCRE
    
    def get_unique_effective_cCRE_idx(self, num):
        unique_effective_cCRE, counts = np.unique(self.effective_cCRE['effective_cCRE_idx'], return_counts=True)
        top_effective_cCRE = unique_effective_cCRE[np.argsort(counts)[-num:]]
        return unique_effective_cCRE, top_effective_cCRE
    
    def set_proximal_trsd(self, proximal_trsd):
        self.proximal_trsd = proximal_trsd