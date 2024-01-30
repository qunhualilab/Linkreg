import numpy as np
import pandas as pd
from gene_class import gene_class

def find_ends_index(locations, tss, distance):
    start_index = np.searchsorted(np.array(locations[:, 0]), tss - distance)
    end_index = np.searchsorted(np.array(locations[:, 1]), tss + distance, side='right')
    return start_index, end_index
def find_values_within_distance(locations, tracks, loc, distance):
    start_index, end_index = find_ends_index(locations, loc, distance)
    values = np.array([track[start_index:end_index].reshape(-1) for track in tracks]).transpose(1, 0)
    return values, locations[start_index:end_index]
def process(expression_file, track_files, distance=500000):
    expression = np.array(pd.read_csv(expression_file, delimiter='\t'))
    tracks = [np.array(pd.read_csv(track_file, delimiter='\t')) for track_file in track_files]
    cell_num = len(expression[0])-5
    genes = []
    locations = np.array(tracks[0][:, :3])
    for chr_str in ['chr'+str(i) for i in range(1, 23)]+['chrX']:
        expression_ch = expression[expression[:, 0]==chr_str]
        tracks_ch = [track[track[:, 0]==chr_str, 1:] for track in tracks]
        locations = np.array(tracks_ch[0][:, :2])
        tracks_ch = [track[:, 2:] for track in tracks_ch]
        for row in expression_ch:
            gene = gene_class(ID=row[3], cell_num=cell_num, expression=list(row[5:]), gene_body=list(row[1:3]), strand=row[4], chromosome=row[0])
            # cCRE, cCRE_loc
            cCRE, cCRE_loc = find_values_within_distance(locations, tracks_ch, gene.tss, distance)
            gene.cCRE_loc = cCRE_loc
            gene.cCRE = cCRE
            genes.append(gene)
    return genes
