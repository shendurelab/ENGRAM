import nupack as nup
import numpy as np
import pandas as pd

model1 = nup.Model(material='rna', celsius=37,sodium=0.5, magnesium=0.0)

peg = 'gcGGCCCAGACUGAGCACGUGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGGACCGAGUCGGUCCUCUGCCAUCANNNNNCGUGCUCAGUCUGACUGCCGUAUAggca'
nts = ['A','U','C','G']
bcs = [a+b+c+d+e for a in nts for b in nts for c in nts for d in nts for e in nts]
table = []
for bc in bcs:
    new_peg = nup.Strand(peg.replace('NNNNN',bc),name='new_peg')
    set1 = nup.ComplexSet(strands=[new_peg])
    result = nup.complex_analysis(complexes=set1, model=model1, compute=['pairs', 'mfe'])
    table.append([bc,result['(new_peg)'].mfe[0].energy]+\
                  list(np.diag(np.rot90(result['(new_peg)'].pairs.to_array()[19:24,108:113])))) # probability of barcode pairing to spacer

pd.DataFrame(table,columns=['barcode','MFE','P1','P2','P3','P4','P5',]).to_csv("your path to save" + '/pegRNA_Free_G.csv',index=False,sep='\t')