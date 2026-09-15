#!/usr/bin/env python3
"""CDS synthetic validation: generates circular/random datasets and evaluates discrimination.
"""
import numpy as np,json,argparse
from cds import circularity_detection_score
from sklearn.metrics import roc_auc_score,confusion_matrix

def gen_dataset(n_samples=200,n_genes=5000,circular=True,seed=42):
 rng=np.random.RandomState(seed)
 X=rng.lognormal(4,1,(n_samples,n_genes))
 genes=[f"GENE_{i}" for i in range(n_genes)]
 genes[0],genes[1]="GZMA","PRF1"
 if circular:
  y=(X[:,0]+X[:,1]>np.median(X[:,0]+X[:,1])).astype(int)
 else:
  y=rng.binomial(1,0.5,n_samples)
 return X,y,genes,['GZMA','PRF1']

def run_validation(n_datasets=100):
 scores,labels=[],[]
 for i in range(n_datasets):
  circ=(i<n_datasets//2)
  X,y,genes,ep=gen_dataset(circular=circ,seed=i+42)
  r=circularity_detection_score(X,y,ep,genes)
  scores.append(r["CDS"])
  labels.append(1 if circ else 0)
 scores,labels=np.array(scores),np.array(labels)
 auroc=roc_auc_score(labels,scores)
 pred=(scores>=0.70).astype(int)
 tn,fp,fn,tp=confusion_matrix(labels,pred).ravel()
 sens=tp/(tp+fn)if(tp+fn)>0 else 0
 spec=tn/(tn+fp)if(tn+fp)>0 else 0
 return{"AUROC":round(float(auroc),4),"sensitivity":round(float(sens),4),"specificity":round(float(spec),4),"mean_cds_circular":round(float(scores[labels==1].mean()),4),"mean_cds_random":round(float(scores[labels==0].mean()),4)}

def main():
 p=argparse.ArgumentParser(description="Validate CDS")
 p.add_argument("--n-datasets",type=int,default=100)
 p.add_argument("--output",type=str,default=None)
 a=p.parse_args()
 results=run_validation(a.n_datasets)
 print(f"AUROC={results['AUROC']:.4f} Sens={results['sensitivity']:.4f} Spec={results['specificity']:.4f}")
 if a.output:
  json.dump(results,open(a.output,"w"),indent=2)
if __name__=="__main__":
 exit(main())
