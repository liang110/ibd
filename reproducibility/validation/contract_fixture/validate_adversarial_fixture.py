from pathlib import Path
import pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss

ROOT=Path(__file__).parent
pred=pd.read_csv(ROOT/"adversarial_predictions.csv")
expected=pd.read_csv(ROOT/"adversarial_expected_outputs.csv")
y=pred.outcome_ibd.to_numpy(int); runs=[pred.run_1.to_numpy(float),pred.run_2.to_numpy(float)]
def metrics(y,p): return {"AUROC":roc_auc_score(y,p),"AP":average_precision_score(y,p),"Brier":brier_score_loss(y,p)}
def part(p):
 z=pred.assign(prob=p).groupby("participant_id",sort=False).agg(outcome_ibd=("outcome_ibd","first"),prob=("prob","mean"))
 return z.outcome_ibd.to_numpy(),z.prob.to_numpy()
vals=[]
for r in runs: vals.append({"A":metrics(y,r),"C":metrics(*part(r))})
p=np.mean(runs,axis=0)
calc={"A":{m:np.mean([v["A"][m] for v in vals]) for m in ["AUROC","AP","Brier"]},"B":metrics(y,p),"C":{m:np.mean([v["C"][m] for v in vals]) for m in ["AUROC","AP","Brier"]},"D":metrics(*part(p))}
for _,row in expected.iterrows():
 got=calc[row.object][row.metric]
 assert abs(got-row.value)<1e-12,(row.object,row.metric,got,row.value)
print("adversarial fixture PASS",len(expected),"values")
