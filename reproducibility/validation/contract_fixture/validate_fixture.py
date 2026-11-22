from pathlib import Path
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss

root=Path(__file__).parent
f=pd.read_csv(root/"predictions.csv")
if f.groupby("participant_id").outcome_ibd.nunique().gt(1).any(): raise SystemExit("participant label conflict")
part=f.groupby("participant_id").agg(outcome_ibd=("outcome_ibd","first"),run_1=("run_1","mean"),run_2=("run_2","mean")).reset_index()
def metrics(y,p): return {"AUROC":roc_auc_score(y,p),"AP":average_precision_score(y,p),"Brier":brier_score_loss(y,p)}
out={}
out["A"]={m:sum(metrics(f.outcome_ibd,f[r])[m] for r in ["run_1","run_2"])/2 for m in ["AUROC","AP","Brier"]}
out["B"]=metrics(f.outcome_ibd,f[["run_1","run_2"]].mean(axis=1))
out["C"]={m:sum(metrics(part.outcome_ibd,part[r])[m] for r in ["run_1","run_2"])/2 for m in ["AUROC","AP","Brier"]}
out["D"]=metrics(part.outcome_ibd,part[["run_1","run_2"]].mean(axis=1))
e=pd.read_csv(root/"expected_outputs.csv")
for _,r in e.iterrows():
    if abs(out[r.object][r.metric]-r.value)>1e-12: raise SystemExit(f"mismatch {r.object} {r.metric}")
print("fixture PASS")
