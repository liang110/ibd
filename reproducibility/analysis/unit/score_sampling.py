"""Evaluate how sampling weights and eligible specimens change post-fit IBD scores."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

QUANTILES = [0.025, 0.975]
SPLIT_METRICS = [roc_auc_score, average_precision_score, brier_score_loss]


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_scorer(bundle):
    path = bundle / "validation/replay_registered.py"
    spec = importlib.util.spec_from_file_location("registered_scorer", path)
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    return scorer


def test_error(checks, key, errors, tolerance=1e-9):
    maximum = float(np.abs(errors).max())
    checks[key] = max(checks.get(key, 0), maximum)
    if not np.isfinite(maximum) or maximum > tolerance:
        raise ValueError(f"{key}: maximum discrepancy {maximum} exceeds {tolerance}")


def prepare_panel(scorer, frame, weights):
    samples, people, groups, sw, pw, draws = scorer.direction_inputs(frame, weights)
    y = samples.outcome_ibd.to_numpy()
    py = people.outcome_ibd.to_numpy()
    counts = samples.groupby("bootstrap_cluster").sample_id.transform("size").to_numpy()
    return samples, people, groups, sw, pw, draws, y, py, sw / counts


def make_bridge(scorer, panel, prob, result, checks, baseline):
    samples, people, groups, sw, pw, draws, y, py, fw = panel
    per_fit = [scorer.score(y, run, fw) for run in prob]
    e = np.mean(per_fit, axis=0)
    # Recompute E independently with scikit-learn on the target and two stored nonuniform draws.
    for i in (0, 1, len(draws)):
        ref = np.mean([[metric(y, run, sample_weight=fw[i]) for metric in SPLIT_METRICS] for run in prob], axis=0)
        test_error(checks, "independent_E_full_and_two_draws", ref - e[:,i])
    a, c = result["A"], result["C"]
    variances = np.array([np.var(prob[:,g], axis=1, ddof=0).mean() for g in groups])
    weighted_variance = (pw * variances).sum(axis=1)/pw.sum(axis=1)
    test_error(checks, "within_person_brier_variance_all_draws", e[2]-c[2]-weighted_variance)
    test_error(checks, "telescoping_all_draws", (e-a)+(c-e)-(c-a))
    if len(samples) == len(people):
        test_error(checks, "singleton_bridge", e-a)
        test_error(checks, "singleton_A_C", c-a)
        test_error(checks, "singleton_B_D", result["D"]-result["B"])
    if baseline:
        test_error(checks, "constant_E_AUROC", e[0]-0.5)
        test_error(checks, "constant_E_AP", e[1]-(pw*py).sum(axis=1)/pw.sum(axis=1))
    return {"A":a,"E":e,"C":c,"E_minus_A":e-a,"C_minus_E":c-e,"C_minus_A":c-a}


def stats(values):
    lo, hi = np.quantile(values[1:], QUANTILES, method="linear")
    return float(values[0]), float(lo), float(hi)


def run(args):
    bundle = args.bundle.resolve()
    scorer = load_scorer(bundle)
    base = bundle / "evidence"
    files = {"predictions":base/"full_reproducibility_inputs/prediction_registry.csv",
             "weights":args.weights or base/"full_reproducibility_inputs/shared_bootstrap_weights.csv",
             "ESP_weights":base/"nielsen_esp_audit/esp_only_shared_weights.csv", "order":args.order}
    read = lambda f: pd.read_csv(f, float_precision="round_trip", keep_default_na=False)
    pred = read(files["predictions"])
    orders = read(files["order"])
    weights = read(files["weights"])
    esp_weights = read(files["ESP_weights"])
    if pred.duplicated(["direction","model","split","sample_id"]).any() or not pred.probability.between(0,1).all():
        raise ValueError("Invalid input predictions")
    if not pred.outcome_ibd.isin([0,1]).all():
        raise ValueError("Invalid outcome")
    args.output.mkdir(parents=True, exist_ok=True)
    bridge_rows, bridge_draws, windows, contrasts, structure = [], [], [], [], []
    standard_points, standard_draws, esp_points = [], [], []
    checks = {}
    for direction, whole in pred.groupby("direction"):
        target = direction.split("_to_")[1]
        if args.key_cases:
            choice = {"PRISM_to_HMP2":"shallow_mlp","HMP2_to_NielsenHB_2014":"logistic_l2"}.get(direction)
            if choice is None:
                continue
            whole = whole.loc[whole.model.eq(choice)]
        subsets = [("full",whole,weights)]
        if target == "NielsenHB_2014" and not args.key_cases:
            subsets.append(("ESP",whole.loc[whole.origin.eq("ESP")],esp_weights))
        for scope, frame, global_weights in subsets:
            draw_table = global_weights.loc[global_weights.direction.eq(direction)]
            if draw_table.empty:
                raise ValueError(f"No bootstrap weights for {direction}")
            order = orders.loc[orders.target.eq(target)].set_index("sample_id").retained_order
            panels = {"all": frame}
            if target != "PRISM":
                ordering = frame.sample_id.map(order)
                if ordering.isna().any():
                    raise ValueError(f"Missing acquisition order for {direction}")
                panels.update(first=frame.loc[ordering.eq(1)],first_two=frame.loc[ordering.le(2)])
            all_scores = {}
            for window, selected in panels.items():
                panel = prepare_panel(scorer, selected, draw_table)
                samples, people, groups, sw, pw, draws, y, py, fw = panel
                for label, cases in samples.groupby("outcome_ibd"):
                    sizes = cases.groupby("bootstrap_cluster").size()
                    structure.append(dict(direction=direction,scope=scope,target=target,window=window,outcome_ibd=label,
                                          n_samples=len(cases), n_participants=len(sizes),n_repeated=int((sizes>1).sum()),
                                          n1=int((sizes==1).sum()), n2=int((sizes==2).sum()), n3plus=int((sizes>=3).sum())))
                for model, block in selected.groupby("model"):
                    prob = block.pivot(index="split",columns="sample_id",values="probability").reindex(columns=samples.sample_id).to_numpy()
                    if prob.shape[0] != 10 or not np.isfinite(prob).all():
                        raise ValueError("Incomplete 10-run panel")
                    result = scorer.four_scores(prob, y, py, groups, sw, pw)
                    meta = dict(scope=scope,direction=direction,model=model,n_samples=len(y),n_participants=len(py),
                                cases_samples=int(y.sum()),cases_participants=int(py.sum()),accepted_draws=len(draws))
                    if window == "all":
                        b = make_bridge(scorer,panel,prob,result,checks,model=="prevalence_baseline")
                        for i, metric in enumerate(scorer.METRICS):
                            row = dict(meta,metric=metric)
                            for name,v in b.items():
                                row[name],row[name+"_lo"],row[name+"_hi"] = stats(v[i])
                            bridge_rows.append(row)
                            bridge_draws.append(pd.DataFrame(dict(scope=scope,direction=direction,model=model,
                                metric=metric,replicate=draws,**{name:v[i,1:] for name,v in b.items()})))
                        if not args.key_cases:
                            for cell,v in result.items():
                                for i,metric in enumerate(scorer.METRICS):
                                    key=dict(direction=direction,model=model,metric=metric,cell=cell)
                                    pt,lo,hi=stats(v[i])
                                    if scope=="full":
                                        standard_points.append(dict(key,value=pt))
                                        standard_draws.append(pd.DataFrame(dict(key,replicate=draws,value=v[i,1:])))
                                    else:
                                        esp_points.append(dict(key,value=pt,bootstrap_mean=np.mean(v[i,1:]),ci_low=lo,ci_high=hi))
                    if window == "first":
                        test_error(checks,"first_specimen_A_C",result["A"]-result["C"])
                        test_error(checks,"first_specimen_B_D",result["B"]-result["D"])
                    all_scores[model,window]=result
                    if target=="PRISM":
                        continue
                    for i,metric in enumerate(scorer.METRICS):
                        for cell,v in result.items():
                            pt,lo,hi=stats(v[i])
                            windows.append(dict(meta,window=window,metric=metric,object=cell,point=pt,ci_low=lo,ci_high=hi))
            if len(panels)>1:
                for model in frame.model.unique():
                    for left,right in [("first_two","first"),("all","first"),("all","first_two")]:
                        for cell in scorer.OBJECTS:
                            diff=all_scores[model,left][cell]-all_scores[model,right][cell]
                            for i,metric in enumerate(scorer.METRICS):
                                pt,lo,hi=stats(diff[i])
                                contrasts.append(dict(scope=scope,direction=direction,model=model,metric=metric,object=cell,
                                    left_window=left,right_window=right,point=pt,ci_low=lo,ci_high=hi,n_participants=len(py),accepted_draws=len(draws)))
            print(f"Scored {direction} {scope}",flush=True)
    if not args.key_cases:
        keys=["direction","model","metric","cell"]
        ds=read(base/"orthogonal_audit/orthogonal_cell_summary.csv").merge(pd.DataFrame(standard_points),on=keys,validate="1:1")
        if len(ds)!=504: raise ValueError("Incomplete full-target check")
        test_error(checks,"all_504_registered_points",ds.value-ds.plugin_point_estimate)
        dr=read(base/"orthogonal_audit/orthogonal_cell_draws.csv").merge(pd.concat(standard_draws),on=keys+["replicate"],validate="1:1")
        if len(dr)!=252000: raise ValueError("Incomplete registered draws")
        test_error(checks,"all_252000_registered_draws",dr.value_x-dr.value_y)
        esp=read(base/"nielsen_esp_audit/esp_only_four_object_summary.csv").rename(columns={"object":"cell"}).merge(pd.DataFrame(esp_points),on=keys,validate="1:1")
        if len(esp)!=168: raise ValueError("Incomplete Spain comparison")
        # Retain the actual precision/centre of archived ESP records for audit.
        esp.to_csv(args.output/"esp_archive_comparison.csv",index=False)
        for field in ["ci_low","ci_high"]:
            checks["esp_archive_max_"+field+"_difference"]=float((esp[field+"_x"]-esp[field+"_y"]).abs().max())
            test_error(checks,"ESP_"+field,esp[field+"_x"]-esp[field+"_y"])
        checks["esp_archive_max_plugin_difference"]=float((esp.value_x-esp.value_y).abs().max())
        checks["esp_archive_max_bootmean_difference"]=float((esp.value_x-esp.bootstrap_mean).abs().max())
        test_error(checks,"ESP_bootstrap_mean_all_168",esp.value_x-esp.bootstrap_mean)
    results={"bridge_summary":pd.DataFrame(bridge_rows),"bridge_draws":pd.concat(bridge_draws),
             "sampling_structure":pd.DataFrame(structure)}
    if windows:
        results.update(window_summary=pd.DataFrame(windows),window_paired=pd.DataFrame(contrasts))
    for name,t in results.items():
        filename=name+(".csv.gz" if name=="bridge_draws" else ".csv")
        compression={"method":"gzip","mtime":0} if filename.endswith(".gz") else None
        t.to_csv(args.output/filename,index=False,lineterminator="\n",compression=compression)
    metadata=dict(rows={n:len(t) for n,t in results.items()},checks=checks,
                  inputs={k:checksum(p) for k,p in files.items()},
                  policy="Linear percentile limits of paired target-participant draws. Integer participant multiplicities are shared across windows and objects within each target scope.")
    (args.output/"sampling_checks.json").write_text(json.dumps(metadata,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle",required=True,type=Path)
    parser.add_argument("--order",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--weights",type=Path,help="Alternate paired participant multiplicity file for key-case stability")
    parser.add_argument("--key-cases",action="store_true",help="Full-target decomposition and three specimen rules for HMP2-to-Nielsen L2 and Franzosa-to-HMP2 shallow MLP")
    run(parser.parse_args())
