"""Score a declared empirical target using probabilities, participant draws and optional specimen order."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from jsonschema import Draft202012Validator

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parents[1]
SCHEMA = json.loads((HERE/"scoring_declaration.schema.json").read_text(encoding="utf-8"))
VALIDATOR = Draft202012Validator(SCHEMA)
PAIRS = [("B","A"),("D","C"),("C","A"),("D","B"),("D","A"),("E","A"),("C","E")]


def load_scoring():
    path = BUNDLE/"analysis/unit/score_sampling.py"
    spec = importlib.util.spec_from_file_location("sampling",path)
    sampling = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sampling)
    return sampling, sampling.load_scorer(BUNDLE)


def validate_declaration(declaration):
    error = next(VALIDATOR.iter_errors(declaration), None)
    if error:
        raise ValueError(f"Scoring declaration at {error.json_path}: {error.message}")
    if declaration["direction"].split("_to_")[-1] != declaration["target"]:
        raise ValueError("Direction and target identifiers disagree")


def read_input(path):
    # Header-based string types preserve identifiers such as 001 in every input.
    cols = pd.read_csv(path, nrows=0).columns
    ids = {"direction","target","model","sample_id","bootstrap_cluster","participant_id"}
    table = pd.read_csv(path, float_precision="round_trip", keep_default_na=False,
                        dtype={field:"str" for field in cols if field in ids})
    aliases = {"run_id":"split","participant_id":"bootstrap_cluster","y_true":"outcome_ibd",
               "y_prob":"probability","replicate_id":"replicate","multiplicity":"weight"}
    if any(old in cols and new in cols for old,new in aliases.items()):
        raise ValueError("Both a canonical and alternate name appear for the same field")
    return table.rename(columns=aliases)


def select_records(frame, order, declaration):
    required = ["sample_id","bootstrap_cluster","split","outcome_ibd","probability"]
    frame = frame.loc[(frame.direction == declaration["direction"]) & (frame.model == declaration["model"]), required]
    if frame.empty or frame.duplicated(["sample_id","split"]).any() or frame.isna().any().any():
        raise ValueError("Empty, missing or duplicate prediction records")
    if frame.sample_id.eq("").any() or frame.bootstrap_cluster.eq("").any():
        raise ValueError("Missing specimen or participant identifiers")
    if set(frame.split) != set(declaration["expected_runs"]):
        raise ValueError("Incomplete source-run panel")
    if not frame.outcome_ibd.isin([0,1]).all() or not frame.probability.between(0,1).all():
        raise ValueError("Invalid positive-class probability or label")
    if not frame.groupby("sample_id")[["bootstrap_cluster","outcome_ibd"]].nunique().eq(1).all().all():
        raise ValueError("Conflicting specimen ownership or labels")
    if not frame.groupby("bootstrap_cluster").outcome_ibd.nunique().eq(1).all():
        raise ValueError("Conflicting participant labels")
    if not frame.groupby("sample_id").split.size().eq(len(declaration["expected_runs"])).all():
        raise ValueError("Incomplete run-by-specimen panel")
    eligibility = declaration["eligibility"]
    if eligibility == "all_retained":
        return frame
    if order is None:
        raise ValueError("Restricted specimen rules require the specimen-order table")
    ranks = order.loc[order.target.eq(declaration["target"]),["sample_id","retained_order"]]
    if ranks.sample_id.duplicated().any():
        raise ValueError("Duplicated specimen-order keys")
    samples = frame[["sample_id","bootstrap_cluster"]].drop_duplicates().merge(ranks,on="sample_id",how="left",validate="one_to_one")
    if not samples.retained_order.notna().all() or samples.duplicated(["bootstrap_cluster","retained_order"]).any():
        raise ValueError("Missing or tied within-participant order")
    for participant, group in samples.groupby("bootstrap_cluster"):
        if not np.array_equal(np.sort(group.retained_order), np.arange(1,len(group)+1)):
            raise ValueError(f"Incomplete eligibility order for participant {participant}")
    ranks_to_keep = 1 if eligibility == "earliest_retained" else 2
    return frame.loc[frame.sample_id.isin(samples.loc[samples.retained_order.le(ranks_to_keep),"sample_id"])]


def evaluate(declaration, records, weights, order):
    validate_declaration(declaration)
    frame = select_records(records, order, declaration)
    local = weights.loc[weights.direction.eq(declaration["direction"])]
    if local.empty or local.replicate.nunique() != declaration["expected_draws"]:
        raise ValueError("Requested draw set does not match the file")
    if local.duplicated(["replicate","bootstrap_cluster"]).any() or local.isna().any().any():
        raise ValueError("Missing or duplicate participant multiplicities")
    sampling, scorer = load_scoring()
    panel = sampling.prepare_panel(scorer,frame,local)
    samples, people, groups, sw, pw, draws, y, py, balanced = panel
    probabilities = frame.pivot(index="split",columns="sample_id",values="probability").reindex(columns=samples.sample_id).to_numpy()
    result = scorer.four_scores(probabilities,y,py,groups,sw,pw)
    objects = [row["object"] for row in declaration["rules"]]
    if "E" in objects:
        result["E"] = np.mean([scorer.score(y,run,balanced) for run in probabilities],axis=0)
    score_rows, contrast_rows = [], []
    common = dict(direction=declaration["direction"],model=declaration["model"],eligibility=declaration["eligibility"],
                  n_specimens=len(samples),n_participants=len(people),n_runs=len(probabilities),accepted_draws=len(draws))
    def append(rows, tag, values):
        for metric in declaration["metrics"]:
            v = values[scorer.METRICS.index(metric)]
            point, low, high = sampling.stats(v)
            rows.append(dict(**common,**tag,metric=metric,point=point,bootstrap_mean=float(np.mean(v[1:])),ci_low=low,ci_high=high))
    for rule in declaration["rules"]:
        append(score_rows,rule,result[rule["object"]])
    for left,right in PAIRS:
        if left in objects and right in objects:
            append(contrast_rows,dict(left_object=left,right_object=right),result[left]-result[right])
    empty_columns = [*common,"left_object","right_object","metric","point","bootstrap_mean","ci_low","ci_high"]
    return pd.DataFrame(score_rows), pd.DataFrame(contrast_rows,columns=empty_columns)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    default = BUNDLE/"evidence/full_reproducibility_inputs"
    parser.add_argument("--predictions",type=Path,default=default/"prediction_registry.csv")
    parser.add_argument("--weights",type=Path,default=default/"shared_bootstrap_weights.csv")
    parser.add_argument("--order",type=Path,default=BUNDLE/"analysis/unit/input/retained_sampling_order.csv")
    parser.add_argument("--declaration",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--contrasts-output",type=Path)
    args=parser.parse_args()
    declaration=json.loads(args.declaration.read_text(encoding="utf-8"))
    validate_declaration(declaration)
    order=read_input(args.order) if declaration["eligibility"]!="all_retained" else None
    scores, contrasts=evaluate(declaration,read_input(args.predictions),read_input(args.weights),order)
    destination=args.contrasts_output or args.output.with_name(args.output.stem+"_contrasts.csv")
    if destination.resolve()==args.output.resolve():
        raise ValueError("Choose distinct paths for scores and contrasts")
    for table,path in [(scores,args.output),(contrasts,destination)]:
        path.parent.mkdir(parents=True,exist_ok=True)
        table.to_csv(path,index=False)


if __name__=="__main__":
    main()
