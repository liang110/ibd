"""Verify the executable research-use declarations and their input checks on the Nielsen example."""
import argparse
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pandas as pd

from score_declaration import BUNDLE, HERE, evaluate, read_input, validate_declaration


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--unit",type=Path,required=True)
    opts=parser.parse_args()
    example=json.loads((HERE/"nielsen_example.json").read_text(encoding="utf-8"))
    src=BUNDLE/"evidence/full_reproducibility_inputs"
    rows=read_input(src/"prediction_registry.csv")
    rows=rows.loc[rows.direction.eq(example["direction"]) & rows.model.eq(example["model"])]
    weights=read_input(src/"shared_bootstrap_weights.csv")
    order=read_input(BUNDLE/"analysis/unit/input/retained_sampling_order.csv")
    windows=pd.read_csv(opts.unit/"window_summary.csv",float_precision="round_trip")
    windows=windows.loc[windows.direction.eq(example["direction"]) & windows.model.eq(example["model"]) & windows.scope.eq("full")]
    bridge=pd.read_csv(opts.unit/"bridge_summary.csv",float_precision="round_trip")
    bridge=bridge.loc[bridge.direction.eq(example["direction"]) & bridge.model.eq(example["model"]) & bridge.scope.eq("full")].set_index("metric")
    checks=[]
    for rule,when,window,n in [("all_retained","after_collection_of_all_retained_specimens","all",396),
                               ("earliest_retained","at_first_retained_collection","first",318),
                               ("first_two_retained","after_up_to_two_retained_collections","first_two",396)]:
        config=dict(example,eligibility=rule,assessment_time=when)
        scores, contrasts=evaluate(config,rows,weights,order)
        reference=windows.loc[windows.window.eq(window)]
        measured=scores.loc[scores.object.ne("E")].merge(reference,on=["object","metric"],validate="one_to_one")
        errs=[np.abs(measured[field+"_x"]-measured[field+"_y"]).max() for field in ("point","ci_low","ci_high")]
        e = scores.loc[scores.object.eq("E")].set_index("metric")
        if window=="first":
            ref=scores.loc[scores.object.eq("C")].set_index("metric")
            errs.extend(np.abs(e[col]-ref[col]).max() for col in ("point","ci_low","ci_high"))
        else:
            # Each Nielsen participant has at most two retained specimens.
            errs.extend(np.abs(e[col]-bridge[field]).max() for col,field in [("point","E"),("ci_low","E_lo"),("ci_high","E_hi")])
        if len(measured)!=12 or len(e)!=3 or len(contrasts)!=21 or not scores.n_specimens.eq(n).all() or max(errs)>1e-12:
            raise ValueError(f"Declared choices differ from scientific output for {window}")
        checks.append(dict(check=rule,outcome="agrees_with_reference",records=len(measured)+len(e),maximum_absolute_error=max(errs)))
    def reject(name,action):
        try:
            action()
        except ValueError:
            checks.append(dict(check=name,outcome="rejected",records=1,maximum_absolute_error=None))
        else:
            raise ValueError(f"Invalid input accepted: {name}")
    bad=deepcopy(example)
    bad["rules"][0]["target_weighting"]="equal_participants"
    reject("incompatible A and E weighting",lambda:validate_declaration(bad))
    reject("undefined scoring declaration version",lambda:validate_declaration(dict(example,schema_version="0")))
    reject("first assessment using subsequent specimens",lambda:validate_declaration(dict(example,assessment_time="at_first_retained_collection")))
    reject("unstable endpoint requested",lambda:validate_declaration(dict(example,label_policy="time_varying")))
    reject("missing expected run",lambda:evaluate(example,rows.loc[rows.split.ne(1)],weights,None))
    reject("duplicate run and specimen",lambda:evaluate(example,pd.concat([rows,rows.iloc[:1]]),weights,None))
    label_conflict=rows.copy()
    label_conflict.iloc[0,label_conflict.columns.get_loc("outcome_ibd")]=1-int(label_conflict.outcome_ibd.iloc[0])
    reject("inconsistent labels",lambda:evaluate(example,label_conflict,weights,None))
    first=dict(example,eligibility="earliest_retained",assessment_time="at_first_retained_collection")
    reject("missing order sidecar",lambda:evaluate(first,rows,weights,None))
    omit_sample=rows.sample_id.iloc[0]
    reject("absent retained order",lambda:evaluate(first,rows,weights,order.loc[order.sample_id.ne(omit_sample)]))
    tied=order.copy()
    tied.loc[tied.target.eq("NielsenHB_2014"),"retained_order"]=1
    reject("ambiguous retained order",lambda:evaluate(first,rows,weights,tied))
    reject("mismatched draw count",lambda:evaluate(dict(example,expected_draws=2000),rows,weights,None))
    reversed_labels=dict(example,positive_label=0)
    reject("undeclared reverse positive encoding",lambda:validate_declaration(reversed_labels))
    pd.DataFrame(checks).to_csv(opts.output,index=False)


if __name__=="__main__":
    main()
