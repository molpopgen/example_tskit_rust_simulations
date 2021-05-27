import glob

import numpy as np
import tskit

print("N p recrate div")
for treefile in glob.glob("output/*.trees"):
    ts = tskit.load(treefile)
    params = eval(ts.provenance(0).record)
    div = ts.diversity([[i for i in ts.samples()]], mode="branch")
    print(params["N"], params["psurvival"], params["recrate"], div[0])
