import glob

import numpy as np
import tskit

print("p div")
for treefile in glob.glob("output/*.trees"):
    ts = tskit.load(treefile)
    params = eval(ts.provenance(0).record)
    div = ts.diversity([[i for i in ts.samples()]], mode="branch")
    print(params["psurvival"], div[0])
