import glob

import numpy as np
import tskit

print("div")
for treefile in glob.glob("moran/*.trees"):
    ts = tskit.load(treefile)
    # for t in ts.trees(): print(t.num_roots)
    params = eval(ts.provenance(0).record)
    samples = [i for i in ts.samples()]
    div = ts.diversity([samples], mode="branch")
    print(div[0])
