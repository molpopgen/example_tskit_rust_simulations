import glob

import numpy as np
import tskit

print("p div")
for treefile in glob.glob("*.trees"):
    ts = tskit.load(treefile)
    params = eval(ts.provenance(0).record)
    samples = np.where(ts.tables.nodes.flags == 1)
    div = ts.diversity([samples[0]], mode="branch")
    mroots = 0
    for t in ts.trees():
        mroots = max(mroots, t.num_roots)
    print(params["psurvival"], div[0], mroots)
