import tskit
import numpy as np
import glob

print("p div")
for treefile in glob.glob('*.trees'):
    ts = tskit.load(treefile)
    s = treefile.split('_')
    psurvival = float(s[1])
    samples = np.where(ts.tables.nodes.flags == 1)
    div = ts.diversity([samples[0]], mode="branch")
    mroots = 0
    for t in ts.trees():
        mroots = max(mroots, t.num_roots)
    print(psurvival, div[0], mroots)

