#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt

inp = sys.argv[1]
out = sys.argv[2]

df = pd.read_csv(inp, delim_whitespace=True, header=None)
df.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, df.shape[1]-1)]

plt.figure(figsize=(6, 5))
plt.scatter(df["PC1"], df["PC2"], s=20)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA of soybean lines")
plt.tight_layout()
plt.savefig(out, dpi=300)
