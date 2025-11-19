import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(sys.argv[1], sep=r"\s+")
df = df.dropna(subset=["P"])
df = df[df["P"] > 0]
df = df.sort_values(["CHR", "BP"])
df["pos"] = range(len(df))

plt.figure(figsize=(14,6))
plt.scatter(df["pos"], -np.log10(df["P"]), c=df["CHR"] % 2, s=4, cmap="coolwarm")
plt.axhline(y=3, color='red', linestyle='--')
plt.xlabel("Genomic position")
plt.ylabel("-log10(P)")
plt.title("Manhattan plot (GLM)")
plt.tight_layout()
plt.savefig(sys.argv[2])