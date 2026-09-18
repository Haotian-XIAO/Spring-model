from __future__ import annotations
import numpy as np

def gini(x: np.ndarray) -> float:
    x=np.sort(np.asarray(x,float)); n=x.size
    if not n or x.sum()==0: return float("nan")
    return float((2*np.dot(np.arange(1,n+1),x)/(n*x.sum()))-(n+1)/n)

def cell_neighbors(hex_cells: np.ndarray) -> set[tuple[int,int]]:
    edge_owner={}; pairs=set()
    for ci,cell in enumerate(hex_cells):
        for j in range(len(cell)):
            edge=tuple(sorted((int(cell[j]),int(cell[(j+1)%len(cell)]))))
            if edge in edge_owner: pairs.add((edge_owner[edge],ci))
            else: edge_owner[edge]=ci
    return pairs

def metrics(snapshot: dict[str,np.ndarray], lower: float=.75, upper: float=1.25, reference_area: float | None = None) -> dict[str,float]:
    a=np.asarray(snapshot["cell_area"],float)
    n=a / float(reference_area) if reference_area is not None else np.asarray(snapshot["normalized_cell_area"],float)
    pairs=cell_neighbors(snapshot["hex_cells"])
    contrast=float(np.mean([abs(n[i]-n[j]) for i,j in pairs])) if pairs else float("nan")
    preserved=(n>=lower)&(n<=upper); adjacency=cell_neighbors(snapshot["hex_cells"])
    parent=np.arange(len(a))
    def root(i):
        while parent[i]!=i: parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in adjacency:
        if preserved[i] and preserved[j]:
            ri,rj=root(i),root(j)
            if ri!=rj: parent[rj]=ri
    sizes={}
    for i in np.flatnonzero(preserved): sizes[root(i)]=sizes.get(root(i),0)+1
    largest=max(sizes.values(),default=0)
    return {
      "cell_count":float(len(a)),"total_cell_area":float(a.sum()),"mean_cell_area":float(a.mean()),
      "median_cell_area":float(np.median(a)),"cell_area_sd":float(a.std()),
      "cell_area_cv":float(a.std()/a.mean()),"cell_area_iqr":float(np.percentile(a,75)-np.percentile(a,25)),
      "cell_area_gini":gini(a),"normalized_area_mean":float(n.mean()),"normalized_area_median":float(np.median(n)),
      "reference_area":float(reference_area) if reference_area is not None else float("nan"),
      "normalized_area_q05":float(np.percentile(n,5)),"normalized_area_q95":float(np.percentile(n,95)),
      "normalized_area_cv":float(n.std()/n.mean()),"neighbor_normalized_area_contrast":contrast,
      "operational_preserved_cell_fraction":float(np.mean((n>=lower)&(n<=upper))),
      "operational_preserved_area_fraction":float(a[(n>=lower)&(n<=upper)].sum()/a.sum()),
      "operational_degraded_low_fraction":float(np.mean(n<lower)),
      "operational_enlarged_fraction":float(np.mean(n>upper)),
      "operational_preserved_component_count":float(len(sizes)),
      "operational_largest_preserved_component_fraction":float(largest/max(int(preserved.sum()),1)),
    }
