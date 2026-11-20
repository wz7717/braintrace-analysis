"""Bo count preparation and donor-local feature eligibility for the R9 source analysis."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd


@dataclass
class SourceData:
    x: np.ndarray
    symbols: np.ndarray
    stable_row_symbols: np.ndarray
    detected_by_donor: dict
    samples: np.ndarray
    donors: np.ndarray
    networks: np.ndarray
    regions: np.ndarray

    def split(self, held_out):
        train = np.flatnonzero(~np.isin(self.donors, held_out))
        query = np.flatnonzero(np.isin(self.donors, held_out))
        detected = sum(self.detected_by_donor[d] for d in np.unique(self.donors[train]))
        eligible = np.unique(self.stable_row_symbols[10 * detected >= len(train)])
        return train, query, np.searchsorted(self.symbols, eligible)


def read_source(repo):
    from core.reference_projection import read_gene_map
    from scripts.run_bo2023_projected_vsd_formal_lomo import read_metadata

    data = Path(repo) / "external_data/Bo2023"
    counts_path = data / "mfas5_819samples_28415genes_featurecounts_counts.txt"
    with counts_path.open(encoding="utf8") as stream:
        samples = np.array(stream.readline().rstrip().split("\t"))
    counts = pd.read_csv(counts_path, sep="\t", index_col=0, header=None, skiprows=1)
    counts.columns = samples
    metadata = read_metadata(
        data / "Information of sequenced samples_update_full878_filter819.xlsx",
        "mfas5_819samples_phenSet4", "Region", "SaleemNetworks", "MonkeyID",
    ).set_index("sample_id").loc[samples]
    donors = metadata.monkey_id.to_numpy(str)
    gene_map = read_gene_map(Path(__file__).with_name("mfas5_ensembl97_gene_map.csv"))
    symbols = gene_map.set_index("gene_id").loc[counts.index, "gene_symbol"]
    sums = counts.groupby(symbols, sort=True).sum()
    raw = counts.to_numpy()
    library = raw.sum(axis=0)
    observed = {d: (raw[:, donors == d] > 0).sum(axis=1) for d in np.unique(donors)}
    expression = np.log1p(sums.to_numpy(dtype=np.float64) / library * 1e6).astype(np.float32)
    return SourceData(expression, sums.index.to_numpy(str), symbols.to_numpy(str), observed,
                      samples, donors, metadata.endpoint_label.to_numpy(str), metadata.region_id.to_numpy(str))


class Fisher:
    """Cache within-class residual sums and means once per training fold."""
    def __init__(self, expression, labels, train):
        self.labels = sorted(set(labels[train]))
        self.members = {label: train[labels[train] == label] for label in self.labels}
        self.means = {}
        self.residuals = {}
        for label, indices in self.members.items():
            values = expression[:, indices].astype(np.float64)
            mean = values.mean(axis=1)
            self.means[label] = mean
            self.residuals[label] = ((values - mean[:, None]) ** 2).sum(axis=1)

    def panel(self, universe, labels, size):
        counts = [len(self.members[label]) for label in labels]
        means = np.column_stack([self.means[label] for label in labels])
        center = np.average(means, axis=1, weights=counts)
        between = np.average((means - center[:, None]) ** 2, axis=1, weights=counts)
        within = sum(self.residuals[label] for label in labels) / max(sum(counts) - len(counts), 1)
        scores = between / (within + 1e-8)
        return universe[np.argsort(scores[universe])[::-1][:size]]

