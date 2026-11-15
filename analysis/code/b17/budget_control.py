"""Owned-run cumulative byte ledger and frozen wall/RSS limits."""
import json
from pathlib import Path

PRODUCTION_CAPS = {'whole_seconds':14400, 'build_seconds':3600, 'query_seconds':300,
                   'rss_bytes':8589934592, 'cache_bytes':2147483648, 'output_bytes':1073741824}
CACHE_NAMES = {'source_profiles.f64le','genes.txt','reference.tsv','reference_meta.tsv','SOURCE_PREPARED.json',
               'QUERY_AXIS.json','SHARED_PREPARED_READY.json','PEARSON_REFERENCE.json','SINGLER_REFERENCE.json',
               'trained.rds','marker_genes.tsv'}


class ByteLedger:
    def __init__(self, root, preparation_root=None):
        self.root = Path(root)
        self.preparation_root = Path(preparation_root) if preparation_root else None
        self.current = {}
        # Charge the visible preparation files below, plus a conservative 1 MiB
        # for pre-clock source rewrites during this small wiring revision.
        self.preparation_revision_allowance = 1048576 if preparation_root else 0
        self.written = {'cache':0, 'output':self.preparation_revision_allowance}
        self.files = {}

    def role(self, path):
        name = path.name
        return 'cache' if name in CACHE_NAMES or name.endswith('.f64le') or name.endswith('_query.tsv') else 'output'

    def note_write(self, path, size, role=None):
        path = Path(path)
        role = role or self.role(path)
        if path in self.files and self.files[path] != role:
            raise RuntimeError('A physical path cannot change budget category')
        self.files[path] = role
        self.written[role] += size
        self.current[path] = size

    def sample(self):
        paths = list(self.root.rglob('*'))
        if self.preparation_root:
            paths.extend(self.preparation_root.rglob('*'))
        for path in paths:
            if not path.is_file():
                continue
            try:
                size = path.stat().st_size
            except FileNotFoundError:
                continue
            role = self.files.get(path, self.role(path) if path.is_relative_to(self.root) else 'output')
            self.files[path] = role
            self.written[role] += max(0, size - self.current.get(path, 0))
            self.current[path] = size
        # Missing/deleted files retain all charged writes. Copies have distinct paths.
        return dict(self.written)

    def write_json(self, path, value):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        data = json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False).encode('utf-8')
        path.write_bytes(data)
        self.note_write(path, len(data), 'output')

    def snapshot(self):
        return {'cumulative_bytes': self.sample(), 'physical_paths_seen':len(self.files),
                'deleted_files_refund_bytes':0,
                'preparation_revision_allowance_bytes':self.preparation_revision_allowance,
                'scope':'Owned run and preparation; sampled child growth, exact controller/delivery writes, conservative pre-clock revision allowance',
                'file_roles':{str(path.relative_to(self.root)) if path.is_relative_to(self.root) else str(path):role for path,role in self.files.items()}}
