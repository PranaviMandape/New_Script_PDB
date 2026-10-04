#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %% [markdown]
# # PDB crystal structure + AlphaFold collection pipeline (SBDD)
#
# Give it **any protein name** (or gene symbol / synonym / UniProt accession). It will
#
# 1. resolve the protein identity (UniProt) and collect aliases,
# 2. search **RCSB PDB** (identity-validated, *Homo sapiens*, *X-RAY DIFFRACTION* only),
# 3. **always and independently** search the **AlphaFold DB**,
# 4. write one Excel workbook with `PDB_Data` and `AlphaFold_Data` (plus audit sheets).
#
# No LLM / paid API is used. Everything is deterministic and comes from official sources
# (UniProt REST, RCSB Search + Data API, RCSB/wwPDB validation files, AlphaFold DB API).
#
# **Run (script):** `python protein_structure_pipeline.py "Cyclin-dependent kinase 4"`
# **Install:** `pip install requests openpyxl tqdm pypdf`

# %%
"""
protein_structure_pipeline.py
=============================

Architecture (each step is its own function / group of functions)
-----------------------------------------------------------------
 1. Input + normalisation ............ normalize_protein_name()
 2. Identity resolution (UniProt) .... resolve_protein_identity(), find_protein_synonyms()
 3. RCSB search (paginated) .......... search_rcsb()
 4. PDB identity validation .......... validate_target_identity()
 5. PDB filtering .................... filter_pdb_structures() / process_pdb_entry()
 6. PDB metadata ..................... retrieve_pdb_metadata(), retrieve_residue_counts(),
                                       retrieve_depositor_r_values()
 7. Validation data .................. retrieve_validation_data(), retrieve_average_b()
 8. Ligands / activity ............... retrieve_ligands(), retrieve_ligand_activity()
 9. AlphaFold search ................. search_alphafold()
10. AlphaFold validation ............. retrieve_alphafold_data()
11. Data validation ................. validate_records()
12. Excel generation ................ create_excel_workbook()
"""

# %%
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import logging
import re
import sys
import threading
import time
import unicodedata
import xml.etree.ElementTree as ET
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import requests
from requests.adapters import HTTPAdapter
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

try:  # progress bars are optional: the pipeline works without tqdm
    from tqdm.auto import tqdm
except Exception:  # pragma: no cover
    class tqdm:  # minimal stand-in
        def __init__(self, iterable=None, total=None, desc="", **_kw):
            self.iterable, self.total, self.desc, self.n = iterable, total, desc, 0

        def __iter__(self):
            for item in self.iterable:
                yield item
                self.update(1)

        def update(self, k=1):
            self.n += k
            if self.total and (self.n == self.total or self.n % 25 == 0):
                print(f"  {self.desc}: {self.n}/{self.total}")

        def close(self):
            pass

        @staticmethod
        def write(msg, file=None):
            print(msg, file=file or sys.stdout)

# %% [markdown]
# ## 1. Configuration
# Everything you may want to tune lives here.

# %%
# ----------------------------- user-tunable configuration -----------------------------
TARGET_ORGANISM = "Homo sapiens"
TARGET_TAXONOMY_ID = 9606
EXPERIMENTAL_METHOD = "X-RAY DIFFRACTION"
MISSING_VALUE = "N/A"

MAX_RETRIES = 3                 # attempts per HTTP request (never retried forever)
REQUEST_TIMEOUT = 30            # seconds
BACKOFF_BASE_SECONDS = 1.5      # exponential back-off: 1.5s, 3s, 6s ...
MAX_WORKERS = 16                # concurrent HTTP workers
SEARCH_PAGE_SIZE = 1000         # rows per RCSB Search API page

CACHE_DIR = Path(".structure_cache")
CACHE_TTL_DATA_SECONDS = 7 * 24 * 3600     # entry/entity/chem-comp/validation data
CACHE_TTL_SEARCH_SECONDS = 24 * 3600       # search results change with weekly releases
CACHE_TTL_UNIPROT_SECONDS = 24 * 3600

MAX_ALIAS_QUERIES = 15          # max name-phrase queries sent to RCSB
MIN_ALIAS_LENGTH = 3
MIN_LIGAND_HEAVY_ATOMS = 6      # components with <= (this - 1) heavy atoms are treated as ions/solvents
PEPTIDE_MAX_LENGTH = 30         # polymer partner <= this many residues => "peptide"
PDF_MAX_PAGES_SCANNED = 40      # validation PDF pages scanned for "Average B, all atoms"
IDENTITY_MIN_CONFIDENT_SCORE = 60

# wwPDB R-free shown in the RCSB validation slider / report overview. The wwPDB validation
# pipeline recalculates R-free (DCC); the depositor value is kept in a separate column.
WWPDB_RFREE_FIELDS = ["DCC_Rfree"]

# Binding-affinity provenance to report. () = every source RCSB integrates
# (BindingDB, PDBbind, Binding MOAD); e.g. ("BindingDB",) restricts to BindingDB.
AFFINITY_SOURCES: Tuple[str, ...] = ()

# AlphaFold DB confidence bands (official pLDDT legend)
PLDDT_BANDS = [(90, "Very high"), (70, "Confident"), (50, "Low"), (0, "Very low")]

# ----------------------------- endpoints -----------------------------
USER_AGENT = "SBDD-Structure-Collector/1.0 (+python-requests)"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core"
RCSB_HEADER_CIF_URL = "https://files.rcsb.org/header/{pdb_id}.cif"
RCSB_STRUCTURE_URL = "https://www.rcsb.org/structure/{pdb_id}"
VALIDATION_PDF_URLS = (
    "https://files.rcsb.org/pub/pdb/validation_reports/{mid}/{pid}/{pid}_full_validation.pdf",
    "https://files.rcsb.org/pub/pdb/validation_reports/{mid}/{pid}/{pid}_full_validation.pdf.gz",
    "https://files.wwpdb.org/pub/pdb/validation_reports/{mid}/{pid}/{pid}_full_validation.pdf.gz",
)
VALIDATION_XML_URLS = (
    "https://files.rcsb.org/pub/pdb/validation_reports/{mid}/{pid}/{pid}_validation.xml.gz",
    "https://files.wwpdb.org/pub/pdb/validation_reports/{mid}/{pid}/{pid}_validation.xml.gz",
)
ALPHAFOLD_PREDICTION_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

# ----------------------------- ligand exclusion lists -----------------------------
# Components that are NOT treated as ligands of interest (water, ions, salts, buffers,
# solvents, cryo-/crystallisation agents, common additives, detergents, simple sugars).
_WATER = {"HOH", "DOD", "WAT"}
_IONS = {
    "NA", "K", "CL", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO", "CD", "HG",
    "BR", "IOD", "F", "LI", "RB", "CS", "SR", "BA", "AL", "PB", "PT", "AU", "AG", "YB", "GD",
    "TB", "SM", "LA", "XE", "KR", "AR", "NH4", "UNX", "UNL",
}
_SALTS_BUFFERS = {
    "SO4", "PO4", "NO3", "SCN", "BCT", "CO3", "ACT", "ACY", "FMT", "CIT", "FLC", "TLA", "TAR",
    "MLI", "MLA", "OXL", "TRS", "EPE", "MES", "MOP", "BTB", "CAC", "TAM", "IMD", "TFA", "PIN",
    "HEZ", "CXS", "AZI", "PEP",
}
_SOLVENTS_CRYSTALLISATION = {
    "EDO", "GOL", "PEG", "PGE", "PG4", "PG5", "PG6", "1PE", "2PE", "P6G", "PE3", "PE4", "PE5",
    "PE8", "PEU", "PGO", "15P", "7PE", "XPE", "MPD", "MRD", "BU1", "BU3", "PDO", "DMS", "DMF",
    "EOH", "IPA", "MOH", "ACN", "BME", "DTT", "DTV", "TCE", "DIO", "1BO", "SUC", "TRE",
}
_DETERGENTS_LIPIDS = {"OLC", "LMT", "BOG", "DMU", "C8E", "LDA", "OLA", "PLM", "MYR", "STE", "BNG", "UNL"}
_GLYCANS = {"NAG", "NDG", "MAN", "BMA", "FUC", "GAL", "GLC", "XYP", "SIA", "FUL", "A2G"}
EXCLUDED_COMPONENT_IDS = (
    _WATER | _IONS | _SALTS_BUFFERS | _SOLVENTS_CRYSTALLISATION | _DETERGENTS_LIPIDS | _GLYCANS
)

# %% [markdown]
# ## 2. Logging, caching and HTTP layer
# * one `requests.Session` per thread (connection pooling),
# * timeout + retry with exponential back-off + `Retry-After` support,
# * disk cache with a TTL (stale entries are refetched; `--refresh` bypasses the cache),
# * **no data** (HTTP 204/404) is clearly separated from **failure** (`APIError`).

# %%
log = logging.getLogger("structure_pipeline")


class _TqdmHandler(logging.Handler):
    def emit(self, record):
        try:
            tqdm.write(self.format(record), file=sys.stdout)
        except Exception:
            self.handleError(record)


def setup_logging(verbose: bool = False) -> None:
    log.handlers.clear()
    handler = _TqdmHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    log.addHandler(handler)
    log.setLevel(logging.DEBUG if verbose else logging.WARNING)
    log.propagate = False
    logging.getLogger("pypdf").setLevel(logging.ERROR)


class APIError(Exception):
    """Raised when a request could not be completed (network error, 5xx, 4xx, ...)."""


@dataclass
class HttpResult:
    status: int
    content: bytes

    @property
    def text(self) -> str:
        return self.content.decode("utf-8", errors="replace")

    def json(self) -> Any:
        # Decimal keeps the exact digits the API emitted (e.g. 0.231, 28.0)
        return json.loads(self.text, parse_float=Decimal)


class DiskCache:
    """Tiny file cache with TTL. Stores (status, bytes) per key."""

    def __init__(self, directory: Path, enabled: bool = True, refresh: bool = False):
        self.dir = Path(directory)
        self.enabled = enabled
        self.refresh = refresh
        if enabled:
            self.dir.mkdir(parents=True, exist_ok=True)

    def _path(self, key: str) -> Path:
        h = hashlib.sha1(key.encode("utf-8")).hexdigest()
        return self.dir / h[:2] / h

    def get(self, key: str, ttl: Optional[float]) -> Optional[Tuple[int, bytes]]:
        if not self.enabled or self.refresh or not ttl:
            return None
        p = self._path(key)
        try:
            if not p.exists() or (time.time() - p.stat().st_mtime) > ttl:
                return None
            raw = p.read_bytes()
            head, _, body = raw.partition(b"\n")
            return int(head), body
        except Exception:
            return None

    def put(self, key: str, status: int, content: bytes) -> None:
        if not self.enabled:
            return
        p = self._path(key)
        try:
            p.parent.mkdir(parents=True, exist_ok=True)
            tmp = p.with_suffix(f".{threading.get_ident()}.tmp")
            tmp.write_bytes(str(status).encode() + b"\n" + content)
            tmp.replace(p)
        except Exception as exc:  # cache problems must never break the run
            log.debug("cache write failed: %s", exc)

    # derived (already-parsed) values, e.g. values extracted from a CIF header or PDF
    def get_json(self, key: str, ttl: Optional[float]) -> Optional[dict]:
        hit = self.get("derived|" + key, ttl)
        if hit is None:
            return None
        try:
            return json.loads(hit[1].decode("utf-8"))
        except Exception:
            return None

    def put_json(self, key: str, value: dict) -> None:
        self.put("derived|" + key, 200, json.dumps(value).encode("utf-8"))


class HttpClient:
    def __init__(self, cache: DiskCache, workers: int = MAX_WORKERS):
        self.cache = cache
        self.workers = workers
        self._local = threading.local()

    def _session(self) -> requests.Session:
        s = getattr(self._local, "session", None)
        if s is None:
            s = requests.Session()
            adapter = HTTPAdapter(pool_connections=8, pool_maxsize=max(8, self.workers * 2), max_retries=0)
            s.mount("https://", adapter)
            s.mount("http://", adapter)
            s.headers.update({"User-Agent": USER_AGENT})
            self._local.session = s
        return s

    @staticmethod
    def _cache_key(method, url, params, json_body) -> str:
        return "|".join([
            method, url,
            json.dumps(params, sort_keys=True) if params else "",
            json.dumps(json_body, sort_keys=True) if json_body else "",
        ])

    def request(self, method: str, url: str, *, params=None, json_body=None, headers=None,
                ttl: Optional[float] = None, timeout: float = REQUEST_TIMEOUT) -> HttpResult:
        """Return HttpResult for 200/204/404. Raise APIError for anything else after retries."""
        key = self._cache_key(method, url, params, json_body)
        hit = self.cache.get(key, ttl)
        if hit is not None:
            return HttpResult(*hit)

        last_error = "unknown error"
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                resp = self._session().request(
                    method, url, params=params, json=json_body, headers=headers, timeout=timeout)
            except requests.RequestException as exc:  # timeouts, connection resets, DNS ...
                last_error = f"{type(exc).__name__}: {exc}"
                self._sleep(attempt, None)
                continue

            if resp.status_code in (200, 204, 404):
                if ttl:
                    self.cache.put(key, resp.status_code, resp.content)
                return HttpResult(resp.status_code, resp.content)
            if resp.status_code in (429, 500, 502, 503, 504):
                last_error = f"HTTP {resp.status_code}"
                self._sleep(attempt, resp.headers.get("Retry-After"))
                continue
            raise APIError(f"HTTP {resp.status_code} for {method} {url}")
        raise APIError(f"{method} {url} failed after {MAX_RETRIES} attempts ({last_error})")

    @staticmethod
    def _sleep(attempt: int, retry_after: Optional[str]) -> None:
        if attempt >= MAX_RETRIES:
            return
        delay = BACKOFF_BASE_SECONDS * (2 ** (attempt - 1))
        if retry_after:
            try:
                delay = max(delay, min(float(retry_after), 60.0))
            except ValueError:
                pass
        time.sleep(delay)

    def get(self, url, **kw) -> HttpResult:
        return self.request("GET", url, **kw)

    def post_json(self, url, body, **kw) -> HttpResult:
        return self.request("POST", url, json_body=body, **kw)

# %% [markdown]
# ## 3. Small utilities
# Name normalisation, case-insensitive dict access, Decimal handling (so source precision survives).

# %%
def normalize_protein_name(name: str) -> str:
    """Trim, unicode-normalise, drop wrapping quotes, collapse whitespace."""
    s = unicodedata.normalize("NFKC", name or "")
    s = s.strip().strip("\"'`")
    return re.sub(r"\s+", " ", s).strip()


def norm_key(s: Any) -> str:
    """Aggressive comparison key: 'Cyclin-dependent kinase 4' == 'cyclin dependent kinase 4' == 'CDK 4'-style."""
    return re.sub(r"[^0-9a-z]+", "", unicodedata.normalize("NFKC", str(s)).casefold())


def ci_get(d: Any, *names: str, default=None):
    """Case/underscore/hyphen-insensitive dictionary lookup (robust to schema naming drift)."""
    if not isinstance(d, dict):
        return default
    wanted = {norm_key(n) for n in names}
    for k, v in d.items():
        if norm_key(k) in wanted and v is not None:
            return v
    return default


def to_decimal(x: Any) -> Optional[Decimal]:
    if x is None or isinstance(x, bool):
        return None
    if isinstance(x, Decimal):
        return x if x.is_finite() else None
    if isinstance(x, int):
        return Decimal(x)
    if isinstance(x, float):
        return Decimal(repr(x))
    s = str(x).strip()
    if s in ("", "?", ".", MISSING_VALUE):
        return None
    try:
        d = Decimal(s)
    except InvalidOperation:
        return None
    return d if d.is_finite() else None


def decimals_of(d: Decimal) -> int:
    exp = d.as_tuple().exponent
    return max(0, -exp) if isinstance(exp, int) else 0


def fmt_decimal(d: Decimal) -> str:
    return format(d, "f")


def excel_number_format(d: Decimal, cap: int = 8) -> str:
    n = min(decimals_of(d), cap)
    return "0" if n == 0 else "0." + "0" * n


def choose_source_value(api_val: Optional[Decimal], text_val: Optional[Decimal]) -> Optional[Decimal]:
    """
    The Data API gives the official numeric value but may drop trailing zeros ("0.20" -> 0.2).
    The mmCIF header keeps the original text. Use the text version ONLY if it is numerically
    identical to the API value (or the API has none); otherwise keep the API value.
    """
    if text_val is not None and (api_val is None or text_val == api_val):
        return text_val
    return api_val


PDB_ID_RE = re.compile(r"^(?:[0-9][A-Za-z0-9]{3}|pdb_[0-9]{4}[A-Za-z0-9]{4})$")
UNIPROT_ACC_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def parse_json(res: HttpResult) -> Any:
    return res.json()


def sanitize_filename(name: str, max_len: int = 150) -> str:
    s = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", name)
    s = re.sub(r"\s+", "_", s.strip()).strip("._ ")
    if not s:
        s = "protein"
    if s.split(".")[0].upper() in {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
                                   *(f"LPT{i}" for i in range(1, 10))}:
        s = "_" + s
    return s[:max_len]

# %% [markdown]
# ## 4. Protein identity resolution (UniProt)
# Input → normalise → UniProt (human) → official name, gene symbol, aliases, accession(s).
# Candidates are scored by **exact normalised match** against recommended/alternative names and gene
# names/synonyms; reviewed (Swiss-Prot) entries get a small bonus. Nothing is guessed silently:
# low-confidence matches are flagged and listed.

# %%
@dataclass
class ProteinIdentity:
    input_name: str
    normalized: str
    status: str = "UNRESOLVED"            # RESOLVED | UNRESOLVED | API_ERROR
    message: str = ""
    official_name: Optional[str] = None
    gene_symbol: Optional[str] = None
    gene_synonyms: List[str] = field(default_factory=list)
    aliases: List[str] = field(default_factory=list)
    accession: Optional[str] = None
    secondary_accessions: List[str] = field(default_factory=list)
    uniprot_entry_name: Optional[str] = None
    reviewed_label: Optional[str] = None
    organism: Optional[str] = None
    taxon_id: Optional[int] = None
    match_score: int = 0
    match_basis: str = ""
    alternatives: List[str] = field(default_factory=list)

    @property
    def all_accessions(self) -> List[str]:
        accs = [self.accession] if self.accession else []
        return list(dict.fromkeys(accs + self.secondary_accessions))

    @property
    def alias_keys(self) -> set:
        keys = {norm_key(a) for a in self.aliases + [self.official_name or "", self.normalized]}
        keys.discard("")
        return keys

    @property
    def gene_keys(self) -> set:
        keys = {norm_key(g) for g in ([self.gene_symbol] if self.gene_symbol else []) + self.gene_synonyms}
        keys.discard("")
        return keys


_NAME_WEIGHTS = {
    "gene_primary": 100, "recommended": 100, "recommended_short": 92, "gene_synonym": 88,
    "alternative": 82, "alternative_short": 80, "orf": 70, "submitted": 60, "accession": 100,
}


def extract_uniprot_names(rec: dict) -> List[Tuple[str, str]]:
    """All (name, kind) pairs of a UniProtKB JSON record."""
    out: List[Tuple[str, str]] = []
    pdsc = rec.get("proteinDescription") or {}
    rn = pdsc.get("recommendedName") or {}
    if (rn.get("fullName") or {}).get("value"):
        out.append((rn["fullName"]["value"], "recommended"))
    for s in rn.get("shortNames") or []:
        if s.get("value"):
            out.append((s["value"], "recommended_short"))
    for alt in pdsc.get("alternativeNames") or []:
        if (alt.get("fullName") or {}).get("value"):
            out.append((alt["fullName"]["value"], "alternative"))
        for s in alt.get("shortNames") or []:
            if s.get("value"):
                out.append((s["value"], "alternative_short"))
    for sub in pdsc.get("submissionNames") or []:
        if (sub.get("fullName") or {}).get("value"):
            out.append((sub["fullName"]["value"], "submitted"))
    for i, g in enumerate(rec.get("genes") or []):
        if (g.get("geneName") or {}).get("value"):
            out.append((g["geneName"]["value"], "gene_primary" if i == 0 else "gene_synonym"))
        for s in g.get("synonyms") or []:
            if s.get("value"):
                out.append((s["value"], "gene_synonym"))
        for s in (g.get("orfNames") or []) + (g.get("orderedLocusNames") or []):
            if s.get("value"):
                out.append((s["value"], "orf"))
    return out


def _is_reviewed(rec: dict) -> bool:
    return str(rec.get("entryType", "")).lower().startswith("uniprotkb reviewed")


def score_uniprot_record(rec: dict, qkey: str) -> Tuple[int, str]:
    best, basis = 0, ""
    for name, kind in extract_uniprot_names(rec):
        k = norm_key(name)
        if k == qkey and _NAME_WEIGHTS[kind] > best:
            best, basis = _NAME_WEIGHTS[kind], f"exact {kind.replace('_', ' ')}: {name}"
    if best == 0 and len(qkey) >= 4:
        for name, _kind in extract_uniprot_names(rec):
            if qkey in norm_key(name):
                best, basis = 20, f"partial name match: {name}"
                break
    if best and _is_reviewed(rec):
        best += 5
    return best, basis


def find_protein_synonyms(rec: dict) -> Tuple[Optional[str], Optional[str], List[str], List[str]]:
    """Return (official name, gene symbol, gene synonyms, all unique aliases) from a UniProt record."""
    names = extract_uniprot_names(rec)
    official = next((n for n, k in names if k == "recommended"), None) or \
        next((n for n, k in names if k == "submitted"), None)
    gene = next((n for n, k in names if k == "gene_primary"), None)
    gene_syn = [n for n, k in names if k == "gene_synonym"]
    aliases = list(dict.fromkeys(n for n, k in names if k != "orf"))
    return official, gene, gene_syn, aliases


def _lucene_quote(s: str) -> str:
    return '"' + s.replace("\\", " ").replace('"', " ") + '"'


def resolve_protein_identity(http: HttpClient, input_name: str,
                             forced_accession: Optional[str] = None) -> ProteinIdentity:
    normalized = normalize_protein_name(input_name)
    ident = ProteinIdentity(input_name=input_name, normalized=normalized)
    if not normalized and not forced_accession:
        ident.message = "Empty protein name."
        return ident

    qkey = norm_key(normalized)
    acc_query = forced_accession or (normalized.upper() if UNIPROT_ACC_RE.match(normalized.upper()) else None)

    if acc_query:
        queries = [f"accession:{acc_query.upper()}"]
    else:
        q = _lucene_quote(normalized)
        queries = [
            f"(organism_id:{TARGET_TAXONOMY_ID}) AND (gene_exact:{q} OR protein_name:{q} OR gene:{q})",
            f"(organism_id:{TARGET_TAXONOMY_ID}) AND ({q})",
        ]

    records: List[dict] = []
    try:
        for query in queries:
            res = http.get(UNIPROT_SEARCH_URL,
                           params={"query": query, "format": "json", "size": 25},
                           ttl=CACHE_TTL_UNIPROT_SECONDS)
            if res.status == 200:
                records = parse_json(res).get("results", [])
            if records:
                break
    except APIError as exc:
        ident.status, ident.message = "API_ERROR", f"UniProt lookup failed: {exc}"
        return ident

    if acc_query:  # accession given: take that record, still insist on target organism
        records = [r for r in records if r.get("primaryAccession", "").upper() == acc_query.upper()
                   or acc_query.upper() in [a.upper() for a in r.get("secondaryAccessions", [])]]
    records = [r for r in records
               if (r.get("organism") or {}).get("taxonId") in (TARGET_TAXONOMY_ID, None)]
    if not records:
        ident.message = f"No {TARGET_ORGANISM} UniProt entry found for '{normalized}'."
        return ident

    scored = []
    for r in records:
        s, basis = score_uniprot_record(r, qkey if not acc_query else norm_key(acc_query))
        if acc_query:
            s, basis = 100 + (5 if _is_reviewed(r) else 0), f"UniProt accession {acc_query.upper()}"
        scored.append((s, basis, r))
    scored.sort(key=lambda t: (-t[0], 0 if _is_reviewed(t[2]) else 1))
    score, basis, best = scored[0]

    official, gene, gene_syn, aliases = find_protein_synonyms(best)
    ident.status = "RESOLVED"
    ident.official_name, ident.gene_symbol, ident.gene_synonyms = official, gene, gene_syn
    ident.aliases = list(dict.fromkeys(aliases + ([normalized] if normalized else [])))
    ident.accession = best.get("primaryAccession")
    ident.secondary_accessions = list(best.get("secondaryAccessions") or [])
    ident.uniprot_entry_name = best.get("uniProtkbId")
    ident.reviewed_label = ("Reviewed (Swiss-Prot)" if _is_reviewed(best)
                            else "Unreviewed (TrEMBL)" if best.get("entryType") else None)
    org = best.get("organism") or {}
    ident.organism, ident.taxon_id = org.get("scientificName"), org.get("taxonId")
    ident.match_score, ident.match_basis = score, basis
    ident.alternatives = [
        f"{r.get('primaryAccession')} ({(extract_uniprot_names(r) or [('?', '')])[0][0]}; score {s})"
        for s, _b, r in scored[1:6]]
    if score < IDENTITY_MIN_CONFIDENT_SCORE:
        ident.message = (f"LOW-CONFIDENCE identity match (score {score}, {basis or 'no exact name match'}). "
                         f"Re-run with --uniprot <ACCESSION> if this is not the protein you meant.")
    return ident

# %% [markdown]
# ## 5. RCSB search (paginated)
# Two independent routes are merged ("synonym-merged candidates"):
# * **UniProt accession(s)** → `reference_sequence_identifiers` (most precise),
# * **aliases / gene symbols** → polymer-entity description phrase + gene name (human only).
#
# The result is only a *candidate* list. Every candidate is validated afterwards.

# %%
def _terminal(attribute: str, operator: str, value: Any) -> dict:
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attribute, "operator": operator, "value": value}}


def _group(op: str, nodes: List[dict]) -> dict:
    return nodes[0] if len(nodes) == 1 else {"type": "group", "logical_operator": op, "nodes": nodes}


def build_accession_query(accessions: List[str]) -> dict:
    nodes = [_group("and", [
        _terminal("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                  "exact_match", "UniProt"),
        _terminal("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                  "exact_match", acc)]) for acc in accessions]
    return _group("or", nodes)


def build_name_query(aliases: List[str], gene_symbols: List[str]) -> Optional[dict]:
    nodes = []
    for a in aliases:
        nodes.append(_terminal("rcsb_polymer_entity.pdbx_description", "contains_phrase", a))
    for g in gene_symbols:
        nodes.append(_terminal("rcsb_entity_source_organism.rcsb_gene_name.value", "exact_match", g))
    if not nodes:
        return None
    return _group("and", [_group("or", nodes),
                          _terminal("rcsb_entity_source_organism.ncbi_scientific_name", "exact_match",
                                    TARGET_ORGANISM)])


def _run_search(http: HttpClient, query: dict) -> List[str]:
    """Run one RCSB Search API query and return ALL matching entry IDs (every page)."""
    ids: List[str] = []
    start, total = 0, None
    while True:
        body = {"query": query, "return_type": "entry",
                "request_options": {"paginate": {"start": start, "rows": SEARCH_PAGE_SIZE},
                                    "results_verbosity": "compact"}}
        res = http.post_json(RCSB_SEARCH_URL, body, ttl=CACHE_TTL_SEARCH_SECONDS)
        if res.status == 204 or not res.content.strip():
            break                                        # 204 = valid query, zero hits
        data = parse_json(res)
        total = int(data.get("total_count", 0))
        batch = data.get("result_set") or []
        ids += [(it if isinstance(it, str) else it.get("identifier", "")).upper() for it in batch]
        start += SEARCH_PAGE_SIZE
        if not batch or start >= total:
            break
    unique = list(dict.fromkeys(i for i in ids if i))
    if total is not None and len(unique) < total:        # pagination safety net
        log.warning("Paginated search returned %d/%d IDs; retrying with return_all_hits.",
                    len(unique), total)
        body = {"query": query, "return_type": "entry",
                "request_options": {"return_all_hits": True, "results_verbosity": "compact"}}
        res = http.post_json(RCSB_SEARCH_URL, body, ttl=CACHE_TTL_SEARCH_SECONDS)
        if res.status == 200:
            data = parse_json(res)
            unique = list(dict.fromkeys(
                (it if isinstance(it, str) else it.get("identifier", "")).upper()
                for it in data.get("result_set") or []))
    return unique


@dataclass
class SearchOutcome:
    candidates: Dict[str, List[str]] = field(default_factory=dict)   # pdb_id -> routes that found it
    failed_routes: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    @property
    def failed(self) -> bool:
        return bool(self.failed_routes)


def search_rcsb(http: HttpClient, identity: ProteinIdentity) -> SearchOutcome:
    out = SearchOutcome()
    routes: List[Tuple[str, dict]] = []
    if identity.all_accessions:
        routes.append((f"UniProt:{','.join(identity.all_accessions[:5])}",
                       build_accession_query(identity.all_accessions[:5])))

    aliases = [a for a in identity.aliases if len(a) >= MIN_ALIAS_LENGTH][:MAX_ALIAS_QUERIES]
    genes = [g for g in ([identity.gene_symbol] if identity.gene_symbol else []) + identity.gene_synonyms
             if g][:MAX_ALIAS_QUERIES]
    if not identity.all_accessions and identity.normalized:   # unresolved: still try the raw input
        aliases = list(dict.fromkeys(aliases + [identity.normalized]))
    nq = build_name_query(aliases, genes)
    if nq:
        routes.append(("names/genes", nq))

    for label, query in routes:
        try:
            for pid in _run_search(http, query):
                out.candidates.setdefault(pid, []).append(label)
        except APIError as exc:
            out.failed_routes.append(label)
            out.errors.append(f"RCSB search route '{label}' failed: {exc}")
            log.warning(out.errors[-1])
    return out

# %% [markdown]
# ## 6. PDB entry processing
# For every candidate (in parallel): entry metadata → method check → polymer entities →
# **target identity validation** (UniProt accession match, or exact name/gene match when the entity has
# no UniProt mapping) → organism check → residue counts, R-values, validation, ligands, activity, chains.
#
# Method is checked first because it costs one request and avoids fetching entities for EM/NMR entries.

# %%
@dataclass
class EntityInfo:
    entity_id: str
    description: str
    polymer_type: str
    length: Optional[int]
    accessions: set
    genes: List[str]
    organisms: List[Tuple[str, Optional[int]]]
    chains: List[str]
    is_target: bool = False
    match_basis: str = ""

    @property
    def is_human(self) -> bool:
        return any((tid == TARGET_TAXONOMY_ID) or (name or "").casefold() == TARGET_ORGANISM.casefold()
                   for name, tid in self.organisms)


def extract_entity_info(ej: dict, entity_id: str) -> EntityInfo:
    rpe = ej.get("rcsb_polymer_entity") or {}
    ep = ej.get("entity_poly") or {}
    ids = ej.get("rcsb_polymer_entity_container_identifiers") or {}

    accessions = set()
    for ref in ids.get("reference_sequence_identifiers") or []:
        if str(ref.get("database_name", "")).lower() == "uniprot" and ref.get("database_accession"):
            accessions.add(str(ref["database_accession"]).split("-")[0].upper())
    for u in ids.get("uniprot_ids") or []:
        accessions.add(str(u).split("-")[0].upper())

    sources = ej.get("rcsb_entity_source_organism") or []
    organisms = [(s.get("ncbi_scientific_name"), s.get("ncbi_taxonomy_id")) for s in sources]
    genes = [g.get("value") for s in sources for g in (s.get("rcsb_gene_name") or []) if g.get("value")]

    asym, auth = ids.get("asym_ids") or [], ids.get("auth_asym_ids") or []
    if auth and len(asym) == len(auth):
        chains = [a if a == b else f"{a}[auth {b}]" for a, b in zip(asym, auth)]
    else:
        chains = list(auth or asym)

    length = ep.get("rcsb_sample_sequence_length")
    if length is None and ep.get("pdbx_seq_one_letter_code_can"):
        length = len(re.sub(r"\s+", "", ep["pdbx_seq_one_letter_code_can"]))
    return EntityInfo(
        entity_id=str(entity_id),
        description=rpe.get("pdbx_description") or "",
        polymer_type=str(ep.get("rcsb_entity_polymer_type") or ep.get("type") or ""),
        length=int(length) if length is not None else None,
        accessions=accessions, genes=genes, organisms=organisms, chains=chains)


def validate_target_identity(info: EntityInfo, identity: ProteinIdentity) -> Tuple[bool, str]:
    """
    Does THIS deposited macromolecule (not the title / description of the entry) correspond to the target?
      1. entity has UniProt mapping  -> match only if it intersects the target accession set.
      2. entity has NO UniProt mapping -> require exact normalised name match OR gene-symbol match.
    A partial name hit such as "P19INK4D CDK4/6 INHIBITOR" never matches.
    """
    if info.polymer_type and info.polymer_type.lower() != "protein":
        return False, ""
    targets = set(a.upper() for a in identity.all_accessions)
    if info.accessions:
        hit = info.accessions & targets
        return (True, f"UniProt:{sorted(hit)[0]}") if hit else (False, "")
    if norm_key(info.description) in identity.alias_keys:
        return True, f"entity name: {info.description}"
    if info.genes and {norm_key(g) for g in info.genes} & identity.gene_keys:
        return True, "gene symbol"
    return False, ""


def retrieve_pdb_metadata(http: HttpClient, pdb_id: str) -> Optional[dict]:
    res = http.get(f"{RCSB_DATA_URL}/entry/{pdb_id}", ttl=CACHE_TTL_DATA_SECONDS)
    return parse_json(res) if res.status == 200 else None


def retrieve_polymer_entity(http: HttpClient, pdb_id: str, entity_id: str) -> Optional[dict]:
    res = http.get(f"{RCSB_DATA_URL}/polymer_entity/{pdb_id}/{entity_id}", ttl=CACHE_TTL_DATA_SECONDS)
    return parse_json(res) if res.status == 200 else None


def entry_methods(entry: dict) -> List[str]:
    methods = [str(e.get("method", "")).strip().upper() for e in entry.get("exptl") or [] if e.get("method")]
    if not methods:
        m = ci_get(entry.get("rcsb_entry_info") or {}, "experimental_method")
        if m:
            methods = [str(m).upper()]
    return list(dict.fromkeys(methods))


def retrieve_residue_counts(entry: dict) -> Tuple[Optional[int], Optional[int]]:
    """'Macromolecule Content': deposited and modeled residue counts exactly as reported by RCSB."""
    info = entry.get("rcsb_entry_info") or {}
    dep = ci_get(info, "deposited_polymer_monomer_count")
    mod = ci_get(info, "deposited_modeled_polymer_monomer_count")
    return (int(dep) if dep is not None else None, int(mod) if mod is not None else None)


def _first_refine_value(entry: dict, *names: str) -> Optional[Decimal]:
    for row in entry.get("refine") or []:
        v = to_decimal(ci_get(row, *names))
        if v is not None:
            return v
    return None


# ---- minimal mmCIF header reader (only to recover the ORIGINAL text of numbers) ----
def _cif_tokens(text: str) -> Iterable[Tuple[str, str]]:
    lines = text.splitlines()
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if line.startswith(";"):                       # multi-line text field
            buf = [line[1:]]
            i += 1
            while i < n and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            i += 1
            yield ("str", "\n".join(buf))
            continue
        pos, L = 0, len(line)
        while pos < L:
            c = line[pos]
            if c.isspace():
                pos += 1
            elif c == "#":
                break
            elif c in ("'", '"'):
                j = pos + 1
                while j < L and not (line[j] == c and (j + 1 == L or line[j + 1].isspace())):
                    j += 1
                yield ("str", line[pos + 1:j])
                pos = j + 1
            else:
                j = pos
                while j < L and not line[j].isspace():
                    j += 1
                yield ("tok", line[pos:j])
                pos = j
        i += 1


def parse_cif_items(text: str, wanted: Iterable[str]) -> Dict[str, str]:
    """Return {item_name: first non-null raw text value} for the wanted mmCIF items."""
    wanted_l = {w.lower(): w for w in wanted}
    found: Dict[str, str] = {}
    tokens = list(_cif_tokens(text))
    i, n = 0, len(tokens)

    def is_keyword(t):
        return t[0] == "tok" and (t[1].startswith("_") or t[1].lower() == "loop_"
                                  or t[1].lower().startswith("data_"))

    while i < n:
        kind, val = tokens[i]
        if kind == "tok" and val.lower() == "loop_":
            i += 1
            headers = []
            while i < n and tokens[i][0] == "tok" and tokens[i][1].startswith("_"):
                headers.append(tokens[i][1].lower())
                i += 1
            start = i
            while i < n and not is_keyword(tokens[i]):
                i += 1
            if any(h in wanted_l for h in headers):
                vals, w = [t[1] for t in tokens[start:i]], len(headers)
                for r in range(0, len(vals) - w + 1, w):
                    for h, v in zip(headers, vals[r:r + w]):
                        if h in wanted_l and wanted_l[h] not in found and v not in ("?", "."):
                            found[wanted_l[h]] = v
        elif kind == "tok" and val.startswith("_"):
            if i + 1 < n and val.lower() in wanted_l and wanted_l[val.lower()] not in found:
                v = tokens[i + 1][1]
                if v not in ("?", "."):
                    found[wanted_l[val.lower()]] = v
            i += 2
        else:
            i += 1
    return found


CIF_ITEMS = ["_refine.ls_d_res_high", "_reflns.d_resolution_high",
             "_refine.ls_R_factor_R_free", "_refine.ls_R_factor_R_work"]


def retrieve_cif_source_text(http: HttpClient, pdb_id: str) -> Dict[str, str]:
    """Original text of resolution / depositor R-values from the (small) mmCIF header. Cached."""
    cached = http.cache.get_json(f"cif|{pdb_id}", CACHE_TTL_DATA_SECONDS)
    if cached is not None:
        return cached
    res = http.get(RCSB_HEADER_CIF_URL.format(pdb_id=pdb_id), ttl=None)
    values = parse_cif_items(res.text, CIF_ITEMS) if res.status == 200 else {}
    http.cache.put_json(f"cif|{pdb_id}", values)
    return values


def retrieve_depositor_r_values(entry: dict, cif: Dict[str, str]) -> Tuple[Optional[Decimal], Optional[Decimal]]:
    """DEPOSITOR R-free / R-work (refine.ls_R_factor_*). Never the DCC / wwPDB recalculated values."""
    rfree = choose_source_value(_first_refine_value(entry, "ls_R_factor_R_free"),
                                to_decimal(cif.get("_refine.ls_R_factor_R_free")))
    rwork = choose_source_value(_first_refine_value(entry, "ls_R_factor_R_work"),
                                to_decimal(cif.get("_refine.ls_R_factor_R_work")))
    return rfree, rwork


def retrieve_resolution(entry: dict, cif: Dict[str, str]) -> Optional[Decimal]:
    res_list = ci_get(entry.get("rcsb_entry_info") or {}, "resolution_combined") or []
    api_val = to_decimal(res_list[0]) if res_list else None
    text_val = to_decimal(cif.get("_refine.ls_d_res_high")) or to_decimal(cif.get("_reflns.d_resolution_high"))
    return choose_source_value(api_val, text_val)

# %% [markdown]
# ## 7. wwPDB validation data
# Priority: (1) RCSB Data API `pdbx_vrpt_summary*` (official wwPDB validation data, structured),
# (2) official `*_validation.xml.gz`, (3) validation PDF text for *Average B*. **No OCR is used.**

# %%
VALIDATION_FIELDS = {
    "ww_rfree": WWPDB_RFREE_FIELDS,
    "clashscore": ["clashscore"],
    "rama": ["percent_ramachandran_outliers", "percent_rama_outliers"],
    "sidechain": ["percent_rotamer_outliers", "percent_rota_outliers"],
    "rsrz": ["percent_RSRZ_outliers"],
}


def _vrpt_blocks(entry: dict) -> List[dict]:
    blocks: List[dict] = []
    for k, v in entry.items():
        if k.lower().startswith("pdbx_vrpt_summary"):
            if isinstance(v, dict):
                blocks.append(v)
            elif isinstance(v, list):
                blocks += [x for x in v if isinstance(x, dict)]
    return blocks


def _validation_from_entry(entry: dict) -> Dict[str, Optional[Decimal]]:
    blocks = _vrpt_blocks(entry)
    out: Dict[str, Optional[Decimal]] = {}
    for key, names in VALIDATION_FIELDS.items():
        val = None
        for b in blocks:
            val = to_decimal(ci_get(b, *names))
            if val is not None:
                break
        out[key] = val
    return out


def _validation_from_xml(http: HttpClient, pdb_id: str) -> Dict[str, Optional[Decimal]]:
    """Fallback: official validation XML (only the <Entry> element is read)."""
    cached = http.cache.get_json(f"vxml|{pdb_id}", CACHE_TTL_DATA_SECONDS)
    if cached is None:
        cached = {}
        pid = pdb_id.lower()
        for tmpl in VALIDATION_XML_URLS:
            res = http.get(tmpl.format(mid=pid[1:3], pid=pid), ttl=None, timeout=60)
            if res.status != 200:
                continue
            raw = res.content
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            for _ev, el in ET.iterparse(io.BytesIO(raw), events=("start",)):
                if el.tag.split("}")[-1] == "Entry":
                    cached = dict(el.attrib)
                    break
            break
        http.cache.put_json(f"vxml|{pdb_id}", cached)
    out: Dict[str, Optional[Decimal]] = {}
    for key, names in VALIDATION_FIELDS.items():
        out[key] = to_decimal(ci_get(cached, *names))
    return out


def retrieve_validation_data(http: HttpClient, entry: dict, pdb_id: str) -> Dict[str, Optional[Decimal]]:
    data = _validation_from_entry(entry)
    if any(v is None for v in data.values()):
        try:
            xml_data = _validation_from_xml(http, pdb_id)
            for k, v in data.items():
                if v is None:
                    data[k] = xml_data.get(k)
        except (APIError, ET.ParseError, OSError) as exc:
            log.debug("Validation XML fallback failed for %s: %s", pdb_id, exc)
    return data


AVG_B_RE = re.compile(
    r"Average\s*B\s*,?\s*all\s*atoms\s*\([^)]{1,15}\)\s*[:\-]?\s*([0-9]+(?:\.[0-9]+)?)", re.I)


def _average_b_from_pdf(data: bytes) -> Optional[str]:
    from pypdf import PdfReader  # imported lazily: only needed for this field
    reader = PdfReader(io.BytesIO(data))
    for i, page in enumerate(reader.pages):
        if i >= PDF_MAX_PAGES_SCANNED:
            break
        try:
            text = page.extract_text() or ""
        except Exception:
            continue
        m = AVG_B_RE.search(text)
        if m:
            return m.group(1)
    return None


def retrieve_average_b(http: HttpClient, pdb_id: str) -> Optional[Decimal]:
    """'Average B, all atoms (A^2)' from section 4 of the wwPDB full validation report (kept as text)."""
    cached = http.cache.get_json(f"avgb|{pdb_id}", CACHE_TTL_DATA_SECONDS)
    if cached is not None:
        return to_decimal(cached.get("avg_b"))
    pid = pdb_id.lower()
    value: Optional[str] = None
    for tmpl in VALIDATION_PDF_URLS:
        res = http.get(tmpl.format(mid=pid[1:3], pid=pid), ttl=None, timeout=90)
        if res.status != 200:
            continue
        raw = res.content
        if raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
        value = _average_b_from_pdf(raw)
        break
    http.cache.put_json(f"avgb|{pdb_id}", {"avg_b": value})
    return to_decimal(value)

# %% [markdown]
# ## 8. Ligands, activity, structure type, chains

# %%
@dataclass
class Ligand:
    comp_id: str
    name: Optional[str]


class ChemCompStore:
    """Thread-safe, de-duplicated chemical-component lookups (shared across all entries)."""

    def __init__(self, http: HttpClient):
        self.http = http
        self._data: Dict[str, dict] = {}
        self._locks: Dict[str, threading.Lock] = defaultdict(threading.Lock)
        self._guard = threading.Lock()

    def get(self, comp_id: str) -> dict:
        with self._guard:
            lock = self._locks[comp_id]
        with lock:
            if comp_id not in self._data:
                res = self.http.get(f"{RCSB_DATA_URL}/chemcomp/{comp_id}", ttl=CACHE_TTL_DATA_SECONDS)
                self._data[comp_id] = parse_json(res) if res.status == 200 else {}
            return self._data[comp_id]


def heavy_atom_count(chem: dict) -> Optional[int]:
    v = ci_get(ci_get(chem, "rcsb_chem_comp_info") or {}, "atom_count_heavy")
    if v is not None:
        return int(v)
    formula = ci_get(ci_get(chem, "chem_comp") or {}, "formula")
    if not formula:
        return None
    total = 0
    for tok in str(formula).split():
        m = re.fullmatch(r"([A-Z][a-z]?)(\d*)", tok)
        if m and m.group(1) not in ("H", "D"):
            total += int(m.group(2) or 1)
    return total or None


def is_relevant_ligand(comp_id: str, chem: dict) -> bool:
    """Exclude water, ions, salts, buffers, solvents, crystallisation agents, additives, tiny components."""
    if comp_id.upper() in EXCLUDED_COMPONENT_IDS:
        return False
    heavy = heavy_atom_count(chem)
    return not (heavy is not None and heavy < MIN_LIGAND_HEAVY_ATOMS)


def retrieve_ligands(http: HttpClient, chem_store: ChemCompStore, pdb_id: str, entry: dict) -> List[Ligand]:
    ids = ci_get(entry.get("rcsb_entry_container_identifiers") or {}, "non_polymer_entity_ids") or []
    ligands: Dict[str, Ligand] = {}
    for eid in ids:
        res = http.get(f"{RCSB_DATA_URL}/nonpolymer_entity/{pdb_id}/{eid}", ttl=CACHE_TTL_DATA_SECONDS)
        if res.status != 200:
            continue
        nj = parse_json(res)
        comp_id = (ci_get(ci_get(nj, "rcsb_nonpolymer_entity_container_identifiers") or {}, "nonpolymer_comp_id")
                   or ci_get(ci_get(nj, "pdbx_entity_nonpoly") or {}, "comp_id"))
        if not comp_id or comp_id in ligands:
            continue
        chem = chem_store.get(comp_id)
        if not is_relevant_ligand(comp_id, chem):
            continue
        name = (ci_get(ci_get(chem, "chem_comp") or {}, "name")
                or ci_get(ci_get(nj, "pdbx_entity_nonpoly") or {}, "name")
                or ci_get(ci_get(nj, "rcsb_nonpolymer_entity") or {}, "pdbx_description"))
        ligands[comp_id] = Ligand(comp_id, str(name) if name else None)
    return list(ligands.values())


def retrieve_ligand_activity(entry: dict, ligands: List[Ligand]) -> str:
    """Binding-affinity annotations (RCSB, e.g. BindingDB) aggregated per ligand/type/unit/source."""
    wanted = {l.comp_id.upper() for l in ligands}
    groups: Dict[Tuple[str, str, str, str], List[Tuple[Decimal, str]]] = defaultdict(list)
    order: List[Tuple[str, str, str, str]] = []
    for row in entry.get("rcsb_binding_affinity") or []:
        comp = str(ci_get(row, "comp_id") or "").upper()
        val = to_decimal(ci_get(row, "value"))
        if comp not in wanted or val is None:
            continue
        source = str(ci_get(row, "provenance_code") or "")
        if AFFINITY_SOURCES and not any(s.lower() in source.lower() for s in AFFINITY_SOURCES):
            continue
        key = (comp, str(ci_get(row, "type") or "affinity"), str(ci_get(row, "unit") or ""), source)
        if key not in groups:
            order.append(key)
        groups[key].append((val, str(ci_get(row, "symbol") or "=")))
    parts = []
    for key in sorted(order, key=lambda k: ([l.comp_id.upper() for l in ligands].index(k[0]), k[1])):
        comp, typ, unit, source = key
        vals = groups[key]
        n = len(vals)
        src = f" [{source}]" if source else ""
        unit_s = f" {unit}" if unit else ""
        if n == 1:
            v, sym = vals[0]
            q = "" if sym in ("=", "") else sym
            parts.append(f"{comp}: {typ}: {q}{fmt_decimal(v)}{unit_s} from 1 assay{src}")
        else:
            nums = [v for v, _s in vals]
            lo, hi = min(nums), max(nums)
            quals = sorted({s for _v, s in vals if s not in ("=", "")})
            qnote = f" (includes qualifier {', '.join(quals)})" if quals else ""
            rng = (f"{fmt_decimal(lo)}{unit_s}" if lo == hi else
                   f"min {fmt_decimal(lo)}, max {fmt_decimal(hi)}{unit_s}")
            parts.append(f"{comp}: {typ}: {rng} from {n} assays{qnote}{src}")
    return "; ".join(parts) if parts else MISSING_VALUE


def classify_structure_type(has_ligand: bool, other_entities: List[EntityInfo]) -> str:
    if has_ligand:
        return "Protein-Ligand Complex"
    types = {e.polymer_type.lower() for e in other_entities}
    if "dna" in types or "na-hybrid" in types:
        return "Protein-DNA Complex"
    if "rna" in types:
        return "Protein-RNA Complex"
    peptides = [e for e in other_entities if e.polymer_type.lower() == "protein"
                and e.length is not None and e.length <= PEPTIDE_MAX_LENGTH]
    if peptides:
        return "Protein-Peptide Complex"
    if any(e.polymer_type.lower() == "protein" for e in other_entities):
        return "Protein-Protein Complex"
    return "Apo"


def retrieve_pdb_chains(target_entities: List[EntityInfo]) -> str:
    chains = []
    for e in target_entities:
        chains += e.chains
    return ", ".join(dict.fromkeys(chains)) if chains else MISSING_VALUE

# %% [markdown]
# ## 9. Per-entry pipeline + filtering

# %%
@dataclass
class PDBRecord:
    protein_name: str
    pdb_id: str
    method: str
    resolution: Optional[Decimal] = None
    deposited: Optional[int] = None
    modeled: Optional[int] = None
    residue_diff: Optional[int] = None
    r_free: Optional[Decimal] = None
    r_work: Optional[Decimal] = None
    r_diff: Optional[Decimal] = None
    ww_rfree: Optional[Decimal] = None
    clashscore: Optional[Decimal] = None
    rama: Optional[Decimal] = None
    sidechain: Optional[Decimal] = None
    rsrz: Optional[Decimal] = None
    avg_b: Optional[Decimal] = None
    ligands: List[Ligand] = field(default_factory=list)
    activity: str = MISSING_VALUE
    structure_type: str = ""
    chains: str = MISSING_VALUE
    match_basis: str = ""
    warnings: List[str] = field(default_factory=list)


@dataclass
class EntryOutcome:
    pdb_id: str
    status: str                     # QUALIFIED | EXCLUDED | ERROR
    category: str = ""              # method | identity | organism | not_found | api_error
    reason: str = ""
    record: Optional[PDBRecord] = None


@dataclass
class PipelineContext:
    http: HttpClient
    chem_store: ChemCompStore
    identity: ProteinIdentity
    skip_average_b: bool = False
    skip_cif_precision: bool = False


def process_pdb_entry(ctx: PipelineContext, pdb_id: str) -> EntryOutcome:
    http = ctx.http
    try:
        entry = retrieve_pdb_metadata(http, pdb_id)
        if entry is None:
            return EntryOutcome(pdb_id, "EXCLUDED", "not_found",
                                "Entry not available in the RCSB Data API (possibly obsolete/withdrawn).")

        methods = entry_methods(entry)
        if EXPERIMENTAL_METHOD not in methods:
            return EntryOutcome(pdb_id, "EXCLUDED", "method",
                                f"Experimental method is {', '.join(methods) or 'unknown'}, "
                                f"not {EXPERIMENTAL_METHOD}.")

        entity_ids = ci_get(entry.get("rcsb_entry_container_identifiers") or {}, "polymer_entity_ids") or []
        entities: List[EntityInfo] = []
        for eid in entity_ids:
            ej = retrieve_polymer_entity(http, pdb_id, str(eid))
            if ej is not None:
                entities.append(extract_entity_info(ej, str(eid)))
        for e in entities:
            e.is_target, e.match_basis = validate_target_identity(e, ctx.identity)

        targets = [e for e in entities if e.is_target]
        if not targets:
            names = "; ".join(f"entity {e.entity_id}: {e.description or '?'}" for e in entities) or "no polymer entities"
            return EntryOutcome(pdb_id, "EXCLUDED", "identity",
                                f"No deposited macromolecule matches the target protein ({names}).")
        human_targets = [e for e in targets if e.is_human]
        if not human_targets:
            orgs = sorted({n or "unknown" for e in targets for n, _t in e.organisms}) or ["unknown"]
            return EntryOutcome(pdb_id, "EXCLUDED", "organism",
                                f"Target macromolecule is not from {TARGET_ORGANISM} (organism: {', '.join(orgs)}).")

        rec = PDBRecord(protein_name=ctx.identity.input_name, pdb_id=pdb_id,
                        method="; ".join(methods),
                        match_basis="; ".join(dict.fromkeys(e.match_basis for e in human_targets)))
        rec.deposited, rec.modeled = retrieve_residue_counts(entry)

        cif: Dict[str, str] = {}
        if not ctx.skip_cif_precision:
            try:
                cif = retrieve_cif_source_text(http, pdb_id)
            except APIError as exc:
                rec.warnings.append(f"mmCIF header unavailable ({exc}); source precision may be reduced.")
        rec.resolution = retrieve_resolution(entry, cif)
        rec.r_free, rec.r_work = retrieve_depositor_r_values(entry, cif)

        val = retrieve_validation_data(http, entry, pdb_id)
        rec.ww_rfree, rec.clashscore, rec.rama = val["ww_rfree"], val["clashscore"], val["rama"]
        rec.sidechain, rec.rsrz = val["sidechain"], val["rsrz"]
        if all(v is None for v in val.values()):
            rec.warnings.append("Validation data unavailable.")

        if not ctx.skip_average_b:
            try:
                rec.avg_b = retrieve_average_b(http, pdb_id)
            except ImportError:
                rec.warnings.append("pypdf not installed: Average B not collected.")
            except Exception as exc:
                rec.warnings.append(f"Average B could not be extracted ({type(exc).__name__}: {exc}).")
            if rec.avg_b is None and not any("Average B" in w for w in rec.warnings):
                rec.warnings.append("Average B not found in validation report.")

        rec.ligands = retrieve_ligands(http, ctx.chem_store, pdb_id, entry)
        rec.activity = retrieve_ligand_activity(entry, rec.ligands)
        rec.structure_type = classify_structure_type(
            bool(rec.ligands), [e for e in entities if not e.is_target])
        rec.chains = retrieve_pdb_chains(human_targets)
        return EntryOutcome(pdb_id, "QUALIFIED", record=rec)

    except APIError as exc:
        return EntryOutcome(pdb_id, "ERROR", "api_error",
                            f"Unable to retrieve data because of an API/network error: {exc}")
    except Exception as exc:  # one bad entry must never stop the run
        log.exception("Unexpected error for %s", pdb_id)
        return EntryOutcome(pdb_id, "ERROR", "api_error", f"Unexpected error: {type(exc).__name__}: {exc}")


def finalize_calculated_fields(rec: PDBRecord) -> None:
    """Programmatic (never typed-in) derived fields; Decimal arithmetic avoids float artefacts."""
    if rec.deposited is not None and rec.modeled is not None:
        rec.residue_diff = int(rec.deposited) - int(rec.modeled)
    if rec.r_free is not None and rec.r_work is not None:
        rec.r_diff = rec.r_free - rec.r_work


@dataclass
class PDBResult:
    state: str = "NOT_RUN"          # OK | NO_STRUCTURES | API_ERROR | PARTIAL
    records: List[PDBRecord] = field(default_factory=list)
    excluded: List[EntryOutcome] = field(default_factory=list)
    errors: List[EntryOutcome] = field(default_factory=list)
    message: str = ""
    candidate_count: int = 0
    funnel: Dict[str, int] = field(default_factory=dict)
    search_errors: List[str] = field(default_factory=list)


def filter_pdb_structures(outcomes: List[EntryOutcome], candidate_count: int) -> PDBResult:
    res = PDBResult(candidate_count=candidate_count)
    for o in outcomes:
        if o.status == "QUALIFIED" and o.record:
            finalize_calculated_fields(o.record)
            res.records.append(o.record)
        elif o.status == "EXCLUDED":
            res.excluded.append(o)
        else:
            res.errors.append(o)
    res.records.sort(key=lambda r: r.pdb_id)
    c = lambda cat: sum(1 for o in res.excluded if o.category == cat)  # noqa: E731
    xray = candidate_count - len(res.errors) - c("not_found") - c("method")
    ident_ok = xray - c("identity")
    res.funnel = {
        "Initial search candidates": candidate_count,
        "Entries that could not be retrieved (API error)": len(res.errors),
        "Entries not found / obsolete": c("not_found"),
        "Excluded: method is not X-RAY DIFFRACTION": c("method"),
        "X-ray entries evaluated": xray,
        "Excluded: target identity not validated": c("identity"),
        "Target identity validated": ident_ok,
        f"Excluded: target not {TARGET_ORGANISM}": c("organism"),
        "Final qualifying structures": len(res.records),
    }
    return res

# %% [markdown]
# ## 10. AlphaFold (always searched, independent of PDB)

# %%
@dataclass
class AFRecord:
    protein_name: str
    af_id: str
    entry_type: str
    source: str
    status: str
    plddt_text: str
    accession: str = ""


@dataclass
class AFResult:
    state: str = "NOT_RUN"          # OK | NO_STRUCTURES | API_ERROR | NO_IDENTITY
    records: List[AFRecord] = field(default_factory=list)
    candidate_count: int = 0
    rejected: List[str] = field(default_factory=list)
    message: str = ""


def plddt_category(value: Decimal) -> str:
    for threshold, label in PLDDT_BANDS:
        if value >= threshold:
            return label
    return PLDDT_BANDS[-1][1]


def search_alphafold(http: HttpClient, identity: ProteinIdentity) -> Tuple[List[dict], str, str]:
    """
    Query the official AlphaFold DB API with the resolved UniProt accession.
    Returns (raw entries, state, message) with state in OK | NO_STRUCTURES | API_ERROR | NO_IDENTITY.
    """
    if identity.status == "API_ERROR":
        return [], "API_ERROR", identity.message
    if not identity.accession:
        return [], "NO_IDENTITY", f"No {TARGET_ORGANISM} UniProt accession could be resolved."
    try:
        res = http.get(ALPHAFOLD_PREDICTION_URL.format(accession=identity.accession),
                       ttl=CACHE_TTL_DATA_SECONDS)
    except APIError as exc:
        return [], "API_ERROR", f"Unable to retrieve data because of an API/network error: {exc}"
    if res.status in (204, 404):
        return [], "NO_STRUCTURES", "AlphaFold DB has no entry for this accession."
    data = parse_json(res)
    if isinstance(data, dict):
        data = [data]
    return list(data or []), ("OK" if data else "NO_STRUCTURES"), ""


def retrieve_alphafold_data(raw_entries: List[dict], identity: ProteinIdentity) -> AFResult:
    """Validate identity + organism for every entry and map to output rows."""
    out = AFResult(candidate_count=len(raw_entries))
    targets = {a.upper() for a in identity.all_accessions}
    for e in raw_entries:
        af_id = ci_get(e, "modelEntityId", "entryId")      # new name first, deprecated name second
        acc = str(ci_get(e, "uniprotAccession") or "").split("-")[0].upper()
        if not af_id:
            out.rejected.append("entry without an AlphaFold identifier")
            continue
        if acc and acc not in targets:
            out.rejected.append(f"{af_id}: UniProt accession {acc} is not the target")
            continue
        tax = ci_get(e, "taxId")
        org = ci_get(e, "organismScientificName")
        if (tax is not None and int(tax) != TARGET_TAXONOMY_ID) or \
           (tax is None and org and str(org).casefold() != TARGET_ORGANISM.casefold()):
            out.rejected.append(f"{af_id}: organism is {org or tax}, not {TARGET_ORGANISM}")
            continue

        explicit_type = ci_get(e, "entryType", "complexType", "modelType")
        tool = str(ci_get(e, "toolUsed") or "")
        if explicit_type:
            entry_type = str(explicit_type)
        elif "multimer" in tool.lower():
            entry_type = "Multimer"
        elif "monomer" in tool.lower():
            entry_type = "Monomer"
        else:
            entry_type = MISSING_VALUE

        source_parts = [str(p) for p in (ci_get(e, "providerId"), tool) if p]
        source = " | ".join(source_parts) if source_parts else MISSING_VALUE

        reviewed = ci_get(e, "isUniProtReviewed", "isReviewed")
        if reviewed is not None:
            status = "Reviewed (Swiss-Prot)" if reviewed else "Unreviewed (TrEMBL)"
        else:
            status = identity.reviewed_label or MISSING_VALUE

        plddt = to_decimal(ci_get(e, "globalMetricValue"))
        plddt_text = f"{fmt_decimal(plddt)} ({plddt_category(plddt)})" if plddt is not None else MISSING_VALUE
        out.records.append(AFRecord(identity.input_name, str(af_id), entry_type, source, status,
                                    plddt_text, acc))
    out.records.sort(key=lambda r: r.af_id)
    out.state = "OK" if out.records else "NO_STRUCTURES"
    return out

# %% [markdown]
# ## 11. Consistency checks before export

# %%
def validate_records(records: List[PDBRecord], identity: ProteinIdentity) -> List[str]:
    issues: List[str] = []
    tol = Decimal("0.0000001")
    for r in records:
        p = r.pdb_id
        if not PDB_ID_RE.match(p):
            issues.append(f"{p}: invalid PDB identifier format")
        if r.deposited is not None and r.modeled is not None:
            if r.deposited < r.modeled:
                issues.append(f"{p}: deposited residue count ({r.deposited}) < modeled ({r.modeled})")
            if r.residue_diff != r.deposited - r.modeled:
                issues.append(f"{p}: residue count difference inconsistent")
        if r.r_free is not None and r.r_work is not None:
            if r.r_diff is None or abs((r.r_free - r.r_work) - r.r_diff) > tol:
                issues.append(f"{p}: R-value difference inconsistent")
        if r.method and EXPERIMENTAL_METHOD not in r.method:
            issues.append(f"{p}: experimental method {r.method} does not include {EXPERIMENTAL_METHOD}")
        if not r.chains or r.chains == MISSING_VALUE:
            issues.append(f"{p}: no target chains identified (identity check suspect)")
        if not r.match_basis:
            issues.append(f"{p}: no identity-match basis recorded")
        for lig in r.ligands:
            if not lig.comp_id:
                issues.append(f"{p}: ligand without ID")
        if r.ligands and r.structure_type != "Protein-Ligand Complex":
            issues.append(f"{p}: ligands present but structure type is {r.structure_type}")
        if not r.ligands and r.structure_type == "Protein-Ligand Complex":
            issues.append(f"{p}: structure type is Protein-Ligand Complex but no ligand recorded")
    return issues

# %% [markdown]
# ## 12. Excel workbook (openpyxl)
# Source values are written as numbers with an Excel number format that reproduces the source
# precision (`0.20` shows as `0.20`, `28.0` as `28.0`). Calculated columns are computed in Python.

# %%
PDB_COLUMNS = [
    "Protein Name", "PDB Hyperlink", "Experimental Method", "Resolution (Å)", "Deposited Residue Count",
    "Modeled Residue Count", "Residue Count Difference", "R-Value Free (Depositor)",
    "R-Value Work (Depositor)", "R-Value Difference", "wwPDB R-Free", "Clashscore",
    "Ramachandran Outliers (%)", "Sidechain Outliers (%)", "RSRZ Outliers (%)",
    "Average B, all atoms (Å²)", "Ligand ID", "Ligand Name", "Ligand Activity", "Structure Type",
    "PDB Chains",
]
PDB_WIDTHS = [30, 12, 20, 13, 14, 14, 14, 14, 14, 13, 12, 11, 15, 15, 14, 16, 18, 48, 52, 22, 22]
AF_COLUMNS = ["Protein Name", "AlphaFold ID", "Entry Type", "Source", "Status", "Average pLDDT"]
AF_WIDTHS = [30, 28, 16, 46, 24, 22]
HEADER_ROW = 4

_HDR_FONT = Font(bold=True, color="FFFFFF")
_HDR_FILL = PatternFill("solid", fgColor="305496")
_WRAP = Alignment(wrap_text=True, vertical="top")


def _write_header(ws, columns: List[str], widths: List[int]) -> None:
    for c, (name, w) in enumerate(zip(columns, widths), start=1):
        cell = ws.cell(row=HEADER_ROW, column=c, value=name)
        cell.font, cell.fill = _HDR_FONT, _HDR_FILL
        cell.alignment = Alignment(wrap_text=True, vertical="center", horizontal="center")
        ws.column_dimensions[get_column_letter(c)].width = w
    ws.row_dimensions[HEADER_ROW].height = 32


def _put_number(cell, value: Optional[Decimal], integer: bool = False) -> None:
    if value is None:
        cell.value = MISSING_VALUE
        return
    if integer or decimals_of(value) == 0 and value == value.to_integral_value():
        cell.value = int(value)
        cell.number_format = "0"
    else:
        cell.value = float(value)
        cell.number_format = excel_number_format(value)


def _put_text(cell, text: Optional[str]) -> None:
    cell.value = text if text not in (None, "") else MISSING_VALUE
    cell.alignment = _WRAP


def write_pdb_sheet(ws, protein_name: str, res: PDBResult) -> None:
    ws["A1"] = f"Total qualifying X-ray crystal structures found for {protein_name}: {len(res.records)}"
    ws["A1"].font = Font(bold=True, size=12)
    ws["A2"] = res.message or ("Filters: identity-validated target macromolecule, "
                               f"{TARGET_ORGANISM}, {EXPERIMENTAL_METHOD}.")
    ws["A2"].font = Font(italic=True)
    _write_header(ws, PDB_COLUMNS, PDB_WIDTHS)

    row = HEADER_ROW + 1
    if not res.records:
        ws.cell(row=row, column=1, value=res.message or "No qualifying crystal structure was found.")
    for r in res.records:
        ws.cell(row=row, column=1, value=r.protein_name).alignment = _WRAP
        link = ws.cell(row=row, column=2, value=r.pdb_id)
        link.hyperlink = RCSB_STRUCTURE_URL.format(pdb_id=r.pdb_id)
        link.style = "Hyperlink"
        _put_text(ws.cell(row=row, column=3), r.method)
        _put_number(ws.cell(row=row, column=4), r.resolution)
        _put_number(ws.cell(row=row, column=5), None if r.deposited is None else Decimal(r.deposited), True)
        _put_number(ws.cell(row=row, column=6), None if r.modeled is None else Decimal(r.modeled), True)
        _put_number(ws.cell(row=row, column=7), None if r.residue_diff is None else Decimal(r.residue_diff), True)
        _put_number(ws.cell(row=row, column=8), r.r_free)
        _put_number(ws.cell(row=row, column=9), r.r_work)
        _put_number(ws.cell(row=row, column=10), r.r_diff)
        _put_number(ws.cell(row=row, column=11), r.ww_rfree)
        _put_number(ws.cell(row=row, column=12), r.clashscore)
        _put_number(ws.cell(row=row, column=13), r.rama)
        _put_number(ws.cell(row=row, column=14), r.sidechain)
        _put_number(ws.cell(row=row, column=15), r.rsrz)
        _put_number(ws.cell(row=row, column=16), r.avg_b)
        _put_text(ws.cell(row=row, column=17), "; ".join(l.comp_id for l in r.ligands) or None)
        _put_text(ws.cell(row=row, column=18),
                  "; ".join(f"{l.comp_id}: {l.name or MISSING_VALUE}" for l in r.ligands) or None)
        _put_text(ws.cell(row=row, column=19), r.activity)
        _put_text(ws.cell(row=row, column=20), r.structure_type)
        _put_text(ws.cell(row=row, column=21), r.chains)
        row += 1
    ws.freeze_panes = ws.cell(row=HEADER_ROW + 1, column=3)
    ws.auto_filter.ref = f"A{HEADER_ROW}:{get_column_letter(len(PDB_COLUMNS))}{max(row - 1, HEADER_ROW + 1)}"


def write_af_sheet(ws, protein_name: str, res: AFResult) -> None:
    ws["A1"] = f"AlphaFold entries found for {protein_name}: {len(res.records)}"
    ws["A1"].font = Font(bold=True, size=12)
    ws["A2"] = res.message or f"Filters: identity-validated UniProt accession, {TARGET_ORGANISM}."
    ws["A2"].font = Font(italic=True)
    _write_header(ws, AF_COLUMNS, AF_WIDTHS)
    row = HEADER_ROW + 1
    if not res.records:
        ws.cell(row=row, column=1, value=res.message or "No AlphaFold structure is available.")
    for r in res.records:
        for c, v in enumerate([r.protein_name, r.af_id, r.entry_type, r.source, r.status, r.plddt_text], start=1):
            ws.cell(row=row, column=c, value=v).alignment = _WRAP
        row += 1
    ws.freeze_panes = ws.cell(row=HEADER_ROW + 1, column=1)
    ws.auto_filter.ref = f"A{HEADER_ROW}:{get_column_letter(len(AF_COLUMNS))}{max(row - 1, HEADER_ROW + 1)}"


def create_excel_workbook(path: Path, identity: ProteinIdentity, pdb: PDBResult, af: AFResult,
                          issues: List[str], retrieved_at: datetime) -> Path:
    wb = Workbook()
    ws_pdb = wb.active
    ws_pdb.title = "PDB_Data"
    write_pdb_sheet(ws_pdb, identity.input_name, pdb)
    write_af_sheet(wb.create_sheet("AlphaFold_Data"), identity.input_name, af)

    ws = wb.create_sheet("Summary")
    rows = [
        ("Original protein input", identity.input_name),
        ("Resolved protein name", identity.official_name or MISSING_VALUE),
        ("Gene symbol", identity.gene_symbol or MISSING_VALUE),
        ("UniProt accession", identity.accession or MISSING_VALUE),
        ("UniProt status", identity.reviewed_label or MISSING_VALUE),
        ("Identity match basis", f"{identity.match_basis} (score {identity.match_score})"
         if identity.match_basis else MISSING_VALUE),
        ("Identity warning", identity.message or "none"),
        ("Aliases used", "; ".join(identity.aliases) or MISSING_VALUE),
        ("Organism filter", TARGET_ORGANISM),
        ("Experimental method filter", EXPERIMENTAL_METHOD),
        ("PDB state", pdb.state),
        ("PDB candidate count (initial search)", pdb.candidate_count),
        ("Final qualifying PDB count", len(pdb.records)),
        ("AlphaFold state", af.state),
        ("AlphaFold candidate count", af.candidate_count),
        ("AlphaFold result count", len(af.records)),
        ("Date/time of retrieval", retrieved_at.strftime("%Y-%m-%d %H:%M:%S")),
        ("", ""),
        ("PDB filtering funnel", ""),
        *[(f"  {k}", v) for k, v in pdb.funnel.items()],
    ]
    for i, (k, v) in enumerate(rows, start=1):
        ws.cell(row=i, column=1, value=k).font = Font(bold=True)
        ws.cell(row=i, column=2, value=v).alignment = _WRAP
    ws.column_dimensions["A"].width = 42
    ws.column_dimensions["B"].width = 100

    ws = wb.create_sheet("Excluded_PDB")
    for c, h in enumerate(["PDB ID", "Outcome", "Reason"], start=1):
        cell = ws.cell(row=1, column=c, value=h)
        cell.font, cell.fill = _HDR_FONT, _HDR_FILL
    r = 2
    for o in sorted(pdb.excluded + pdb.errors, key=lambda o: o.pdb_id):
        outcome = ("API/network error" if o.status == "ERROR"
                   else "Structure found but excluded (did not satisfy identity/organism/method criteria)")
        for c, v in enumerate([o.pdb_id, outcome, o.reason], start=1):
            ws.cell(row=r, column=c, value=v).alignment = _WRAP
        r += 1
    for col, w in zip("ABC", (12, 60, 110)):
        ws.column_dimensions[col].width = w
    ws.freeze_panes = "A2"

    ws = wb.create_sheet("Warnings")
    for c, h in enumerate(["Scope", "Message"], start=1):
        cell = ws.cell(row=1, column=c, value=h)
        cell.font, cell.fill = _HDR_FONT, _HDR_FILL
    r = 2
    msgs = [("Identity", identity.message)] if identity.message else []
    msgs += [("PDB search", m) for m in pdb.search_errors]
    msgs += [(rec.pdb_id, w) for rec in pdb.records for w in rec.warnings]
    msgs += [("AlphaFold", m) for m in af.rejected]
    msgs += [("Consistency check", m) for m in issues]
    for scope, msg in msgs or [("-", "No warnings.")]:
        ws.cell(row=r, column=1, value=scope)
        ws.cell(row=r, column=2, value=msg).alignment = _WRAP
        r += 1
    ws.column_dimensions["A"].width = 20
    ws.column_dimensions["B"].width = 130

    try:
        wb.save(path)
    except PermissionError:   # file open in Excel
        path = path.with_name(f"{path.stem}_{datetime.now():%Y%m%d_%H%M%S}{path.suffix}")
        wb.save(path)
    return path

# %% [markdown]
# ## 13. Orchestration
# `run_pipeline()` runs the PDB branch and the AlphaFold branch **independently**.
# AlphaFold is never conditional on the PDB result.

# %%
@dataclass
class RuntimeOptions:
    outdir: Path = Path(".")
    workers: int = MAX_WORKERS
    use_cache: bool = True
    refresh_cache: bool = False
    skip_average_b: bool = False
    skip_cif_precision: bool = False
    uniprot_override: Optional[str] = None
    verbose: bool = False


def _banner(title: str) -> None:
    print("=" * 50 + f"\n{title}\n" + "=" * 50)


def _section(title: str) -> None:
    print("\n" + "-" * 50 + f"\n{title}\n" + "-" * 50)


def run_pdb_branch(ctx: PipelineContext, opts: RuntimeOptions) -> PDBResult:
    identity = ctx.identity
    _section("RCSB PDB SEARCH")
    if identity.status == "API_ERROR":
        res = PDBResult(state="API_ERROR", search_errors=[identity.message],
                        message="Unable to retrieve data because of an API/network error.")
        print(f"{res.message} (protein identity could not be resolved: {identity.message})")
        return res
    print("Searching RCSB (all result pages)...")
    search = search_rcsb(ctx.http, identity)
    ids = sorted(search.candidates)
    print(f"Candidates retrieved: {len(ids)}")

    if search.failed and not ids:
        res = PDBResult(state="API_ERROR", candidate_count=0, search_errors=search.errors,
                        message="Unable to retrieve data because of an API/network error.")
        print(res.message)
        return res
    if not ids:
        res = PDBResult(state="NO_STRUCTURES", candidate_count=0,
                        message=f"No crystal structures are available for {identity.input_name}.")
        print(res.message)
        return res

    outcomes: List[EntryOutcome] = []
    print("Validating identity / organism / method and collecting metadata, validation, ligands, activity...")
    with ThreadPoolExecutor(max_workers=opts.workers) as pool:
        futures = {pool.submit(process_pdb_entry, ctx, pid): pid for pid in ids}
        bar = tqdm(total=len(futures), desc="PDB entries")
        for fut in as_completed(futures):
            outcomes.append(fut.result())
            bar.update(1)
        bar.close()

    res = filter_pdb_structures(outcomes, len(ids))
    res.search_errors = search.errors
    for k, v in res.funnel.items():
        print(f"  {k}: {v}")
    for rec in res.records:
        for w in rec.warnings:
            log.warning("PDB %s: %s Continuing...", rec.pdb_id, w)

    incomplete = bool(res.errors) or search.failed
    if res.records:
        res.state = "PARTIAL" if incomplete else "OK"
        if incomplete:
            res.message = (f"Results may be incomplete: {len(res.errors)} entries and "
                           f"{len(search.failed_routes)} search routes failed because of API/network errors.")
    elif incomplete:
        res.state = "API_ERROR"
        res.message = "Unable to retrieve data because of an API/network error."
    else:
        res.state = "NO_STRUCTURES"
        res.message = f"No crystal structures are available for {identity.input_name}."
    print(f"\nTotal qualifying X-ray crystal structures found for {identity.input_name}: {len(res.records)}")
    if res.state != "OK":
        print(res.message)
    return res


def run_alphafold_branch(http: HttpClient, identity: ProteinIdentity) -> AFResult:
    _section("ALPHAFOLD SEARCH")
    print("Searching AlphaFold DB (always performed, independent of PDB results)...")
    raw, state, msg = search_alphafold(http, identity)
    if state == "API_ERROR":
        res = AFResult(state="API_ERROR",
                       message="Unable to retrieve data because of an API/network error.")
        log.warning(msg)
        print(res.message)
        return res
    if state in ("NO_STRUCTURES", "NO_IDENTITY"):
        res = AFResult(state=state, message="No AlphaFold structure is available.")
        print(f"No AlphaFold structure is available for {identity.input_name}.")
        if state == "NO_IDENTITY":
            print(f"  ({msg})")
        return res
    print(f"AlphaFold entries found: {len(raw)}")
    res = retrieve_alphafold_data(raw, identity)
    for r in res.rejected:
        log.warning("AlphaFold entry rejected: %s", r)
    if not res.records:
        res.message = "No AlphaFold structure is available."
        print(f"No AlphaFold structure is available for {identity.input_name}.")
    else:
        print(f"Qualifying AlphaFold entries: {len(res.records)}")
    return res


def run_pipeline(protein_name: str, opts: Optional[RuntimeOptions] = None) -> Dict[str, Any]:
    opts = opts or RuntimeOptions()
    setup_logging(opts.verbose)
    protein_name = (protein_name or "").strip()
    if not protein_name and not opts.uniprot_override:
        raise ValueError("A protein name is required.")

    cache = DiskCache(CACHE_DIR, enabled=opts.use_cache, refresh=opts.refresh_cache)
    http = HttpClient(cache, workers=opts.workers)
    started = datetime.now()

    _banner("PROTEIN STRUCTURE AUTOMATION")
    print(f"\nInput protein:\n{protein_name}\n\nResolving protein identity...")
    identity = resolve_protein_identity(http, protein_name, opts.uniprot_override)
    if identity.status == "RESOLVED":
        print(f"Official name: {identity.official_name}")
        print(f"Gene symbol: {identity.gene_symbol or MISSING_VALUE}")
        print(f"UniProt accession: {identity.accession} ({identity.reviewed_label or MISSING_VALUE})")
        syn = [a for a in identity.aliases if a not in (identity.official_name, identity.gene_symbol)]
        print(f"Synonyms identified ({len(syn)}): {'; '.join(syn[:12])}{' ...' if len(syn) > 12 else ''}")
        print(f"Match basis: {identity.match_basis} (score {identity.match_score})")
        if identity.message:
            print(f"WARNING: {identity.message}")
            if identity.alternatives:
                print("  Other candidates: " + " | ".join(identity.alternatives))
    else:
        print(f"WARNING: identity not resolved - {identity.message}")

    ctx = PipelineContext(http=http, chem_store=ChemCompStore(http), identity=identity,
                          skip_average_b=opts.skip_average_b, skip_cif_precision=opts.skip_cif_precision)

    # ---- two independent branches: neither depends on the other ----
    pdb_res = run_pdb_branch(ctx, opts)
    af_res = run_alphafold_branch(http, identity)

    issues = validate_records(pdb_res.records, identity)
    for msg in issues:
        log.warning("Consistency check: %s", msg)

    _section("EXCEL OUTPUT")
    opts.outdir.mkdir(parents=True, exist_ok=True)
    filename = f"{sanitize_filename(identity.input_name or protein_name)}_PDB_AlphaFold_Data.xlsx"
    print("Creating workbook...")
    out_path = create_excel_workbook(opts.outdir / filename, identity, pdb_res, af_res, issues, started)
    print("PDB_Data created\nAlphaFold_Data created")
    print(f"\nOutput file:\n{out_path}")
    _banner(f"PROCESS COMPLETED ({(datetime.now() - started).total_seconds():.1f}s)")
    return {"identity": identity, "pdb": pdb_res, "alphafold": af_res, "issues": issues,
            "output_path": out_path}

# %% [CLI]
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Collect PDB X-ray structures + AlphaFold entries for a protein.")
    ap.add_argument("protein", nargs="*", help="protein name / gene symbol / synonym / UniProt accession")
    ap.add_argument("--outdir", default=".", help="output directory (default: current directory)")
    ap.add_argument("--uniprot", help="force a UniProt accession (skips name resolution)")
    ap.add_argument("--workers", type=int, default=MAX_WORKERS, help=f"concurrent requests (default {MAX_WORKERS})")
    ap.add_argument("--no-cache", action="store_true", help="disable the disk cache")
    ap.add_argument("--refresh", action="store_true", help="ignore cached data and refetch everything")
    ap.add_argument("--skip-average-b", action="store_true",
                    help="skip downloading validation PDFs (Average B will be N/A) - much faster")
    ap.add_argument("--skip-cif-precision", action="store_true",
                    help="skip mmCIF headers (trailing zeros of resolution/R-values may be lost)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    name = " ".join(args.protein).strip() or (args.uniprot or "")
    if not name:
        name = input("Enter protein name: ").strip()
    opts = RuntimeOptions(outdir=Path(args.outdir), workers=max(1, args.workers),
                          use_cache=not args.no_cache, refresh_cache=args.refresh,
                          skip_average_b=args.skip_average_b, skip_cif_precision=args.skip_cif_precision,
                          uniprot_override=args.uniprot, verbose=args.verbose)
    try:
        run_pipeline(name, opts)
    except KeyboardInterrupt:
        print("\nInterrupted.")
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
