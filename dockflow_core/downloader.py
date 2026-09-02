"""Target and ligand downloading: RCSB PDB, UniProt/AlphaFold, PubChem, ZINC.

Design notes
------------
* Every network call goes through :class:`HTTPClient`, which adds retries,
  timeouts, a descriptive user agent and an on-disk cache.
* Failures raise :class:`DownloadError` with actionable messages (never bare
  ``requests`` exceptions).
* Structure ids are validated before requests are issued.

Endpoints used
--------------
RCSB download  https://files.rcsb.org/download/{id}.pdb
RCSB ligands   https://files.rcsb.org/ligands/view/{code}_ideal.sdf
RCSB metadata  https://data.rcsb.org/rest/v1/core/entry/{id}
PDBe mappings  https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{acc}
AlphaFold      https://alphafold.ebi.ac.uk/api/prediction/{acc}
PubChem PUG    https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/...
ZINC22         https://zinc22.docking.org/substances/{zinc_id}.json
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import quote

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .models import LigandRecord, ProteinRecord
from .pdbio import is_valid_pdb_id, parse_pdb
from .utils import DockFlowError, ensure_dir, get_logger

logger = get_logger("downloader")

__all__ = [
    "DownloadError",
    "HTTPClient",
    "PDBDownloader",
    "LigandDownloader",
    "TargetResolver",
    "download_pdb",
    "download_ligand_from_pdb",
    "resolve_target",
]

USER_AGENT = "DockFlow-Automator/0.1.0 (github.com/dockflow/DockFlow-Automator)"

_UNIPROT_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}){2}[0-9]")
# Syntax characters that essentially never occur in compound names but are
# ubiquitous in SMILES (=, #, brackets, branches, stereo marks, ring digits).
_SMILES_HINT = re.compile(r"[=#\[\]@\\/()\.]|Cl|Br|[cnpso][0-9]")


class DownloadError(DockFlowError):
    """A network resource could not be retrieved."""


# ---------------------------------------------------------------------------
# HTTP client
# ---------------------------------------------------------------------------
class HTTPClient:
    """A small requests session wrapper with retries and a file cache."""

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        timeout: float = 60.0,
        retries: int = 3,
        offline: bool = False,
        session: requests.Session | None = None,
    ) -> None:
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            ensure_dir(self.cache_dir)
        self.timeout = timeout
        self.offline = offline
        self.session = session or requests.Session()
        retry = Retry(
            total=retries,
            backoff_factor=0.8,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=("GET", "POST", "HEAD"),
        )
        adapter = HTTPAdapter(max_retries=retry, pool_maxsize=4)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
        self.session.headers.update({"User-Agent": USER_AGENT, "Accept": "*/*"})

    # -- core ---------------------------------------------------------------
    def get(
        self,
        url: str,
        params: dict[str, Any] | None = None,
        stream: bool = False,
    ) -> requests.Response:
        if self.offline:
            raise DownloadError(
                f"offline mode is enabled; refusing to fetch {url}. "
                "Re-run without --offline or provide local files."
            )
        logger.debug("GET %s params=%s", url, params)
        try:
            response = self.session.get(url, params=params, timeout=self.timeout, stream=stream)
        except requests.RequestException as exc:
            raise DownloadError(f"network error while fetching {url}: {exc}") from exc
        if response.status_code >= 400:
            response.close()
            raise DownloadError(
                f"HTTP {response.status_code} for {url}"
                + (f" ({response.reason})" if hasattr(response, "reason") else "")
            )
        return response

    def get_text(self, url: str, params: dict[str, Any] | None = None) -> str:
        response = self.get(url, params=params)
        text = response.text
        response.close()
        return text

    def get_json(self, url: str, params: dict[str, Any] | None = None) -> Any:
        text = self.get_text(url, params=params)
        try:
            import json

            return json.loads(text)
        except ValueError as exc:
            raise DownloadError(f"invalid JSON from {url}: {exc}") from exc

    def download_file(self, url: str, destination: str | Path, force: bool = False) -> Path:
        """Download ``url`` to ``destination`` atomically, with caching."""
        dest = Path(destination)
        if dest.is_file() and dest.stat().st_size > 0 and not force:
            logger.debug("cache hit: %s", dest)
            return dest
        dest.parent.mkdir(parents=True, exist_ok=True)
        if self.offline:
            raise DownloadError(f"offline mode: cannot download {url}")
        logger.info("downloading %s -> %s", url, dest)
        response = self.get(url, stream=True)
        tmp = dest.with_suffix(dest.suffix + ".part")
        try:
            with open(tmp, "wb") as fh:
                for chunk in response.iter_content(chunk_size=1 << 16):
                    if chunk:
                        fh.write(chunk)
            tmp.replace(dest)
        finally:
            response.close()
            if tmp.exists():
                tmp.unlink(missing_ok=True)
        if dest.stat().st_size == 0:
            raise DownloadError(f"empty response for {url}")
        return dest


def _json_path(data: Any, *keys: str, default: Any = None) -> Any:
    """Safe nested lookup in parsed JSON."""
    node = data
    for key in keys:
        if not isinstance(node, dict):
            return default
        node = node.get(key)
        if node is None:
            return default
    return node


# ---------------------------------------------------------------------------
# PDB / UniProt downloader
# ---------------------------------------------------------------------------
class PDBDownloader:
    """Fetches experimentally solved structures from RCSB and AlphaFold."""

    RCSB_FILES = "https://files.rcsb.org/download/{pdb_id}.pdb"
    RCSB_CIF = "https://files.rcsb.org/download/{pdb_id}.cif"
    RCSB_LIGAND_SDF = "https://files.rcsb.org/ligands/view/{code}_ideal.sdf"
    RCSB_LIGAND_MODEL_SDF = "https://files.rcsb.org/ligands/view/{code}_model.sdf"
    RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    RCSB_LIGAND_PDB = "https://files.rcsb.org/ligands/view/{code}.pdb"
    PDBE_UNIPROT = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{acc}"
    ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

    def __init__(self, client: HTTPClient | None = None,
                 cache_dir: str | Path | None = None) -> None:
        self.client = client or HTTPClient(cache_dir=cache_dir)

    # -- structures ---------------------------------------------------------
    @staticmethod
    def normalize_pdb_id(pdb_id: str) -> str:
        clean = pdb_id.strip().upper().removeprefix("PDB:")
        if not is_valid_pdb_id(clean):
            raise DownloadError(
                f"invalid PDB id {pdb_id!r}: expected 4 characters, first a digit "
                "(e.g. 1HVR)"
            )
        return clean

    def fetch_structure(
        self,
        pdb_id: str,
        output_dir: str | Path,
        fmt: str = "pdb",
        force: bool = False,
    ) -> ProteinRecord:
        """Download a PDB entry (``pdb`` or ``cif``) and return a record."""
        clean = self.normalize_pdb_id(pdb_id)
        output_dir = ensure_dir(output_dir)
        if fmt not in ("pdb", "cif"):
            raise DownloadError(f"unsupported format {fmt!r} (use 'pdb' or 'cif')")
        url = (self.RCSB_FILES if fmt == "pdb" else self.RCSB_CIF).format(pdb_id=clean)
        path = output_dir / f"{clean.lower()}.{fmt}"
        self.client.download_file(url, path, force=force)
        record = ProteinRecord(identifier=clean, source="pdb", path=path)
        try:
            self._annotate_record(record)
        except DownloadError as exc:
            logger.warning("metadata unavailable for %s: %s", clean, exc)
        try:
            record.ligand_codes = self.entry_ligand_codes(clean, structure_path=path)
        except DownloadError as exc:
            logger.debug("ligand listing from structure failed: %s", exc)
        return record

    def _annotate_record(self, record: ProteinRecord) -> None:
        url = self.RCSB_ENTRY.format(pdb_id=record.identifier)
        data = self.client.get_json(url)
        record.title = _json_path(data, "struct", "title", default="") or ""
        resolutions = _json_path(data, "rcsb_entry_info", "resolution_combined", default=[]) or []
        if resolutions:
            record.resolution = float(resolutions[0])
        uniprots = (
            _json_path(data, "rcsb_entry_container_identifiers", "uniprot_ids", default=[]) or []
        )
        record.uniprot_ids = [str(u) for u in uniprots][:10]

    def entry_ligand_codes(
        self, pdb_id: str, structure_path: str | Path | None = None
    ) -> list[str]:
        """Ligand-like HET residue codes of an entry (parsed locally when possible)."""
        if structure_path and Path(structure_path).is_file():
            atoms = parse_pdb(structure_path)
            codes: list[str] = []
            seen: set[str] = set()
            for atom in atoms:
                name = atom.resname.strip()
                if atom.is_polymer or atom.is_water or not name.isalnum():
                    continue
                if name not in seen:
                    seen.add(name)
                    codes.append(name)
            return codes
        # Fall back to parsing the downloaded structure without saving it.
        clean = self.normalize_pdb_id(pdb_id)
        text = self.client.get_text(self.RCSB_FILES.format(pdb_id=clean))
        atoms = parse_pdb(text)
        codes = []
        seen = set()
        for atom in atoms:
            name = atom.resname.strip()
            if atom.is_polymer or atom.is_water or not name.isalnum():
                continue
            if name not in seen:
                seen.add(name)
                codes.append(name)
        return codes

    # -- ligands ------------------------------------------------------------
    def fetch_ligand(
        self,
        resname: str,
        output_dir: str | Path,
        prefer: str = "ideal",
        force: bool = False,
    ) -> Path:
        """Download a chemical component from RCSB as an SDF file.

        ``ideal`` returns the idealized (0 K minimized) coordinates,
        ``model`` the coordinates as found in the first PDB entry.
        """
        code = resname.strip().upper()
        if not re.fullmatch(r"[A-Z0-9]{1,3}", code):
            raise DownloadError(f"invalid ligand code {resname!r}")
        output_dir = ensure_dir(output_dir)
        if prefer == "model":
            urls = [self.RCSB_LIGAND_MODEL_SDF, self.RCSB_LIGAND_SDF, self.RCSB_LIGAND_PDB]
        else:
            urls = [self.RCSB_LIGAND_SDF, self.RCSB_LIGAND_MODEL_SDF, self.RCSB_LIGAND_PDB]
        last_error: DownloadError | None = None
        for template in urls:
            url = template.format(code=code)
            suffix = ".sdf" if ".sdf" in template else ".pdb"
            path = output_dir / f"{code.lower()}{suffix}"
            try:
                return self.client.download_file(url, path, force=force)
            except DownloadError as exc:
                last_error = exc
                logger.debug("ligand endpoint %s failed: %s", url, exc)
        raise DownloadError(
            f"ligand {code!r} not found on RCSB (tried SDF/PDB endpoints): {last_error}"
        )

    # -- UniProt ------------------------------------------------------------
    @staticmethod
    def normalize_uniprot(acc: str) -> str:
        clean = acc.strip().upper().removeprefix("UNIPROT:").split()[0] if acc.strip() else ""
        if not _UNIPROT_RE.match(clean):
            raise DownloadError(f"invalid UniProt accession {acc!r} (e.g. P29978)")
        return clean

    def pdbs_for_uniprot(self, acc: str) -> list[str]:
        """PDB entries mapped to a UniProt accession (PDBe API)."""
        clean = self.normalize_uniprot(acc)
        url = self.PDBE_UNIPROT.format(acc=clean)
        data = self.client.get_json(url)
        # {"P29978": {"PDB": [{"PDB_id": "1a4x", ...}, ...]}, ...}
        entries: list[str] = []
        if isinstance(data, dict):
            for value in data.values():
                mappings = value.get("PDB", {}) if isinstance(value, dict) else {}
                if isinstance(mappings, dict):
                    # Actual shape: value["PDB"] is a dict of pdb_id -> list
                    for pdb_id in mappings:
                        entries.append(str(pdb_id).upper())
                elif isinstance(mappings, list):
                    for item in mappings:
                        pdb_id = item.get("PDB_id") or item.get("pdb_id")
                        if pdb_id:
                            entries.append(str(pdb_id).upper())
        if not entries and isinstance(data, dict):
            # Some responses nest differently; do a shallow sweep.
            for key in data:
                if re.fullmatch(r"[0-9][a-z0-9]{3}", key):
                    entries.append(key.upper())
        return sorted(set(entries))

    def fetch_alphafold_model(
        self, acc: str, output_dir: str | Path, force: bool = False
    ) -> ProteinRecord:
        """Download the AlphaFold predicted model for a UniProt accession."""
        clean = self.normalize_uniprot(acc)
        output_dir = ensure_dir(output_dir)
        url = self.ALPHAFOLD_API.format(acc=clean)
        data = self.client.get_json(url)
        if not isinstance(data, list) or not data:
            raise DownloadError(f"no AlphaFold prediction for {clean}")
        entry = data[0]
        pdb_url = entry.get("pdbUrl") or entry.get("pdb_url")
        if not pdb_url:
            raise DownloadError(f"AlphaFold response for {clean} lacks a pdbUrl")
        path = output_dir / f"af_{clean.lower()}.pdb"
        self.client.download_file(pdb_url, path, force=force)
        record = ProteinRecord(
            identifier=clean,
            source="uniprot-alphafold",
            path=path,
            title=(entry.get("uniprotDescription") or f"AlphaFold prediction for {clean}"),
            resolution=None,
            uniprot_ids=[clean],
        )
        return record


# ---------------------------------------------------------------------------
# Small-molecule ligand downloader (PubChem, ZINC)
# ---------------------------------------------------------------------------
class LigandDownloader:
    """Fetches small molecules from PubChem and ZINC22."""

    PUBCHEM_PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    ZINC_SUBSTANCE = "https://zinc22.docking.org/substances/{zinc_id}"
    ZINC_SUBSTANCE_JSON = "https://zinc22.docking.org/substances/{zinc_id}.json"

    def __init__(self, client: HTTPClient | None = None,
                 cache_dir: str | Path | None = None) -> None:
        self.client = client or HTTPClient(cache_dir=cache_dir)

    # -- PubChem ------------------------------------------------------------
    @staticmethod
    def _looks_like_smiles(text: str) -> bool:
        """Heuristic: real SMILES almost always contain syntax characters.

        Plain names ("aspirin") and CIDs ("2244") do not.
        """
        if not text or " " in text:
            return False
        return bool(_SMILES_HINT.search(text))

    def pubchem_cid_from_name(self, name: str) -> int:
        text = self.client.get_text(f"{self.PUBCHEM_PUG}/compound/name/{quote(name)}/cids/TXT")
        for token in text.split():
            if token.strip().isdigit():
                return int(token.strip())
        raise DownloadError(f"PubChem found no CID for name {name!r}")

    def pubchem_cid_from_smiles(self, smiles: str) -> int:
        """Resolve a SMILES string to a PubChem CID (identity service)."""
        try:
            text = self.client.get_text(
                f"{self.PUBCHEM_PUG}/compound/identity/smiles/cids/TXT",
                params={"smiles": smiles},
            )
        except DownloadError:
            raise DownloadError(f"PubChem found no CID for SMILES {smiles!r}") from None
        for token in text.split():
            if token.strip().isdigit():
                return int(token.strip())
        raise DownloadError(f"PubChem found no CID for SMILES {smiles!r}")

    def fetch_pubchem(
        self,
        identifier: str,
        output_dir: str | Path,
        prefer_3d: bool = True,
        force: bool = False,
    ) -> LigandRecord:
        """Download a compound from PubChem by name, CID or SMILES.

        3D conformers are requested when available; otherwise the 2D record
        is returned and 3D coordinates are embedded later by the preparator.
        """
        identifier = identifier.strip()
        if not identifier:
            raise DownloadError("empty PubChem identifier")
        if identifier.isdigit():
            cid = int(identifier)
        elif self._looks_like_smiles(identifier):
            cid = self.pubchem_cid_from_smiles(identifier)
        else:
            try:
                cid = self.pubchem_cid_from_name(identifier)
            except DownloadError:
                # a bare formula like "CCO" can still be a SMILES
                cid = self.pubchem_cid_from_smiles(identifier)
        output_dir = ensure_dir(output_dir)
        attempts: list[str] = []
        if prefer_3d:
            attempts.append(f"{self.PUBCHEM_PUG}/compound/cid/{cid}/SDF?record_type=3d")
        attempts.append(f"{self.PUBCHEM_PUG}/compound/cid/{cid}/SDF?record_type=2d")
        path = output_dir / f"pubchem_{cid}.sdf"
        last_error: DownloadError | None = None
        for url in attempts:
            try:
                self.client.download_file(url, path, force=force)
                break
            except DownloadError as exc:
                last_error = exc
                logger.debug("PubChem endpoint failed (%s): %s", url, exc)
        else:
            raise DownloadError(f"PubChem SDF unavailable for CID {cid}: {last_error}")
        properties = self.pubchem_properties(cid, fields=("CanonicalSMILES", "MolecularWeight"))
        record = LigandRecord(
            identifier=f"pubchem_{cid}",
            source="pubchem",
            value=identifier,
            path=path,
            status="downloaded",
        )
        record.num_rotatable_bonds = None
        if properties:
            smiles = properties.get("CanonicalSMILES")
            if smiles:
                record.value = str(smiles)
        return record

    def pubchem_properties(self, cid: int, fields: Sequence[str]) -> dict[str, Any]:
        url = f"{self.PUBCHEM_PUG}/compound/cid/{cid}/property/{','.join(fields)}/JSON"
        try:
            data = self.client.get_json(url)
        except DownloadError as exc:
            logger.debug("PubChem properties failed: %s", exc)
            return {}
        props = _json_path(data, "PropertyTable", "Properties", default=[])
        if isinstance(props, list) and props and isinstance(props[0], dict):
            return props[0]
        return {}

    # -- ZINC ----------------------------------------------------------------
    @staticmethod
    def normalize_zinc_id(zinc_id: str) -> str:
        clean = zinc_id.strip().upper().removeprefix("ZINC:")
        if not re.fullmatch(r"ZINC\d{6,15}", clean):
            raise DownloadError(f"invalid ZINC id {zinc_id!r} (e.g. ZINC000000000001)")
        return clean

    def fetch_zinc(
        self, zinc_id: str, output_dir: str | Path, force: bool = False
    ) -> LigandRecord:
        """Download a substance from ZINC22 (3D SDF preferred, SMILES fallback)."""
        clean = self.normalize_zinc_id(zinc_id)
        output_dir = ensure_dir(output_dir)
        record = LigandRecord(
            identifier=clean.lower(),
            source="zinc",
            value=clean,
            status="pending",
        )
        # 1) Try the direct SDF endpoint.
        try:
            path = output_dir / f"{clean.lower()}.sdf"
            self.client.download_file(self.ZINC_SUBSTANCE.format(zinc_id=clean) + ".sdf", path,
                                      force=force)
            record.path = path
            record.status = "downloaded"
            return record
        except DownloadError as exc:
            logger.debug("ZINC SDF endpoint failed for %s: %s", clean, exc)
        # 2) Resolve the SMILES through the JSON API; the preparator will
        #    embed 3D coordinates from it.
        smiles = self.zinc_smiles(clean)
        record.value = smiles
        record.status = "downloaded"
        path = output_dir / f"{clean.lower()}.smi"
        path.write_text(f"{smiles} {clean}\n", encoding="utf-8")
        record.path = path
        return record

    def zinc_smiles(self, zinc_id: str) -> str:
        clean = self.normalize_zinc_id(zinc_id)
        url = self.ZINC_SUBSTANCE_JSON.format(zinc_id=clean)
        data = self.client.get_json(url)
        if isinstance(data, dict):
            substances = data.get("substances") or [data]
        elif isinstance(data, list):
            substances = data
        else:
            substances = []
        for substance in substances:
            if not isinstance(substance, dict):
                continue
            for key in ("smiles", "smiles-string", "canonical_smiles"):
                value = substance.get(key)
                if value:
                    return str(value)
        raise DownloadError(f"ZINC returned no SMILES for {clean}")


# ---------------------------------------------------------------------------
# Target resolution (used by the pipeline and the GUI)
# ---------------------------------------------------------------------------
@dataclass
class TargetSpec:
    """Where the macromolecular target comes from."""

    pdb_id: str | None = None
    uniprot: str | None = None
    file: str | Path | None = None
    prefer_alphafold: bool = False

    @classmethod
    def parse(cls, value: str | Path | dict[str, Any]) -> TargetSpec:
        if isinstance(value, dict):
            return cls(
                pdb_id=value.get("pdb_id"),
                uniprot=value.get("uniprot"),
                file=value.get("file"),
                prefer_alphafold=bool(value.get("prefer_alphafold", False)),
            )
        text = str(value).strip()
        if not text:
            raise DownloadError("empty target specification")
        if text.lower().startswith("pdb:"):
            return cls(pdb_id=text[4:])
        if text.lower().startswith("uniprot:"):
            return cls(uniprot=text[8:])
        if text.lower().startswith("file:"):
            return cls(file=text[5:])
        candidate = Path(text)
        if candidate.exists() and candidate.is_file():
            return cls(file=candidate)
        if is_valid_pdb_id(text.upper()):
            return cls(pdb_id=text)
        if _UNIPROT_RE.match(text.upper()):
            return cls(uniprot=text)
        raise DownloadError(
            f"cannot interpret target {value!r}: use a PDB id (1HVR), a UniProt "
            f"accession (P29978), a file path, or an explicit pdb:/uniprot:/file: prefix"
        )


class TargetResolver:
    """Resolves a :class:`TargetSpec` into a local structure file."""

    def __init__(self, downloader: PDBDownloader | None = None) -> None:
        self.downloader = downloader or PDBDownloader()

    def resolve(self, spec: TargetSpec | str | Path, output_dir: str | Path) -> ProteinRecord:
        if not isinstance(spec, TargetSpec):
            spec = TargetSpec.parse(spec)
        if spec.file:
            path = Path(spec.file)
            if not path.is_file():
                raise DownloadError(f"target file not found: {path}")
            identifier = path.stem
            codes: list[str] = []
            if path.suffix.lower() in {".pdb", ".ent"}:
                atoms = parse_pdb(path)
                seen: set[str] = set()
                for atom in atoms:
                    name = atom.resname.strip()
                    if atom.is_polymer or atom.is_water or not name.isalnum():
                        continue
                    if name not in seen:
                        seen.add(name)
                        codes.append(name)
            return ProteinRecord(
                identifier=identifier, source="file", path=path, ligand_codes=codes
            )
        if spec.pdb_id:
            return self.downloader.fetch_structure(spec.pdb_id, output_dir)
        if spec.uniprot:
            # Try experimental structures first, fall back to AlphaFold.
            try:
                pdb_ids = self.downloader.pdbs_for_uniprot(spec.uniprot)
            except DownloadError as exc:
                logger.warning("UniProt mapping lookup failed: %s", exc)
                pdb_ids = []
            if pdb_ids and not spec.prefer_alphafold:
                logger.info("UniProt %s mapped to %d PDB entries; using %s",
                            spec.uniprot, len(pdb_ids), pdb_ids[0])
                return self.downloader.fetch_structure(pdb_ids[0], output_dir)
            return self.downloader.fetch_alphafold_model(spec.uniprot, output_dir)
        raise DownloadError("target specification is empty")


# ---------------------------------------------------------------------------
# Convenience top-level functions
# ---------------------------------------------------------------------------
def download_pdb(pdb_id: str, output_dir: str | Path, client: HTTPClient | None = None) -> Path:
    """Download a PDB structure and return its path."""
    record = PDBDownloader(client).fetch_structure(pdb_id, output_dir)
    assert record.path is not None
    return record.path


def download_ligand_from_pdb(
    resname: str, output_dir: str | Path, client: HTTPClient | None = None
) -> Path:
    """Download an RCSB chemical component as SDF and return its path."""
    return PDBDownloader(client).fetch_ligand(resname, output_dir)


def resolve_target(spec: str | Path | TargetSpec, output_dir: str | Path) -> ProteinRecord:
    """Resolve any supported target specification to a local structure."""
    return TargetResolver().resolve(spec, output_dir)

