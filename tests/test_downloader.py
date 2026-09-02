"""Downloader tests with fully mocked HTTP (no network access)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import requests

from dockflow_core.downloader import (
    DownloadError,
    HTTPClient,
    LigandDownloader,
    PDBDownloader,
    TargetResolver,
    TargetSpec,
)


class FakeResponse:
    def __init__(self, status_code=200, text="", content=b""):
        self.status_code = status_code
        self.text = text
        self.content = content if content else text.encode("utf-8")
        self.reason = "OK" if status_code < 400 else "Not Found"

    def iter_content(self, chunk_size=1 << 16):
        yield self.content

    def close(self):
        pass


def patch_client(monkeypatch, handler):
    """Replace HTTPClient.get with ``handler(url, params) -> FakeResponse``."""

    def fake_get(self, url, params=None, stream=False):
        return handler(url, params)

    monkeypatch.setattr(HTTPClient, "get", fake_get)


def test_normalize_pdb_id():
    assert PDBDownloader.normalize_pdb_id("1hvr") == "1HVR"
    with pytest.raises(DownloadError):
        PDBDownloader.normalize_pdb_id("XXXX")
    with pytest.raises(DownloadError):
        PDBDownloader.normalize_pdb_id("12")


def test_normalize_uniprot():
    assert PDBDownloader.normalize_uniprot("p29978") == "P29978"
    with pytest.raises(DownloadError):
        PDBDownloader.normalize_uniprot("NOTREAL")


def test_normalize_zinc_id():
    assert LigandDownloader.normalize_zinc_id("zinc000000000001") == "ZINC000000000001"
    with pytest.raises(DownloadError):
        LigandDownloader.normalize_zinc_id("12345")


def test_fetch_structure(monkeypatch, tmp_path: Path):
    calls = []

    def handler(url, params):
        calls.append(url)
        if "download" in url:
            pdb_id = url.rsplit("/", 1)[-1].split(".")[0]
            assert pdb_id == "1HVR"
            return FakeResponse(text="ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n")
        if "data.rcsb.org" in url:
            return FakeResponse(text=json.dumps({
                "struct": {"title": "HIV-1 PROTEASE COMPLEX"},
                "rcsb_entry_info": {"resolution_combined": [2.0]},
                "rcsb_entry_container_identifiers": {"uniprot_ids": ["P03366"]},
            }))
        raise AssertionError(f"unexpected url {url}")

    patch_client(monkeypatch, handler)
    record = PDBDownloader().fetch_structure("1hvr", tmp_path)
    assert record.identifier == "1HVR"
    assert record.title == "HIV-1 PROTEASE COMPLEX"
    assert record.resolution == 2.0
    assert record.uniprot_ids == ["P03366"]
    assert record.ligand_codes == []  # only one atom, no HETATM
    assert record.path is not None and record.path.is_file()

    # second call hits the file cache (no new download request)
    downloads_before = sum(1 for u in calls if "files.rcsb.org/download" in u)
    PDBDownloader().fetch_structure("1HVR", tmp_path)
    downloads_after = sum(1 for u in calls if "files.rcsb.org/download" in u)
    assert downloads_before == downloads_after == 1


def test_fetch_structure_http_error(monkeypatch, tmp_path: Path):
    def handler(url, params):
        return FakeResponse(status_code=404)

    patch_client(monkeypatch, handler)
    downloader = PDBDownloader()
    with pytest.raises(DownloadError) as excinfo:
        downloader.fetch_structure("1ZZZ", tmp_path)
    # the request failed and the error mentions the offending structure
    assert "1ZZZ" in str(excinfo.value)


def test_http_client_get_rejects_error_status(monkeypatch):
    client = HTTPClient()
    # patch the underlying session so the real status-check logic runs
    monkeypatch.setattr(
        client.session, "get",
        lambda *a, **k: FakeResponse(status_code=404),
    )
    with pytest.raises(DownloadError, match="404"):
        client.get("https://files.rcsb.org/download/1ZZZ.pdb")


def test_fetch_structure_invalid_id(monkeypatch, tmp_path: Path):
    patch_client(monkeypatch, lambda u, p: FakeResponse())
    with pytest.raises(DownloadError, match="invalid PDB id"):
        PDBDownloader().fetch_structure("ZZZZ", tmp_path)


def test_offline_mode(tmp_path: Path):
    client = HTTPClient(offline=True)
    with pytest.raises(DownloadError, match="offline"):
        client.get_text("https://example.org")
    with pytest.raises(DownloadError, match="offline"):
        client.download_file("https://example.org/x.pdb", tmp_path / "x.pdb")


def test_pdbs_for_uniprot(monkeypatch):
    payload = {"P29978": {"PDB": {"1a4x": [{"pdb_id": "1A4X"}], "2hvy": [{}]}}}

    def handler(url, params):
        assert "pdbe/api/mappings" in url
        return FakeResponse(text=json.dumps(payload))

    patch_client(monkeypatch, handler)
    ids = PDBDownloader().pdbs_for_uniprot("P29978")
    assert "1A4X" in ids and "2HVY" in ids


def test_fetch_alphafold(monkeypatch, tmp_path: Path):
    def handler(url, params):
        if "api/prediction" in url:
            return FakeResponse(text=json.dumps([
                {"pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P29978-F1.pdb",
                 "uniprotDescription": "PROTEASE"},
            ]))
        assert "files/AF" in url
        return FakeResponse(text="ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 50.00           C\nEND\n")

    patch_client(monkeypatch, handler)
    record = PDBDownloader().fetch_alphafold_model("P29978", tmp_path)
    assert record.source == "uniprot-alphafold"
    assert record.identifier == "P29978"
    assert record.path is not None and record.path.name == "af_p29978.pdb"


def test_pubchem_by_cid(monkeypatch, tmp_path: Path):
    urls = []

    def handler(url, params):
        urls.append(url)
        if "record_type" in url or url.endswith("/SDF"):
            return FakeResponse(content=b"FAKESDF")
        if "/cids/TXT" in url:
            return FakeResponse(text="2244\n")
        if "/property/" in url:
            return FakeResponse(text=json.dumps({
                "PropertyTable": {"Properties": [
                    {"CanonicalSMILES": "CC(=O)Oc1ccccc1C(=O)O"}
                ]}
            }))
        raise AssertionError(f"unexpected url {url}")

    patch_client(monkeypatch, handler)
    record = LigandDownloader().fetch_pubchem("2244", tmp_path)
    assert record.identifier == "pubchem_2244"
    assert record.path is not None and record.path.is_file()
    assert any("record_type=3d" in u for u in urls)
    assert record.value == "CC(=O)Oc1ccccc1C(=O)O"  # from properties


def test_pubchem_by_name(monkeypatch, tmp_path: Path):
    def handler(url, params):
        if "/compound/name/aspirin/cids/TXT" in url:
            return FakeResponse(text="2244")
        if "record_type" in url:
            return FakeResponse(content=b"FAKESDF")
        return FakeResponse(text="{}\n")

    patch_client(monkeypatch, handler)
    record = LigandDownloader().fetch_pubchem("aspirin", tmp_path)
    assert record.identifier == "pubchem_2244"


def test_pubchem_not_found(monkeypatch, tmp_path: Path):
    def handler(url, params):
        return FakeResponse(status_code=404)

    patch_client(monkeypatch, handler)
    with pytest.raises(DownloadError):
        LigandDownloader().fetch_pubchem("definitely-not-a-molecule-xyzzy", tmp_path)


def test_zinc_smiles_fallback(monkeypatch, tmp_path: Path):
    def handler(url, params):
        if url.endswith(".json"):
            assert "ZINC000000000001" in url
            return FakeResponse(text=json.dumps(
                {"substances": [{"smiles": "CCO"}]}))
        # SDF endpoint fails -> falls back to SMILES
        return FakeResponse(status_code=404)

    patch_client(monkeypatch, handler)
    record = LigandDownloader().fetch_zinc("ZINC000000000001", tmp_path)
    assert record.value == "CCO"
    assert record.path is not None and record.path.suffix == ".smi"


def test_target_spec_parse():
    spec = TargetSpec.parse("1HVR")
    assert spec.pdb_id == "1HVR"
    spec = TargetSpec.parse("uniprot:P29978")
    assert spec.uniprot == "P29978"
    spec = TargetSpec.parse("pdb:1HVR")
    assert spec.pdb_id == "1HVR"
    with pytest.raises(DownloadError):
        TargetSpec.parse("")
    with pytest.raises(DownloadError):
        TargetSpec.parse("!!!nonsense!!!")
    with pytest.raises(DownloadError):
        TargetSpec.parse("ZXYZ")


def test_target_resolver_local_file(receptor_pdb_path: Path):
    record = TargetResolver().resolve(str(receptor_pdb_path), receptor_pdb_path.parent)
    assert record.source == "file"
    assert record.identifier == receptor_pdb_path.stem
    assert "BEN" in record.ligand_codes and "ZN2" in record.ligand_codes


def test_download_file_atomic(monkeypatch, tmp_path: Path):
    def handler(url, params):
        return FakeResponse(content=b"DATA" * 100)

    patch_client(monkeypatch, handler)
    client = HTTPClient()
    path = client.download_file("https://example.org/file.pdb", tmp_path / "file.pdb")
    assert path.read_bytes() == b"DATA" * 100
    assert not path.with_suffix(".pdb.part").exists()


def test_real_session_headers():
    client = HTTPClient()
    assert "DockFlow" in client.session.headers["User-Agent"]
    assert isinstance(client.session, requests.Session)
