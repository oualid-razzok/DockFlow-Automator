"""interactions.json schema validation (audit item 21).

Validates the analyzer's real output against
``docs/schemas/interactions.schema.json`` using a small in-tree JSON-Schema
subset validator (type/required/properties/additionalProperties/items/
$ref/enum/minimum/anyOf/minItems) - deliberately dependency-free so no new
optional dependency is introduced for testing. CI runs this test as the
schema conformance gate.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from dockflow_core.analyzer import analyze_docking_result
from dockflow_core.models import DockingResult, PoseRecord
from dockflow_core.preparator import ReceptorPreparator, ReceptorPrepOptions

SCHEMA_PATH = Path(__file__).resolve().parents[1] / "docs" / "schemas" / "interactions.schema.json"


# ---------------------------------------------------------------------------
# Minimal JSON-Schema subset validator
# ---------------------------------------------------------------------------
def _resolve_ref(root: dict, ref: str) -> dict:
    assert ref.startswith("#/"), f"only local refs supported: {ref}"
    node: Any = root
    for part in ref[2:].split("/"):
        node = node[part]
    return node


def validate_schema(
    instance: Any, schema: dict, root: dict | None = None, path: str = "$"
) -> list[str]:
    """Return the list of violations ([] when the instance conforms)."""
    root = root if root is not None else schema
    if "$ref" in schema:
        return validate_schema(instance, _resolve_ref(root, schema["$ref"]), root, path)
    errors: list[str] = []
    if "anyOf" in schema:
        if not any(not validate_schema(instance, sub, root, path) for sub in schema["anyOf"]):
            errors.append(f"{path}: matched none of anyOf")
        return errors
    expected = schema.get("type")
    if expected == "object":
        if not isinstance(instance, dict):
            return [f"{path}: expected object, got {type(instance).__name__}"]
        for key in schema.get("required", []):
            if key not in instance:
                errors.append(f"{path}: missing required property {key!r}")
        properties = schema.get("properties", {})
        additional = schema.get("additionalProperties", True)
        for key, value in instance.items():
            if key in properties:
                errors += validate_schema(value, properties[key], root, f"{path}.{key}")
            elif isinstance(additional, dict):
                errors += validate_schema(value, additional, root, f"{path}.{key}")
            elif additional is False:
                errors.append(f"{path}: unexpected property {key!r}")
        if "minItems" in schema and len(instance) < schema["minItems"]:
            errors.append(f"{path}: fewer than minItems entries")
    elif expected == "array":
        if not isinstance(instance, list):
            return [f"{path}: expected array, got {type(instance).__name__}"]
        items = schema.get("items")
        if items:
            for index, item in enumerate(instance):
                errors += validate_schema(item, items, root, f"{path}[{index}]")
        if len(instance) < schema.get("minItems", 0):
            errors.append(f"{path}: fewer than minItems entries")
    elif expected == "string":
        if not isinstance(instance, str):
            errors.append(f"{path}: expected string, got {type(instance).__name__}")
    elif expected == "integer":
        if isinstance(instance, bool) or not isinstance(instance, int):
            errors.append(f"{path}: expected integer, got {instance!r}")
    elif expected == "number":
        if isinstance(instance, bool) or not isinstance(instance, (int, float)):
            errors.append(f"{path}: expected number, got {instance!r}")
    elif expected == "boolean":
        if not isinstance(instance, bool):
            errors.append(f"{path}: expected boolean, got {instance!r}")
    elif expected == "null":
        if instance is not None:
            errors.append(f"{path}: expected null, got {instance!r}")
    if "enum" in schema and instance not in schema["enum"]:
        errors.append(f"{path}: {instance!r} not in enum {schema['enum']}")
    if "minimum" in schema and isinstance(instance, (int, float)) \
            and not isinstance(instance, bool):
        if instance < schema["minimum"]:
            errors.append(f"{path}: {instance} < minimum {schema['minimum']}")
    return errors


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
@pytest.fixture
def schema() -> dict:
    return json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))


@pytest.fixture
def real_payload(receptor_pdb_path: Path, docked_pdbqt_path: Path,
                 tmp_path: Path) -> dict:
    prep = ReceptorPreparator(
        ReceptorPrepOptions(engine="none", charge_model="zero", keep_resnames=["BEN"])
    ).prepare(receptor_pdb_path, tmp_path)
    result = DockingResult(
        ligand_name="lig",
        poses=[PoseRecord(model=1, affinity=-9.423)],
        out_path=docked_pdbqt_path,
    )
    analyses = analyze_docking_result(result, prep.pdbqt_path, top_poses=3)
    return {"lig": [analysis.to_dict() for analysis in analyses]}


def test_schema_file_exists_and_parses():
    assert SCHEMA_PATH.is_file()
    data = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    assert data["title"].startswith("DockFlow")
    assert "definitions" in data


def test_analyzer_output_conforms_to_schema(schema, real_payload):
    errors = validate_schema(real_payload, schema)
    assert errors == [], f"schema violations: {errors}"


def test_schema_rejects_missing_required_field(schema, real_payload):
    payload = json.loads(json.dumps(real_payload))
    del payload["lig"][0]["affinity"]
    errors = validate_schema(payload, schema)
    assert any("affinity" in error for error in errors)


def test_schema_rejects_wrong_contact_kind(schema, real_payload):
    payload = json.loads(json.dumps(real_payload))
    if payload["lig"][0]["contacts"]:
        payload["lig"][0]["contacts"][0]["kind"] = "covalent-bond"
    else:  # contacts list empty in minimal fixture - force one
        payload["lig"][0]["contacts"] = [{"kind": "ionic-bond"}]
    errors = validate_schema(payload, schema)
    assert any("enum" in error for error in errors)


def test_schema_rejects_unqualified_legacy_keys(schema, real_payload):
    """The v0.1 key names (hbonds, ligand_efficiency) must fail validation."""
    payload = json.loads(json.dumps(real_payload))
    payload["lig"][0]["hbonds"] = 2
    errors = validate_schema(payload, schema)
    assert any("hbonds" in error for error in errors)


def test_validator_catches_type_mismatches(schema):
    payload = {"lig": [{"pose_index": "1", "affinity": -9.0,
                        "num_contacts": 0, "contacts": [], "residues": []}]}
    errors = validate_schema(payload, schema)
    assert any("pose_index" in error for error in errors)
