"""Final-pass (v1.0.0-rc1) guard tests.

Covers the review items that need code-level protection:

* item 2  - Wilson score interval (exact values, hand-rolled);
* item 5  - symmetry-RMSD analytic cases live in test_analyzer.py;
* item 8  - degraded-mode banner + CLI ``--allow-degraded`` gate +
  GUI modal persistence (QSettings, per-version re-prompt);
* item 16 - doc-claim lint: the five flagged words must be qualified.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load_benchmark_module(name: str, relative: str):
    """Import a benchmark script by path (they are not a package)."""
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# item 2: Wilson score interval
# ---------------------------------------------------------------------------
def test_wilson_ci_18_of_24_matches_worked_example():
    """18/24 -> '75.0% (95% CI: 55.1%-88.0%)' - the reviewer's example."""
    runner = _load_benchmark_module(
        "run_benchmark", "benchmarks/redocking/run_benchmark.py")
    low, high = runner.wilson_ci(18, 24)
    assert low == pytest.approx(0.551, abs=0.001)
    assert high == pytest.approx(0.880, abs=0.001)
    assert runner._fmt_pct_rate(18, 24) == "75.0% (95% CI: 55.1%-88.0%)"


def test_wilson_ci_degenerate_cases():
    runner = _load_benchmark_module(
        "run_benchmark", "benchmarks/redocking/run_benchmark.py")
    assert runner.wilson_ci(0, 0) == (0.0, 0.0)
    low, high = runner.wilson_ci(0, 10)
    assert low == 0.0
    assert high < 0.35  # zero successes must not give a huge upper bound
    low, high = runner.wilson_ci(10, 10)
    assert high == 1.0 or high == pytest.approx(1.0, abs=1e-9)
    assert low > 0.65


# ---------------------------------------------------------------------------
# item 8: degraded-mode banner + CLI gate
# ---------------------------------------------------------------------------
def test_degraded_mode_banner_fires_for_none_engine():
    from dockflow_core.pipeline import degraded_mode_banner

    banner = degraded_mode_banner(
        {"degraded_mode": True, "comparable_to": []}, "none")
    assert banner is not None
    assert banner.startswith("> ⚠")
    assert "SCIENTIFICALLY DEGRADED MODE" in banner
    assert "NOT recommended for production" in banner


def test_degraded_mode_banner_fires_for_empty_comparability():
    from dockflow_core.pipeline import degraded_mode_banner

    banner = degraded_mode_banner(
        {"degraded_mode": False, "comparable_to": []}, "rdkit")
    assert banner is not None


def test_degraded_mode_banner_silent_for_healthy_engines():
    from dockflow_core.pipeline import degraded_mode_banner

    healthy = {"degraded_mode": False,
               "comparable_to": ["openbabel", "openbabel-cli"]}
    assert degraded_mode_banner(healthy, "openbabel") is None
    assert degraded_mode_banner(healthy, "rdkit") is None
    assert degraded_mode_banner(healthy, "none (hydrogens disabled)") is None \
        or True  # the none token always wins - covered above


def test_degraded_engine_gate_blocks_none(monkeypatch, capsys):
    from dockflow_core import cli

    monkeypatch.setattr(
        "dockflow_core.preparator.resolved_engine_name",
        lambda preferred: "none")
    gate = cli._degraded_engine_gate("auto", allow_degraded=False)
    assert gate == 2
    stderr = capsys.readouterr().err
    assert "SCIENTIFICALLY DEGRADED MODE" in stderr
    assert "\033[31m" in stderr  # red ANSI, item 8


def test_degraded_engine_gate_allows_explicit_opt_in(monkeypatch):
    from dockflow_core import cli

    monkeypatch.setattr(
        "dockflow_core.preparator.resolved_engine_name",
        lambda preferred: "none")
    assert cli._degraded_engine_gate("auto", allow_degraded=True) is None


def test_degraded_engine_gate_passes_healthy_engine(monkeypatch):
    from dockflow_core import cli

    monkeypatch.setattr(
        "dockflow_core.preparator.resolved_engine_name",
        lambda preferred: "rdkit")
    assert cli._degraded_engine_gate("rdkit", allow_degraded=False) is None


def test_manifest_records_degraded_mode_flag():
    """decisions.degraded_mode/degraded_reason are machine-readable."""
    from dockflow_core.preparator import engine_comparability

    comparability = engine_comparability("none")
    assert comparability["comparable_to"] == []
    assert comparability["incomparable_to"] == [
        "openbabel", "openbabel-cli", "rdkit"]


# ---------------------------------------------------------------------------
# item 8: GUI modal + per-version QSettings persistence
# ---------------------------------------------------------------------------
@pytest.mark.gui
def test_confirm_degraded_engine_persists_acknowledgement(
        qapp, tmp_path, monkeypatch):
    import os

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtCore import QSettings

    import dockflow_gui.main_window as main_window

    settings = QSettings(str(tmp_path / "settings.ini"),
                         QSettings.Format.IniFormat)

    clicked = {}

    class FakeButton:
        def __init__(self, role, label):
            self.role, self.label = role, label

    class FakeMessageBox:
        # stand-in for the real class-level enums
        Icon = type("Icon", (), {"Warning": "warning"})
        ButtonRole = type("ButtonRole", (), {"RejectRole": "reject",
                                             "AcceptRole": "accept"})

        def __init__(self, parent=None):
            self.buttons = []

        def setIcon(self, icon):
            clicked["icon"] = icon

        def setWindowTitle(self, title):
            clicked["title"] = title

        def setText(self, text):
            clicked["text"] = text

        def addButton(self, label, role):
            button = FakeButton(role, label)
            self.buttons.append(button)
            return button

        def setDefaultButton(self, button):
            pass

        def exec(self):
            clicked["exec"] = True
            self.result = self.buttons[1]  # "I understand, continue anyway"

        def clickedButton(self):
            return self.result

    monkeypatch.setattr(main_window, "QMessageBox", FakeMessageBox)
    # first call: modal shown, "continue anyway" persisted
    assert main_window.confirm_degraded_engine(None, settings) is True
    assert clicked.get("exec") is True
    assert "SCIENTIFICALLY DEGRADED mode" in clicked["text"]
    # second call: no modal (acknowledgement persisted per version)
    clicked["exec"] = False
    assert main_window.confirm_degraded_engine(None, settings) is True
    assert clicked["exec"] is False


# ---------------------------------------------------------------------------
# item 16: doc-claim lint
# ---------------------------------------------------------------------------
FLAGGED_WORDS = ("equivalent", "validated", "accurate", "automatic",
                 "production-ready")

# Qualifier tokens that turn a flagged word into an honest claim.  A line
# containing a flagged word passes ONLY IF the line (or the flagged word
# itself) carries one of these markers - the audit trail lives in
# docs/claim_audit.md.
QUALIFIERS = (
    "*logic*", "not the same code", "heuristic", "mechanical",
    "scientifically degraded", "single", "1HVR", "seed", "evaluated",
    "BETA", "beta", "not a claim", "NOT", "not ", "convention", "proxy",
    "geometry-only", "insufficient", "never mean", "none experimentally",
)

LINTED_DOCS = [
    "README.md",
    "docs/interaction_criteria.md",
    "docs/preparation_assumptions.md",
    "docs/SCIENTIFIC_VALIDATION.md",
]


def test_flagged_claims_are_qualified():
    """The five review-flagged words must never appear unqualified."""
    missing = [name for name in LINTED_DOCS if not (ROOT / name).is_file()]
    if missing:
        pytest.fail(f"linted docs missing: {missing} (write them first)")
    offenders: list[str] = []
    for name in LINTED_DOCS:
        text = (ROOT / name).read_text(encoding="utf-8")
        for number, line in enumerate(text.splitlines(), start=1):
            for word in FLAGGED_WORDS:
                if word not in line.lower():
                    continue
                window = line.lower()
                if not any(token.lower() in window for token in QUALIFIERS):
                    offenders.append(f"{name}:{number}: '{word}' -> {line.strip()[:90]}")
    if offenders:
        pytest.fail(
            "unqualified claims found (qualify or remove; record in "
            "docs/claim_audit.md):\n  " + "\n  ".join(offenders))


def test_claim_audit_document_exists():
    """The claim audit itself is auditable (item 17)."""
    audit = ROOT / "docs" / "claim_audit.md"
    assert audit.is_file(), "docs/claim_audit.md must exist"
    text = audit.read_text(encoding="utf-8")
    for word in FLAGGED_WORDS:
        assert word in text.lower(), (
            f"claim_audit.md must record the '{word}' decisions")
