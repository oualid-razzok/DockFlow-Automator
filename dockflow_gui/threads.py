"""Qt worker threads: keep the UI responsive during long operations."""

from __future__ import annotations

import traceback
from collections.abc import Callable
from typing import Any

from PyQt6.QtCore import QObject, QThread, pyqtSignal


class WorkerSignals(QObject):
    """Signals shared by every worker (Qt signals are thread-safe)."""

    started = pyqtSignal()
    progress = pyqtSignal(int, str)          # percent, message
    log = pyqtSignal(str)                    # raw log line
    step = pyqtSignal(str, str)              # step name, status
    result = pyqtSignal(object)              # worker return value
    error = pyqtSignal(str, str)             # message, traceback text
    finished = pyqtSignal()
    cancelled = pyqtSignal()


class Worker(QThread):
    """Run a callable in a background thread.

    If the callable declares a ``cancel_event`` keyword parameter it is
    automatically fed a :class:`threading.Event` that the user can trigger
    through :meth:`cancel`.
    """

    def __init__(self, fn: Callable[..., Any], *args, parent=None, **kwargs) -> None:
        super().__init__(parent)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        self._cancelled = False
        self._cancel_event = kwargs.get("cancel_event")

    # -- control ------------------------------------------------------------
    def cancel(self) -> None:
        self._cancelled = True
        if self._cancel_event is not None:
            self._cancel_event.set()

    @property
    def cancelled(self) -> bool:
        return self._cancelled

    # -- execution ----------------------------------------------------------
    def run(self) -> None:
        import inspect
        import threading

        self.signals.started.emit()
        kwargs = dict(self.kwargs)
        try:
            signature = inspect.signature(self.fn)
            if "cancel_event" in signature.parameters and \
                    "cancel_event" not in kwargs:
                import threading

                event = threading.Event()
                kwargs["cancel_event"] = event
                self._cancel_event = event
        except (TypeError, ValueError):  # builtins / C callables
            kwargs = dict(self.kwargs)
        try:
            result = self.fn(*self.args, **kwargs)
            if self._cancelled:
                self.signals.cancelled.emit()
            else:
                self.signals.result.emit(result)
        except Exception as exc:  # noqa: BLE001 - report to the UI, never crash
            tb = traceback.format_exc()
            self.signals.error.emit(f"{type(exc).__name__}: {exc}", tb)
        finally:
            self.signals.finished.emit()


class PipelineWorker(QThread):
    """Run a full :class:`dockflow_core.pipeline.DockingPipeline` with signals."""

    def __init__(self, pipeline, parent=None) -> None:
        super().__init__(parent)
        self.pipeline = pipeline
        self.signals = WorkerSignals()

    def cancel(self) -> None:
        self.pipeline.cancel()

    def run(self) -> None:
        from dockflow_core.pipeline import PipelineEvents

        self.signals.started.emit()
        events = PipelineEvents(
            on_step=self._on_step,
            on_progress=self._on_progress,
            on_log=self._on_log,
        )
        self.pipeline.events = events
        try:
            report = self.pipeline.run()
            self.signals.result.emit(report)
        except Exception as exc:  # noqa: BLE001
            import traceback

            self.signals.error.emit(f"{type(exc).__name__}: {exc}",
                                    traceback.format_exc())
        finally:
            self.signals.finished.emit()

    def _on_step(self, step: str, status: str, detail: str | None) -> None:
        self.signals.step.emit(f"{step} {detail or ''}".strip(), status)

    def _on_progress(self, fraction: float, message: str) -> None:
        self.signals.progress.emit(int(fraction * 100), message)

    def _on_log(self, message: str) -> None:
        self.signals.log.emit(message)


class FunctionWorker(Worker):
    """Worker with convenience signal payload for (progress, log) callables.

    ``fn`` may accept ``on_progress(fraction, message)`` and
    ``on_log(message)`` keyword arguments which are bridged into Qt
    signals automatically.
    """

    def run(self) -> None:
        import inspect
        import threading

        self.signals.started.emit()
        kwargs = dict(self.kwargs)
        callbacks = {}
        try:
            signature = inspect.signature(self.fn)
            for name in ("on_progress", "on_log", "cancel_event"):
                if name in signature.parameters and name not in kwargs:
                    if name == "on_progress":
                        callbacks[name] = self._emit_progress
                    elif name == "on_log":
                        callbacks[name] = self._emit_log
                    else:
                        event = threading.Event()
                        kwargs[name] = event
                        self._cancel_event = event
        except (TypeError, ValueError):
            callbacks = {}
        kwargs.update(callbacks)
        try:
            result = self.fn(*self.args, **kwargs)
            if getattr(self, "_cancel_event", None) is not None and \
                    self._cancel_event.is_set():
                self.signals.cancelled.emit()
            else:
                self.signals.result.emit(result)
        except Exception as exc:  # noqa: BLE001
            import traceback

            self.signals.error.emit(f"{type(exc).__name__}: {exc}",
                                    traceback.format_exc())
        finally:
            self.signals.finished.emit()

    def _emit_progress(self, fraction: float, message: str) -> None:
        self.signals.progress.emit(int(max(0.0, min(1.0, fraction)) * 100), message)

    def _emit_log(self, message: str) -> None:
        self.signals.log.emit(message)

    def cancel(self) -> None:
        super().cancel()
        event = getattr(self, "_cancel_event", None)
        if event is not None:
            event.set()
