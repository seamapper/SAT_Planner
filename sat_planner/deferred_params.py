"""
Deferred QLineEdit parameters: apply on Enter / focus loss, highlight while pending.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional

from PyQt6.QtWidgets import QLineEdit

PENDING_LINEEDIT_STYLE = (
    "QLineEdit { background-color: rgb(72, 58, 32); border: 1px solid rgb(200, 150, 60); }"
)
_NORMAL_STYLE_PROPERTY = "_deferred_normal_stylesheet"


def bind_deferred_line_edit(
    widget: QLineEdit,
    applied_store: Dict[int, str],
    on_commit: Optional[Callable[[QLineEdit], None]] = None,
    bound_registry: Optional[List[QLineEdit]] = None,
) -> None:
    """Wire *widget* so side effects run on commit (Enter or editingFinished), not each keystroke."""
    key = id(widget)
    applied_store[key] = widget.text()
    if widget.property(_NORMAL_STYLE_PROPERTY) is None:
        widget.setProperty(_NORMAL_STYLE_PROPERTY, widget.styleSheet() or "")

    if bound_registry is not None and widget not in bound_registry:
        bound_registry.append(widget)

    def _refresh_style() -> None:
        if widget.signalsBlocked():
            return
        normal = widget.property(_NORMAL_STYLE_PROPERTY) or ""
        if widget.text() != applied_store.get(key):
            widget.setStyleSheet(PENDING_LINEEDIT_STYLE)
        else:
            widget.setStyleSheet(normal)

    def _commit() -> None:
        if widget.signalsBlocked():
            return
        if widget.text() == applied_store.get(key):
            return
        applied_store[key] = widget.text()
        widget.setStyleSheet(widget.property(_NORMAL_STYLE_PROPERTY) or "")
        if on_commit is not None:
            on_commit(widget)

    widget.textChanged.connect(_refresh_style)
    widget.editingFinished.connect(_commit)
    widget.returnPressed.connect(_commit)


def deferred_set_line_edit(
    widget: QLineEdit,
    text: str,
    applied_store: Dict[int, str],
    *,
    mark_applied: bool = True,
) -> None:
    """Set line-edit text from code without marking the field pending."""
    widget.blockSignals(True)
    widget.setText(text)
    widget.blockSignals(False)
    if mark_applied:
        applied_store[id(widget)] = text
        widget.setStyleSheet(widget.property(_NORMAL_STYLE_PROPERTY) or "")


def deferred_mark_applied(widget: QLineEdit, applied_store: Dict[int, str]) -> None:
    """Record current text as applied and clear pending styling."""
    applied_store[id(widget)] = widget.text()
    widget.setStyleSheet(widget.property(_NORMAL_STYLE_PROPERTY) or "")


def commit_deferred_line_edit(
    widget: QLineEdit,
    applied_store: Dict[int, str],
    on_commit: Optional[Callable[[QLineEdit], None]] = None,
) -> None:
    """Force-apply a single field if it is still pending."""
    key = id(widget)
    if widget.text() == applied_store.get(key):
        return
    applied_store[key] = widget.text()
    widget.setStyleSheet(widget.property(_NORMAL_STYLE_PROPERTY) or "")
    if on_commit is not None:
        on_commit(widget)


def commit_all_deferred_line_edits(
    widgets: List[QLineEdit],
    applied_store: Dict[int, str],
    commit_handlers: Dict[int, Callable[[QLineEdit], None]],
) -> None:
    """Commit every bound field that still differs from its last applied value."""
    for widget in widgets:
        key = id(widget)
        if widget.text() == applied_store.get(key):
            continue
        applied_store[key] = widget.text()
        widget.setStyleSheet(widget.property(_NORMAL_STYLE_PROPERTY) or "")
        handler = commit_handlers.get(key)
        if handler is not None:
            handler(widget)
