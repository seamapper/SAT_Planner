"""
Qt UI helpers for SAT Planner: message boxes and simple dialogs.
Take parent (QWidget) as first argument; safe to call from app or any widget.
"""
import datetime
import os
import sys

from PyQt6.QtWidgets import (
    QApplication,
    QMessageBox,
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QTextEdit,
    QPushButton,
    QWidget,
)
from PyQt6.QtCore import QTimer
from PyQt6.QtGui import QPalette, QColor

from sat_planner.constants import __version__


def media_path(filename):
    """Return absolute path to a file in the project media directory (source or frozen exe)."""
    if getattr(sys, "frozen", False):
        base = getattr(sys, "_MEIPASS", os.path.dirname(sys.executable))
        return os.path.join(base, "media", filename)
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(os.path.dirname(pkg_dir), "media", filename)


def get_about_program_info(source_path=None):
    """Return About-dialog fields shared by the dialog and map startup placeholder."""
    compile_date = "Unknown"
    try:
        if getattr(sys, "frozen", False):
            path = sys.executable
        else:
            path = source_path or __file__
        if path and os.path.exists(path):
            compile_date = datetime.datetime.fromtimestamp(os.path.getmtime(path)).strftime(
                "%B %d, %Y"
            )
    except Exception:
        compile_date = "Unknown"

    return {
        "title": f"UNH/CCOM-JHC SAT Planner v{__version__}",
        "compile_date": compile_date,
        "author": "Paul Johnson",
        "email": "pjohnson@ccom.unh.edu",
        "institution": (
            "Center for Coastal and Ocean Mapping/Joint Hydrographic Center, "
            "University of New Hampshire"
        ),
        "grant": (
            "This program was developed at the University of New Hampshire, "
            "Center for Coastal and Ocean Mapping - Joint Hydrographic Center "
            "(UNH/CCOM-JHC) under the grant NA25NOSX400C0001-T1-01 from the National "
            "Oceanic and Atmospheric Administration (NOAA)."
        ),
        "license": "This software is released for general use under the BSD 3-Clause License.",
        "logo_path": media_path("CCOM.png"),
    }


def apply_dark_theme(app):
    """
    Apply a dark theme to the Qt GUI (not matplotlib). Uses Fusion style and a dark
    palette so the app looks consistent in both light and dark OS themes.
    """
    app.setStyle("Fusion")
    palette = QPalette()
    # Window and base backgrounds
    palette.setColor(QPalette.ColorRole.Window, QColor(53, 53, 53))
    palette.setColor(QPalette.ColorRole.WindowText, QColor(240, 240, 240))
    palette.setColor(QPalette.ColorRole.Base, QColor(42, 42, 42))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor(53, 53, 53))
    # Buttons
    palette.setColor(QPalette.ColorRole.Button, QColor(53, 53, 53))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor(240, 240, 240))
    # Text inputs and content
    palette.setColor(QPalette.ColorRole.Text, QColor(240, 240, 240))
    palette.setColor(QPalette.ColorRole.PlaceholderText, QColor(127, 127, 127))
    # Highlights (selection)
    palette.setColor(QPalette.ColorRole.Highlight, QColor(42, 130, 218))
    palette.setColor(QPalette.ColorRole.HighlightedText, QColor(255, 255, 255))
    # Disabled state (gray text)
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.WindowText, QColor(127, 127, 127))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, QColor(127, 127, 127))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, QColor(127, 127, 127))
    # Links
    palette.setColor(QPalette.ColorRole.Link, QColor(42, 130, 218))
    app.setPalette(palette)


def show_message(parent, msg_type, title, message):
    """Show a message box. msg_type: 'error', 'warning', 'info', or 'question'."""
    msg = QMessageBox(parent)
    msg.setWindowTitle(title)
    msg.setText(message)
    if msg_type == "error":
        msg.setIcon(QMessageBox.Icon.Critical)
    elif msg_type == "warning":
        msg.setIcon(QMessageBox.Icon.Warning)
    elif msg_type == "info":
        msg.setIcon(QMessageBox.Icon.Information)
    elif msg_type == "question":
        msg.setIcon(QMessageBox.Icon.Question)
    msg.exec()
    return msg


def ask_yes_no(parent, title, message):
    """Show a yes/no question dialog. Returns True if user chose Yes."""
    reply = QMessageBox.question(
        parent, title, message,
        QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
    )
    return reply == QMessageBox.StandardButton.Yes


def ask_ok_cancel(parent, title, message):
    """Show an OK/Cancel question dialog. Returns True if user chose OK."""
    reply = QMessageBox.question(
        parent, title, message,
        QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel,
    )
    return reply == QMessageBox.StandardButton.Ok


def show_statistics_dialog(parent, title, text):
    """Create a custom dialog window with copy functionality for statistics display."""
    dialog = QDialog(parent)
    dialog.setWindowTitle(title)
    dialog.setMinimumSize(600, 500)
    dialog.setModal(True)

    main_layout = QVBoxLayout(dialog)
    main_layout.setContentsMargins(10, 10, 10, 10)

    text_widget = QTextEdit()
    text_widget.setReadOnly(True)
    text_widget.setPlainText(text)
    text_widget.setFont(QApplication.font())
    main_layout.addWidget(text_widget)

    button_frame = QWidget()
    button_layout = QHBoxLayout(button_frame)
    button_layout.setContentsMargins(0, 0, 0, 0)

    copy_btn = QPushButton("Copy to Clipboard")

    def copy_to_clipboard():
        try:
            cursor = text_widget.textCursor()
            if cursor.hasSelection():
                selected_text = cursor.selectedText()
            else:
                selected_text = text_widget.toPlainText()
            clipboard = QApplication.clipboard()
            clipboard.setText(selected_text)
            copy_btn.setText("Copied!")
            QTimer.singleShot(1000, lambda: copy_btn.setText("Copy to Clipboard"))
        except Exception as e:
            show_message(parent, "error", "Copy Error", f"Failed to copy to clipboard: {e}")

    copy_btn.clicked.connect(copy_to_clipboard)
    button_layout.addWidget(copy_btn)
    button_layout.addStretch()

    close_btn = QPushButton("Close")
    close_btn.clicked.connect(dialog.accept)
    button_layout.addWidget(close_btn)

    main_layout.addWidget(button_frame)

    dialog.exec()
