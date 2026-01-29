"""
Qt UI helpers for SAT Planner: message boxes and simple dialogs.
Take parent (QWidget) as first argument; safe to call from app or any widget.
"""
from PyQt6.QtWidgets import QMessageBox


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
