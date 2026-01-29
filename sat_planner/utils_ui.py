"""
Qt UI helpers for SAT Planner: message boxes and simple dialogs.
Take parent (QWidget) as first argument; safe to call from app or any widget.
"""
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
