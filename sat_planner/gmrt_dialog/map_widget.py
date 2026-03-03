# Copyright (c) 2026 Paul Johnson
# SPDX-License-Identifier: BSD-3-Clause
# Embedded in SAT Planner.

"""
Custom widget for displaying map with rectangle drawing and coordinate display.
"""

import math
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import Qt, QRect, QPoint, pyqtSignal
from PyQt6.QtGui import QPixmap, QPainter, QPen, QColor


class MapWidget(QWidget):
    """Custom widget for displaying map with rectangle drawing functionality."""
    bounds_selected = pyqtSignal(tuple)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.pixmap = None
        self.drawing_mode = False
        self.selection_rect = None
        self.drag_start = None
        self.current_bounds = None
        self.image_rect = None
        self.cursor_pos = None
        self.setMinimumSize(600, 400)
        self.setStyleSheet("QWidget { border: 1px solid gray; background-color: #2d2d2d; }")
        self.setMouseTracking(True)

    def set_pixmap(self, pixmap):
        self.pixmap = pixmap
        self.update()

    def set_bounds(self, west, east, south, north):
        self.current_bounds = (west, east, south, north)

    def enable_drawing(self, enabled):
        self.drawing_mode = enabled
        if not enabled:
            self.selection_rect = None
            self.drag_start = None
        self.update()

    def clear_selection(self):
        self.selection_rect = None
        self.drag_start = None
        self.update()

    def get_selection_bounds(self):
        if self.selection_rect is None or self.current_bounds is None or self.map_content_rect is None:
            return None
        west, east, south, north = self.current_bounds
        image_sel_rect = QRect(
            self.selection_rect.x() - self.map_content_rect.x(),
            self.selection_rect.y() - self.map_content_rect.y(),
            self.selection_rect.width(),
            self.selection_rect.height()
        )
        if (image_sel_rect.x() < 0 or image_sel_rect.y() < 0 or
            image_sel_rect.right() > self.map_content_rect.width() or
            image_sel_rect.bottom() > self.map_content_rect.height()):
            return None
        img_w = self.map_content_rect.width()
        img_h = self.map_content_rect.height()
        lon_per_pixel = (east - west) / (img_w - 1)
        sel_west = west + (image_sel_rect.x() * lon_per_pixel)
        sel_east = west + (image_sel_rect.x() + image_sel_rect.width()) * lon_per_pixel

        def lat_to_merc(lat):
            return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
        def merc_to_lat(merc):
            return math.degrees(2 * math.atan(math.exp(merc)) - math.pi / 2)
        merc_north = lat_to_merc(north)
        merc_south = lat_to_merc(south)
        y0 = image_sel_rect.y()
        y1 = image_sel_rect.y() + image_sel_rect.height()
        merc_y0 = merc_north - (y0 / (img_h - 1)) * (merc_north - merc_south)
        merc_y1 = merc_north - (y1 / (img_h - 1)) * (merc_north - merc_south)
        sel_north = merc_to_lat(merc_y0)
        sel_south = merc_to_lat(merc_y1)
        return (sel_west, sel_east, sel_south, sel_north)

    def get_cursor_coordinates(self, mouse_pos):
        if self.current_bounds is None or self.map_content_rect is None:
            return None
        if not self.map_content_rect.contains(mouse_pos):
            return None
        west, east, south, north = self.current_bounds
        image_x = mouse_pos.x() - self.map_content_rect.x()
        image_y = mouse_pos.y() - self.map_content_rect.y()
        img_w = self.map_content_rect.width()
        img_h = self.map_content_rect.height()
        lon_per_pixel = (east - west) / (img_w - 1)
        longitude = west + (image_x * lon_per_pixel)

        def lat_to_merc(lat):
            return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
        def merc_to_lat(merc):
            return math.degrees(2 * math.atan(math.exp(merc)) - math.pi / 2)
        merc_north = lat_to_merc(north)
        merc_south = lat_to_merc(south)
        merc_y = merc_north - (image_y / (img_h - 1)) * (merc_north - merc_south)
        latitude = merc_to_lat(merc_y)
        return (longitude, latitude)

    def mousePressEvent(self, event):
        if self.pixmap is None or self.map_content_rect is None:
            return
        if event.button() == Qt.MouseButton.LeftButton and self.map_content_rect.contains(event.pos()):
            self.drag_start = event.pos()
            self.selection_rect = QRect(self.drag_start, self.drag_start)
            self.update()

    def mouseMoveEvent(self, event):
        self.cursor_pos = event.pos()
        if getattr(self, "map_content_rect", None) is not None and self.map_content_rect.contains(event.pos()):
            self.setCursor(Qt.CursorShape.CrossCursor)
        else:
            self.setCursor(Qt.CursorShape.ArrowCursor)
        self.update()
        if self.drag_start is None or self.pixmap is None or self.map_content_rect is None:
            return
        if not self.map_content_rect.contains(event.pos()):
            return
        constrained_pos = QPoint(
            max(self.map_content_rect.left(), min(event.pos().x(), self.map_content_rect.right())),
            max(self.map_content_rect.top(), min(event.pos().y(), self.map_content_rect.bottom()))
        )
        self.selection_rect = QRect(self.drag_start, constrained_pos).normalized()
        self.update()

    def mouseReleaseEvent(self, event):
        if self.drag_start is None or self.pixmap is None or self.map_content_rect is None:
            return
        if event.button() == Qt.MouseButton.LeftButton:
            constrained_pos = QPoint(
                max(self.map_content_rect.left(), min(event.pos().x(), self.map_content_rect.right())),
                max(self.map_content_rect.top(), min(event.pos().y(), self.map_content_rect.bottom()))
            )
            self.selection_rect = QRect(self.drag_start, constrained_pos).normalized()
            self.drag_start = None
            self.update()
            bounds = self.get_selection_bounds()
            if bounds:
                self.bounds_selected.emit(bounds)

    def paintEvent(self, event):
        painter = QPainter(self)
        if self.pixmap and not self.pixmap.isNull():
            scaled_pixmap = self.pixmap.scaled(
                self.size(),
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation
            )
            x = (self.width() - scaled_pixmap.width()) // 2
            y = (self.height() - scaled_pixmap.height()) // 2
            painter.drawPixmap(x, y, scaled_pixmap)
            self.image_rect = QRect(x, y, scaled_pixmap.width(), scaled_pixmap.height())
            self.map_content_rect = QRect(
                self.image_rect.left() + 33,
                self.image_rect.top() + 13,
                self.image_rect.width() - 66,
                self.image_rect.height() - 26
            )
        else:
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             "Click 'Refresh Map' to load preview")
            self.image_rect = None
            self.map_content_rect = None
        if self.selection_rect is not None:
            pen = QPen(QColor(255, 0, 0), 2)
            painter.setPen(pen)
            painter.drawRect(self.selection_rect)
            painter.setBrush(QColor(255, 0, 0, 50))
            painter.drawRect(self.selection_rect)
        if self.cursor_pos is not None and self.map_content_rect is not None and self.map_content_rect.contains(self.cursor_pos):
            coords = self.get_cursor_coordinates(self.cursor_pos)
            if coords:
                longitude, latitude = coords
                coord_text = f"Lon: {longitude:.4f}°\nLat: {latitude:.4f}°"
                font = painter.font()
                font.setPointSize(9)
                painter.setFont(font)
                text_rect = painter.boundingRect(QRect(), Qt.TextFlag.TextDontClip, coord_text)
                text_rect.setWidth(text_rect.width() + 10)
                text_rect.setHeight(text_rect.height() + 6)
                text_rect.moveTo(self.map_content_rect.left() + 10, self.map_content_rect.bottom() - text_rect.height() - 10)
                painter.setBrush(QColor(0, 0, 0, 180))
                painter.setPen(QPen(QColor(255, 255, 255)))
                painter.drawRect(text_rect)
                painter.setPen(QColor(255, 255, 255))
                painter.drawText(text_rect, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop, coord_text)
