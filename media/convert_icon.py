#!/usr/bin/env python3
"""
Convert CCOM.png to CCOM.ico for use with PyInstaller
"""

import os
from PIL import Image

def convert_png_to_ico():
    """Convert CCOM.png to CCOM.ico with multiple sizes"""
    
    png_path = 'ccom_logo_no_text_dark_background.png'
    ico_path = 'CCOM.ico'
    
    if not os.path.exists(png_path):
        print(f"Error: {png_path} not found")
        return False
    
    try:
        # Open the PNG image
        img = Image.open(png_path)
        
        # Convert to RGBA if not already
        if img.mode != 'RGBA':
            img = img.convert('RGBA')
        
        # Create multiple sizes for the ICO file (Windows standard sizes)
        sizes = [(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)]
        
        # Resize image for each size
        icon_images = []
        for size in sizes:
            resized = img.resize(size, Image.Resampling.LANCZOS)
            icon_images.append(resized)
        
        # Save as ICO file
        icon_images[0].save(ico_path, format='ICO', sizes=[(img.width, img.height) for img in icon_images])
        
        print(f"✓ Successfully converted {png_path} to {ico_path}")
        print(f"✓ Icon file size: {os.path.getsize(ico_path) / 1024:.1f} KB")
        
        return True
        
    except Exception as e:
        print(f"✗ Error converting icon: {e}")
        return False

if __name__ == "__main__":
    print("Converting CCOM.png to CCOM.ico...")
    success = convert_png_to_ico()
    
    if success:
        print("Icon conversion completed successfully!")
    else:
        print("Icon conversion failed!")




