#!/usr/bin/env python3

"""
Collect all fitting result images (fit_result/*.png) into a single PDF
with bin number annotations
"""

import os
import glob
import re
from PIL import Image, ImageDraw, ImageFont
import numpy as np

def extract_bin_number(filename):
    """
    Extract bin number from filename
    Expected format: fit_result/bin_XXX.png or similar
    
    Returns:
        int: bin number, or -1 if not found
    """
    match = re.search(r'bin_(\d+)', filename)
    if match:
        return int(match.group(1))
    return -1

def add_bin_label(image, bin_number, position='top-left', font_size=40):
    """
    Add bin number label to image
    
    Parameters:
    -----------
    image : PIL.Image
        Input image
    bin_number : int
        Bin number to display
    position : str
        Position of label: 'top-left', 'top-right', 'bottom-left', 'bottom-right'
    font_size : int
        Font size for label
        
    Returns:
    --------
    PIL.Image : Image with label added
    """
    img_with_label = image.copy()
    draw = ImageDraw.Draw(img_with_label)
    
    # Try to use a nice font, fall back to default if not available
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", font_size)
    except:
        try:
            font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf", font_size)
        except:
            font = ImageFont.load_default()
    
    label_text = f"Bin {bin_number}"
    
    # Get text bounding box
    bbox = draw.textbbox((0, 0), label_text, font=font)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    # Calculate position
    margin = 20
    img_width, img_height = img_with_label.size
    
    if position == 'top-left':
        text_pos = (margin, margin)
    elif position == 'top-right':
        text_pos = (img_width - text_width - margin, margin)
    elif position == 'bottom-left':
        text_pos = (margin, img_height - text_height - margin)
    elif position == 'bottom-right':
        text_pos = (img_width - text_width - margin, img_height - text_height - margin)
    else:
        text_pos = (margin, margin)
    
    # Draw background rectangle for better visibility
    padding = 10
    rect_coords = [
        text_pos[0] - padding,
        text_pos[1] - padding,
        text_pos[0] + text_width + padding,
        text_pos[1] + text_height + padding
    ]
    draw.rectangle(rect_coords, fill=(255, 255, 255, 230), outline=(0, 0, 0))
    
    # Draw text
    draw.text(text_pos, label_text, fill=(0, 0, 0), font=font)
    
    return img_with_label

def collect_fit_results_to_pdf(fit_result_dir="fit_result", 
                                output_pdf="all_fit_results.pdf",
                                pattern="*.png",
                                label_position='top-left'):
    """
    Collect all fitting result images into a single PDF
    
    Parameters:
    -----------
    fit_result_dir : str
        Directory containing fit result images
    output_pdf : str
        Output PDF filename
    pattern : str
        Filename pattern to match (e.g., "*.png", "bin_*.png")
    label_position : str
        Position of bin label on images
    """
    
    # Find all PNG files
    search_pattern = os.path.join(fit_result_dir, pattern)
    image_files = glob.glob(search_pattern)
    
    if not image_files:
        print(f"No images found matching pattern: {search_pattern}")
        return
    
    print(f"Found {len(image_files)} images")
    
    # Sort by bin number
    image_files_with_bins = []
    for img_file in image_files:
        bin_num = extract_bin_number(img_file)
        image_files_with_bins.append((img_file, bin_num))
    
    # Sort by bin number (-1 for files without bin number, they go to end)
    image_files_with_bins.sort(key=lambda x: (x[1] if x[1] >= 0 else float('inf'), x[0]))
    
    # Process images
    processed_images = []
    for img_file, bin_num in image_files_with_bins:
        try:
            print(f"Processing: {img_file} (Bin {bin_num})")
            img = Image.open(img_file)
            
            # Convert to RGB if necessary (in case of RGBA or other modes)
            if img.mode != 'RGB':
                img = img.convert('RGB')
            
            # Add bin label if bin number was found
            if bin_num >= 0:
                img = add_bin_label(img, bin_num, position=label_position)
            
            processed_images.append(img)
            
        except Exception as e:
            print(f"Error processing {img_file}: {e}")
            continue
    
    if not processed_images:
        print("No valid images to save")
        return
    
    # Save to PDF
    print(f"\nSaving {len(processed_images)} images to {output_pdf}...")
    
    # First image is used to create the PDF, rest are appended
    processed_images[0].save(
        output_pdf,
        save_all=True,
        append_images=processed_images[1:],
        resolution=100.0,
        quality=95
    )
    
    print(f"Successfully created {output_pdf}")
    print(f"Total pages: {len(processed_images)}")

def main():
    """
    Main function with command line argument parsing
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Collect fitting result images into a single PDF with bin labels'
    )
    
    parser.add_argument(
        '-i', '--input-dir',
        default='fit_result',
        help='Directory containing fit result images (default: fit_result)'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='all_fit_results.pdf',
        help='Output PDF filename (default: all_fit_results.pdf)'
    )
    
    parser.add_argument(
        '-p', '--pattern',
        default='*.png',
        help='Filename pattern to match (default: *.png)'
    )
    
    parser.add_argument(
        '--label-position',
        choices=['top-left', 'top-right', 'bottom-left', 'bottom-right'],
        default='top-left',
        help='Position of bin label (default: top-left)'
    )
    
    args = parser.parse_args()
    
    collect_fit_results_to_pdf(
        fit_result_dir=args.input_dir,
        output_pdf=args.output,
        pattern=args.pattern,
        label_position=args.label_position
    )

if __name__ == "__main__":
    main()
