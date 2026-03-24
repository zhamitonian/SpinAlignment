#!/usr/bin/env python3

"""
Create an HTML viewer for fitting result images
Lightweight alternative to PDF for viewing many images
"""

import os
import glob
import re
import base64
from pathlib import Path

def extract_bin_number(filename):
    """Extract bin number from filename"""
    match = re.search(r'bin_(\d+)', filename)
    if match:
        return int(match.group(1))
    return -1

def image_to_base64(image_path):
    """Convert image to base64 for embedding in HTML"""
    with open(image_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

def create_html_viewer(fit_result_dir="fit_result",
                       output_html="fit_results_viewer.html",
                       pattern="*.png",
                       embed_images=False,
                       images_per_row=3):
    """
    Create an HTML viewer for fit result images
    
    Parameters:
    -----------
    fit_result_dir : str
        Directory containing fit result images
    output_html : str
        Output HTML filename
    pattern : str
        Filename pattern to match
    embed_images : bool
        If True, embed images as base64 (larger HTML but standalone)
        If False, use relative paths (smaller HTML but needs image files)
    images_per_row : int
        Number of images per row in the grid
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
    
    image_files_with_bins.sort(key=lambda x: (x[1] if x[1] >= 0 else float('inf'), x[0]))
    
    # Generate HTML
    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Fit Results Viewer</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: Arial, sans-serif;
            background-color: #f5f5f5;
            padding: 20px;
        }
        
        .header {
            background-color: #2c3e50;
            color: white;
            padding: 20px;
            margin-bottom: 20px;
            border-radius: 5px;
            position: sticky;
            top: 0;
            z-index: 100;
            box-shadow: 0 2px 5px rgba(0,0,0,0.2);
        }
        
        .header h1 {
            margin-bottom: 10px;
        }
        
        .controls {
            display: flex;
            gap: 15px;
            align-items: center;
            flex-wrap: wrap;
        }
        
        .controls label {
            display: flex;
            align-items: center;
            gap: 5px;
        }
        
        .controls input, .controls select {
            padding: 5px 10px;
            border-radius: 3px;
            border: 1px solid #ccc;
        }
        
        #searchBox {
            width: 200px;
        }
        
        .stats {
            color: #ecf0f1;
            font-size: 14px;
        }
        
        .gallery {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }
        
        .image-card {
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            transition: transform 0.2s, box-shadow 0.2s;
        }
        
        .image-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 4px 16px rgba(0,0,0,0.2);
        }
        
        .image-card.hidden {
            display: none;
        }
        
        .image-header {
            background-color: #34495e;
            color: white;
            padding: 10px 15px;
            font-weight: bold;
        }
        
        .image-container {
            position: relative;
            cursor: pointer;
        }
        
        .image-container img {
            width: 100%;
            height: auto;
            display: block;
        }
        
        /* Modal for full-size image */
        .modal {
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.9);
        }
        
        .modal-content {
            margin: auto;
            display: block;
            max-width: 95%;
            max-height: 95%;
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
        }
        
        .close {
            position: absolute;
            top: 15px;
            right: 35px;
            color: #f1f1f1;
            font-size: 40px;
            font-weight: bold;
            cursor: pointer;
        }
        
        .close:hover {
            color: #ff0000;
        }
        
        .nav-buttons {
            position: absolute;
            top: 50%;
            transform: translateY(-50%);
            font-size: 40px;
            color: white;
            cursor: pointer;
            padding: 20px;
            user-select: none;
        }
        
        .nav-buttons:hover {
            background-color: rgba(255,255,255,0.1);
        }
        
        .prev {
            left: 20px;
        }
        
        .next {
            right: 20px;
        }
        
        @media (max-width: 768px) {
            .gallery {
                grid-template-columns: 1fr;
            }
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>Fit Results Viewer</h1>
        <div class="controls">
            <label>
                Search Bin:
                <input type="text" id="searchBox" placeholder="e.g., 247">
            </label>
            <label>
                Images per row:
                <select id="columnsSelect">
                    <option value="1">1</option>
                    <option value="2">2</option>
                    <option value="3" selected>3</option>
                    <option value="4">4</option>
                    <option value="5">5</option>
                </select>
            </label>
            <span class="stats">
                Showing <span id="visibleCount">0</span> of <span id="totalCount">0</span> images
            </span>
        </div>
    </div>
    
    <div class="gallery" id="gallery">
        <!-- Images will be inserted here -->
    </div>
    
    <!-- Modal for full-size image -->
    <div id="imageModal" class="modal">
        <span class="close">&times;</span>
        <span class="nav-buttons prev">&#10094;</span>
        <span class="nav-buttons next">&#10095;</span>
        <img class="modal-content" id="modalImage">
    </div>
    
    <script>
        // Image data
        const images = [
"""
    
    # Add image data
    for img_file, bin_num in image_files_with_bins:
        rel_path = os.path.relpath(img_file, os.path.dirname(output_html))
        
        if embed_images:
            print(f"Embedding: {img_file}")
            img_data = image_to_base64(img_file)
            img_src = f"data:image/png;base64,{img_data}"
        else:
            img_src = rel_path
        
        bin_label = f"Bin {bin_num}" if bin_num >= 0 else "Unknown"
        
        html_content += f"""            {{
                src: "{img_src}",
                bin: {bin_num},
                label: "{bin_label}",
                filename: "{os.path.basename(img_file)}"
            }},
"""
    
    html_content += """        ];
        
        const gallery = document.getElementById('gallery');
        const searchBox = document.getElementById('searchBox');
        const columnsSelect = document.getElementById('columnsSelect');
        const visibleCount = document.getElementById('visibleCount');
        const totalCount = document.getElementById('totalCount');
        const modal = document.getElementById('imageModal');
        const modalImg = document.getElementById('modalImage');
        const closeBtn = document.querySelector('.close');
        const prevBtn = document.querySelector('.prev');
        const nextBtn = document.querySelector('.next');
        
        let currentImageIndex = -1;
        let visibleImages = [];
        
        // Initialize gallery
        function renderGallery() {
            gallery.innerHTML = '';
            images.forEach((img, index) => {
                const card = document.createElement('div');
                card.className = 'image-card';
                card.dataset.bin = img.bin;
                card.dataset.index = index;
                
                card.innerHTML = `
                    <div class="image-header">${img.label}</div>
                    <div class="image-container">
                        <img src="${img.src}" alt="${img.label}" loading="lazy">
                    </div>
                `;
                
                card.querySelector('.image-container').addEventListener('click', () => {
                    openModal(index);
                });
                
                gallery.appendChild(card);
            });
            
            updateStats();
        }
        
        // Search functionality
        searchBox.addEventListener('input', (e) => {
            const searchTerm = e.target.value.toLowerCase();
            const cards = document.querySelectorAll('.image-card');
            
            cards.forEach(card => {
                const binNum = card.dataset.bin;
                if (searchTerm === '' || binNum.includes(searchTerm)) {
                    card.classList.remove('hidden');
                } else {
                    card.classList.add('hidden');
                }
            });
            
            updateStats();
        });
        
        // Column selection
        columnsSelect.addEventListener('change', (e) => {
            const columns = e.target.value;
            const minWidth = Math.floor(1200 / columns);
            gallery.style.gridTemplateColumns = `repeat(auto-fit, minmax(${minWidth}px, 1fr))`;
        });
        
        // Modal functions
        function openModal(index) {
            visibleImages = Array.from(document.querySelectorAll('.image-card:not(.hidden)'))
                .map(card => parseInt(card.dataset.index));
            currentImageIndex = visibleImages.indexOf(index);
            
            modal.style.display = 'block';
            modalImg.src = images[index].src;
        }
        
        function closeModal() {
            modal.style.display = 'none';
        }
        
        function showPrevImage() {
            if (currentImageIndex > 0) {
                currentImageIndex--;
                modalImg.src = images[visibleImages[currentImageIndex]].src;
            }
        }
        
        function showNextImage() {
            if (currentImageIndex < visibleImages.length - 1) {
                currentImageIndex++;
                modalImg.src = images[visibleImages[currentImageIndex]].src;
            }
        }
        
        function updateStats() {
            const visible = document.querySelectorAll('.image-card:not(.hidden)').length;
            visibleCount.textContent = visible;
            totalCount.textContent = images.length;
        }
        
        // Event listeners
        closeBtn.addEventListener('click', closeModal);
        prevBtn.addEventListener('click', showPrevImage);
        nextBtn.addEventListener('click', showNextImage);
        
        modal.addEventListener('click', (e) => {
            if (e.target === modal) {
                closeModal();
            }
        });
        
        document.addEventListener('keydown', (e) => {
            if (modal.style.display === 'block') {
                if (e.key === 'Escape') closeModal();
                if (e.key === 'ArrowLeft') showPrevImage();
                if (e.key === 'ArrowRight') showNextImage();
            }
        });
        
        // Initialize
        renderGallery();
    </script>
</body>
</html>
"""
    
    # Write HTML file
    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    file_size = os.path.getsize(output_html) / 1024 / 1024
    print(f"\nSuccessfully created {output_html}")
    print(f"File size: {file_size:.2f} MB")
    print(f"Total images: {len(image_files)}")
    print(f"\nOpen it in a web browser to view the results.")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Create HTML viewer for fitting result images'
    )
    
    parser.add_argument(
        '-i', '--input-dir',
        default='fit_result',
        help='Directory containing fit result images (default: fit_result)'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='fit_results_viewer.html',
        help='Output HTML filename (default: fit_results_viewer.html)'
    )
    
    parser.add_argument(
        '-p', '--pattern',
        default='*.png',
        help='Filename pattern to match (default: *.png)'
    )
    
    parser.add_argument(
        '--embed',
        action='store_true',
        help='Embed images as base64 (creates standalone HTML but larger file)'
    )
    
    parser.add_argument(
        '--columns',
        type=int,
        default=3,
        help='Default number of images per row (default: 3)'
    )
    
    args = parser.parse_args()
    
    create_html_viewer(
        fit_result_dir=args.input_dir,
        output_html=args.output,
        pattern=args.pattern,
        embed_images=args.embed,
        images_per_row=args.columns
    )

if __name__ == "__main__":
    main()
