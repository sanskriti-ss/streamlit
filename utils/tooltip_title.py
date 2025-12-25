"""
Utility module for displaying section titles with hover tooltips showing sample images.
"""
import streamlit as st
import os
from pathlib import Path
from PIL import Image


def display_title_with_tooltip(title_text, sample_image_filename=None, description_text=None):
    """
    Display a section title with an optional hover tooltip showing a sample screenshot.
    Clicking the expand button will show the full-size image.
    
    Args:
        title_text (str): The main title text to display
        sample_image_filename (str, optional): Filename of the sample image in the sample_plots folder
        description_text (str, optional): Additional description text to show below the title
    
    Returns:
        None
    """
    # If a sample image is provided, add a tooltip and click-to-expand functionality
    if sample_image_filename:
        # Get the path to the sample image
        base_dir = Path(__file__).parent.parent
        sample_image_path = base_dir / "sample_plots" / sample_image_filename
        
        if sample_image_path.exists():
            # Generate a unique ID for this specific tooltip to avoid conflicts
            import hashlib
            unique_id = hashlib.md5(f"{title_text}{sample_image_filename}".encode()).hexdigest()[:8]
            
            # Initialize session state for this modal
            modal_key = f"modal_{unique_id}"
            if modal_key not in st.session_state:
                st.session_state[modal_key] = False
            
            # CSS for tooltip styling
            tooltip_css = f"""
            <style>
            .tooltip-container-{unique_id} {{
                position: relative;
                display: inline-block;
                margin-bottom: 0.5em;
            }}
            
            .tooltip-container-{unique_id} .tooltip-title {{
                font-size: 1.75em;
                font-weight: 600;
                cursor: help;
                border-bottom: 2px dotted #666;
            }}
            
            .tooltip-container-{unique_id} .tooltip-image {{
                visibility: hidden;
                position: absolute;
                z-index: 1000;
                left: 50%;
                transform: translateX(-50%);
                top: 100%;
                margin-top: 10px;
                background-color: white;
                border: 2px solid #ddd;
                border-radius: 8px;
                padding: 10px;
                box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
                max-width: 500px;
                opacity: 0;
                transition: opacity 0.3s, visibility 0.3s;
            }}
            
            .tooltip-container-{unique_id}:hover .tooltip-image {{
                visibility: visible;
                opacity: 1;
            }}
            
            .tooltip-container-{unique_id} .tooltip-image img {{
                max-width: 100%;
                height: auto;
                display: block;
                border-radius: 4px;
            }}
            
            .tooltip-label {{
                font-size: 0.8em;
                color: #666;
                font-style: italic;
                margin-left: 10px;
            }}
            
            .expand-hint {{
                font-size: 0.7em;
                color: #999;
                text-align: center;
                margin-top: 5px;
                font-weight: normal;
            }}
            </style>
            """
            
            # Read the image and convert to base64 for embedding
            import base64
            with open(sample_image_path, "rb") as img_file:
                img_data = base64.b64encode(img_file.read()).decode()
            
            # Determine the MIME type based on file extension
            file_ext = sample_image_path.suffix.lower()
            mime_type = "image/png" if file_ext == ".png" else "image/jpeg"
            
            # Create HTML with tooltip
            tooltip_html = f"""
            {tooltip_css}
            <div class="tooltip-container-{unique_id}">
                <span class="tooltip-title">{title_text}</span>
                <span class="tooltip-label">(hover to see example)</span>
                <div class="tooltip-image">
                    <img src="data:{mime_type};base64,{img_data}" alt="Sample plot">
                    <div class="expand-hint">↓ Click button to expand in new tab ↓</div>
                </div>
            </div>
            """
            
            # Create layout with title and expand button
            col1, col2 = st.columns([5, 1])
            with col1:
                st.markdown(tooltip_html, unsafe_allow_html=True)
            with col2:
                st.markdown("<div style='margin-top: 0.8em;'></div>", unsafe_allow_html=True)
                # Use HTML component that opens image in new tab
                st.components.v1.html(
                    f"""
                    <style>
                    .expand-btn {{
                        background-color: #ffffff;
                        border: 1px solid #d3d3d3;
                        border-radius: 4px;
                        padding: 0.25rem 0.75rem;
                        cursor: pointer;
                        font-size: 14px;
                        color: #262730;
                        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
                        font-weight: 400;
                        transition: all 0.2s ease;
                    }}
                    .expand-btn:hover {{
                        background-color: #f0f2f6;
                        border-color: #999;
                        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
                    }}
                    .expand-btn:active {{
                        background-color: #e0e2e6;
                        transform: translateY(1px);
                    }}
                    </style>
                    <button onclick="openImage()" class="expand-btn">
                        Expand
                    </button>
                    <script>
                    function openImage() {{
                        var image = new Image();
                        image.src = 'data:{mime_type};base64,{img_data}';
                        var w = window.open('');
                        w.document.write(image.outerHTML);
                        w.document.close();
                    }}
                    </script>
                    """,
                    height=40
                )
        else:
            # If image doesn't exist, just show the title with reduced size
            title_html = f'<h1 style="font-size: 1.75em; margin-bottom: 0.5em;">{title_text}</h1>'
            st.markdown(title_html, unsafe_allow_html=True)
    else:
        # No tooltip, just show the title with reduced size
        title_html = f'<h1 style="font-size: 1.75em; margin-bottom: 0.5em;">{title_text}</h1>'
        st.markdown(title_html, unsafe_allow_html=True)
    
    # Display optional description text
    if description_text:
        st.write(description_text)
