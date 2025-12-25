"""
Utility module for applying card-like styling to Streamlit components.
"""
import streamlit as st


def apply_card_style():
    """
    Apply global card styling to make containers feel like cards with:
    - Soft shadows
    - Rounded corners
    - Proper padding
    - Clean borders
    - 3D effect when opened
    
    This should be called once at the beginning of a page that uses cards.
    """
    card_css = """
    <style>
    /* Style for expander cards - closed state */
    .streamlit-expanderHeader {
        background-color: white;
        border-radius: 12px;
        border: 1px solid #e0e0e0;
        box-shadow: 0 3px 8px rgba(0,0,0,0.10);
        transition: all 0.3s ease;
    }
    
    .streamlit-expanderHeader:hover {
        box-shadow: 0 6px 16px rgba(0,0,0,0.16);
        border-color: #c0c0c0;
        transform: translateY(-1px);
    }
    
    /* Style for expanded card content - 3D effect with prominent shadows */
    div[data-testid="stExpander"] > details[open] > div {
        background-color: white;
        border-radius: 0 0 12px 12px;
        border: 2px solid #b0b0b0;
        border-top: none;
        box-shadow: 0 12px 32px rgba(0,0,0,0.18), 0 6px 12px rgba(0,0,0,0.12);
        padding: 1.2rem;
    }
    
    /* Enhanced header style when expanded */
    div[data-testid="stExpander"] > details[open] > summary {
        border: 2px solid #b0b0b0;
        border-bottom: 1px solid #e0e0e0;
        box-shadow: 0 4px 12px rgba(0,0,0,0.14);
        border-radius: 12px 12px 0 0 !important;
    }
    
    /* Remove default Streamlit expander styling conflicts */
    div[data-testid="stExpander"] {
        border: none !important;
        box-shadow: none !important;
    }
    
    /* Add spacing between cards */
    div[data-testid="stExpander"] {
        margin-bottom: 1rem;
    }
    </style>
    """
    st.markdown(card_css, unsafe_allow_html=True)


def create_custom_card(content_html, card_class="custom-card"):
    """
    Create a custom card with specified HTML content.
    
    Args:
        content_html (str): The HTML content to display inside the card
        card_class (str): CSS class name for the card (for custom styling if needed)
    
    Returns:
        None (renders the card via st.markdown)
    """
    card_style = f"""
    <style>
    .{card_class} {{
        background: white;
        border-radius: 16px;
        padding: 1.2rem;
        box-shadow: 0 4px 14px rgba(0,0,0,0.06);
        border: 1px solid #eee;
        margin-bottom: 1rem;
        transition: all 0.2s ease;
    }}
    
    .{card_class}:hover {{
        box-shadow: 0 6px 20px rgba(0,0,0,0.1);
        border-color: #ddd;
    }}
    </style>
    """
    
    card_html = f"""
    {card_style}
    <div class="{card_class}">
        {content_html}
    </div>
    """
    
    st.markdown(card_html, unsafe_allow_html=True)


def get_card_container_style():
    """
    Returns inline style string for creating card-like containers.
    Useful for wrapping content in card styling without using expanders.
    
    Returns:
        str: CSS style string
    """
    return """
        background: white;
        border-radius: 16px;
        padding: 1.2rem;
        box-shadow: 0 4px 14px rgba(0,0,0,0.06);
        border: 1px solid #eee;
        margin-bottom: 1rem;
    """
