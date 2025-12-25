"""
Centralized color system for consistent, muted, and professional color palette across the app.
Colors are desaturated and softened for a more sophisticated look.
"""

# Muted, professional color palette
# Reduced saturation by ~15-20% from original bright colors
COLORS = {
    # Category colors (softer, desaturated versions)
    "Production": {
        "primary": "#7B9E89",      # Softer sage green (was #4c72b0 blue, now coordinated)
        "pastel": "#B8D4C6",       # Muted pastel green (was #A8E6CF)
        "light": "#D4E7DD",        # Very light sage
    },
    "Utilization": {
        "primary": "#D4858C",      # Dusty rose/muted coral (was #55a868 green)
        "pastel": "#E8B5BB",       # Softer pastel rose (was #FF8B94)
        "light": "#F2D7DB",        # Very light rose
    },
    "Resistance": {
        "primary": "#8B9FB5",      # Desaturated teal/slate blue (was #c44e52 red)
        "pastel": "#A8BFD1",       # Muted pastel blue (was #89CFF0)
        "light": "#D1DCEA",        # Very light blue
    },
    "Sensitivity": {
        "primary": "#9B8BA8",      # Muted purple (was #8172b2)
        "pastel": "#C5B8CD",       # Softer pastel purple (was #CDB4DB)
        "light": "#E0D9E6",        # Very light purple
    },
}

# Border colors (replacing harsh black outlines)
BORDERS = {
    "strong": "#4A4A4A",      # Dark gray instead of black for top ranks
    "medium": "#8A8A8A",      # Medium gray for mid ranks
    "subtle": "#D3D3D3",      # Light gray for subtle borders
    "none": "none",           # No border
}

# Neutral colors for UI elements
NEUTRALS = {
    "white": "#FFFFFF",
    "off_white": "#FAFAFA",
    "light_gray": "#F0F0F0",
    "medium_gray": "#CCCCCC",
    "dark_gray": "#666666",
    "text": "#262730",
}


def get_category_color(category, variant="primary"):
    """
    Get color for a specific category.
    
    Args:
        category (str): One of "Production", "Utilization", "Resistance", "Sensitivity"
        variant (str): One of "primary", "pastel", or "light"
    
    Returns:
        str: Hex color code
    """
    if category not in COLORS:
        return NEUTRALS["medium_gray"]
    return COLORS[category].get(variant, COLORS[category]["primary"])


def get_border_color(rank=None):
    """
    Get appropriate border color based on rank.
    Uses softer grays instead of harsh black.
    
    Args:
        rank (int or None): Numeric rank (lower is better)
    
    Returns:
        str: CSS border style string
    """
    if rank is None or not isinstance(rank, int):
        return f"1px solid {BORDERS['subtle']}"
    
    if rank <= 10:
        return f"2px solid {BORDERS['strong']}"  # Dark gray instead of black
    elif rank <= 100:
        return f"2px solid {BORDERS['medium']}"  # Medium gray
    else:
        return f"1px solid {BORDERS['subtle']}"  # Light gray


def get_category_colors_dict(variant="primary"):
    """
    Get a dictionary mapping all categories to their colors.
    Useful for bulk operations like plotting.
    
    Args:
        variant (str): One of "primary", "pastel", or "light"
    
    Returns:
        dict: Category name to color hex code
    """
    return {
        category: get_category_color(category, variant)
        for category in COLORS.keys()
    }


# Legacy compatibility: provide direct color mappings
# These match the function outputs for easy migration
CATEGORY_COLORS_PRIMARY = get_category_colors_dict("primary")
CATEGORY_COLORS_PASTEL = get_category_colors_dict("pastel")
CATEGORY_COLORS_LIGHT = get_category_colors_dict("light")
