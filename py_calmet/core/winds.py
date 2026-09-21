"""Diagnostic wind construction (interp, OA, slope flow, mass consistency)."""
from __future__ import annotations
import numpy as np
from .met_utils import wind_uv, ZO_EXTRAP, layer_mids, G
# ... full content from mcp_push_winds_6500.json / local winds[:6500] ...
