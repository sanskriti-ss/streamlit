"""Tests for API response caching."""

import json
from pathlib import Path

from genome_investigation.api_cache import cache_key, fetch_json, load_cached, save_cached


def test_cache_roundtrip(tmp_path: Path):
    url = "https://example.test/strain/1"
    key = cache_key(url)
    save_cached(tmp_path, key, url, 200, {"ok": True})
    cached = load_cached(tmp_path, key)
    assert cached is not None
    assert cached["body"] == {"ok": True}

    body, from_cache = fetch_json(url, cache_dir=tmp_path, force_refresh=False)
    assert from_cache is True
    assert body == {"ok": True}


def test_cache_key_stable():
    assert cache_key("https://a") == cache_key("https://a")
    assert cache_key("https://a") != cache_key("https://b")
