"""File-based cache for external API responses."""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Any, Optional
from urllib import error as urlerror
from urllib import request as urlrequest

DEFAULT_CACHE_DIR = Path(__file__).resolve().parent / "cache" / "api"


def cache_key(url: str, extra: str = "") -> str:
    raw = f"{url}|{extra}".encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def cache_path(cache_dir: Path, key: str) -> Path:
    return cache_dir / f"{key}.json"


def load_cached(cache_dir: Path, key: str) -> Optional[dict]:
    path = cache_path(cache_dir, key)
    if not path.exists():
        return None
    try:
        with path.open(encoding="utf-8") as f:
            payload = json.load(f)
        if isinstance(payload, dict) and "body" in payload:
            return payload
    except (json.JSONDecodeError, OSError):
        return None
    return None


def save_cached(cache_dir: Path, key: str, url: str, status: int, body: Any) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "url": url,
        "status": status,
        "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "body": body,
    }
    with cache_path(cache_dir, key).open("w", encoding="utf-8") as f:
        json.dump(payload, f)


def fetch_json(
    url: str,
    *,
    cache_dir: Path,
    force_refresh: bool = False,
    timeout: int = 45,
    headers: Optional[dict] = None,
    extra_key: str = "",
) -> tuple[Optional[dict], bool]:
    """
    Fetch JSON from URL with disk cache.
    Returns (body, from_cache).
    """
    key = cache_key(url, extra_key)
    if not force_refresh:
        cached = load_cached(cache_dir, key)
        if cached is not None:
            return cached.get("body"), True

    hdrs = {"User-Agent": "bacdive-genome-enrichment/0.1 (+research)"}
    if headers:
        hdrs.update(headers)
    req = urlrequest.Request(url, headers=hdrs)
    try:
        with urlrequest.urlopen(req, timeout=timeout) as resp:
            status = resp.getcode()
            raw = resp.read().decode("utf-8")
    except (urlerror.URLError, TimeoutError) as exc:
        save_cached(cache_dir, key, url, 0, {"error": str(exc)})
        return None, False

    try:
        body = json.loads(raw) if raw.strip() else {}
    except json.JSONDecodeError:
        body = {"error": "invalid_json", "raw_preview": raw[:500]}
        status = status or 500

    save_cached(cache_dir, key, url, status, body)
    return body, False
