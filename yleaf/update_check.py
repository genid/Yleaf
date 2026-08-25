"""Check whether a newer Yleaf release is available on GitHub.

Yleaf is distributed from GitHub rather than PyPI, and installations vary (git
checkout, ``pip install git+``, conda environment, bundled dashboard sidecar),
so this module never touches the installation. It only reports that a newer
release exists and leaves upgrading to the user.

The check is best-effort by design: it is bounded by a short timeout, caches its
answer so repeated runs do not hit the network, and swallows every error. A
machine that is offline, behind a proxy, or rate limited must run exactly as it
did before.

Disable it with ``--no-update-check`` or by setting ``YLEAF_NO_UPDATE_CHECK=1``.
"""

import json
import logging
import os
import time
import urllib.request
from pathlib import Path
from typing import Optional, Tuple

LOG = logging.getLogger("yleaf_logger")

RELEASES_API_URL = "https://api.github.com/repos/genid/Yleaf/releases/latest"
RELEASES_PAGE_URL = "https://github.com/genid/Yleaf/releases/latest"
REQUEST_TIMEOUT = 2  # seconds; never make a user wait on this
CACHE_MAX_AGE = 24 * 60 * 60  # re-check at most once a day
ENV_OPT_OUT = "YLEAF_NO_UPDATE_CHECK"


def _cache_file() -> Path:
    """Cache under XDG_CACHE_HOME; the package directory may be read-only."""
    base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
    return Path(base) / "yleaf" / "update_check.json"


def _read_cache() -> Optional[str]:
    """Cached latest version, or None when absent, stale or unreadable."""
    try:
        cache = json.loads(_cache_file().read_text())
        if time.time() - float(cache["checked_at"]) > CACHE_MAX_AGE:
            return None
        return cache["latest"]
    except Exception:
        return None


def _write_cache(latest: str):
    try:
        path = _cache_file()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"checked_at": time.time(), "latest": latest}))
    except Exception:
        pass  # a cache we cannot write just means we check again next time


def _fetch_latest_version() -> Optional[str]:
    """Tag name of the most recent GitHub release, without the leading 'v'."""
    request = urllib.request.Request(
        RELEASES_API_URL,
        headers={"Accept": "application/vnd.github+json", "User-Agent": "Yleaf"},
    )
    with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT) as response:
        release = json.loads(response.read().decode())
    tag = (release.get("tag_name") or "").strip()
    return tag.lstrip("vV") or None


def _version_key(version: str) -> Tuple[int, ...]:
    """Comparable key for a dotted version. Trailing non-numeric parts (rc, beta)
    are ignored, so 4.2.0-rc1 compares equal to 4.2.0 and never looks newer."""
    parts = []
    for part in version.split("."):
        digits = ""
        for char in part:
            if not char.isdigit():
                break
            digits += char
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


def is_newer(latest: str, current: str) -> bool:
    latest_key, current_key = _version_key(latest), _version_key(current)
    if not latest_key or not current_key:
        return False
    return latest_key > current_key


def get_available_update(current_version: str) -> Optional[str]:
    """Return the newer release version, or None if up to date or unreachable."""
    if os.environ.get(ENV_OPT_OUT, "").strip() not in ("", "0", "false", "False"):
        return None
    latest = _read_cache()
    if latest is None:
        try:
            latest = _fetch_latest_version()
        except Exception as error:
            LOG.debug(f"update check skipped: {error}")
            return None
        if latest is None:
            return None
        _write_cache(latest)
    return latest if is_newer(latest, current_version) else None


def log_update_notice(current_version: str):
    """Log a one-line notice when a newer release exists. Never raises."""
    try:
        latest = get_available_update(current_version)
    except Exception as error:  # belt and braces: startup must never fail here
        LOG.debug(f"update check failed: {error}")
        return
    if latest is not None:
        LOG.info(f"A new Yleaf release is available: v{current_version} -> v{latest}. "
                 f"See {RELEASES_PAGE_URL}")
