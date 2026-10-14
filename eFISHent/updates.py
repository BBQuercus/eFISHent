"""Update notifications and self-update."""

import json
import logging
import shutil
import subprocess
import sys
import time
from pathlib import Path

from . import __version__
from .dependencies import _efishent_home, _ensure_deps_on_path, check_all_dependencies

_PYPI_JSON_URL = "https://pypi.org/pypi/efishent/json"
_UPDATE_CHECK_INTERVAL = 86400  # 24 hours


def _notify_if_newer(latest: str) -> None:
    """Print an update notice if ``latest`` is newer than the running version."""
    if not latest or latest == __version__:
        return

    from packaging.version import Version

    if Version(latest) > Version(__version__):
        print(
            f"\n  Update available: v{__version__} -> v{latest}"
            f"  (run 'efishent --update' to upgrade)\n"
        )


def _cached_latest_version(cache_file: Path) -> str | None:
    """Return the cached latest version, or None if the cache is missing or stale."""
    if not cache_file.exists():
        return None
    with open(cache_file) as f:
        cache = json.load(f)
    last_check = cache.get("timestamp", 0)
    if time.time() - last_check >= _UPDATE_CHECK_INTERVAL:
        return None
    return cache.get("latest_version", "")


def _https_get_json(url: str) -> dict:
    """Fetch a JSON document over https; other URL schemes are refused."""
    import urllib.parse
    import urllib.request

    if urllib.parse.urlparse(url).scheme != "https":
        raise ValueError(f"Refusing to open non-https URL: {url}")

    # Only https handlers are installed, so redirects to other schemes fail too
    opener = urllib.request.OpenerDirector()
    for handler in (
        urllib.request.ProxyHandler(),
        urllib.request.HTTPSHandler(),
        urllib.request.HTTPRedirectHandler(),
        urllib.request.HTTPDefaultErrorHandler(),
        urllib.request.HTTPErrorProcessor(),
        urllib.request.UnknownHandler(),
    ):
        opener.add_handler(handler)

    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with opener.open(req, timeout=5) as resp:
        return json.loads(resp.read().decode())


def _run_update_check(cache_file: Path) -> None:
    """Notify about a newer release, querying PyPI unless a fresh cache exists."""
    try:
        # Check cache first (skip if checked within last 24 hours)
        latest = _cached_latest_version(cache_file)
        if latest is None:
            data = _https_get_json(_PYPI_JSON_URL)
            latest = data.get("info", {}).get("version", "")

            cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_file, "w") as f:
                json.dump({"timestamp": time.time(), "latest_version": latest}, f)

        _notify_if_newer(latest)
    except Exception as e:
        # Never fail the pipeline for an update check
        logging.getLogger("custom-logger").debug(f"Update check failed: {e}")


def check_for_updates() -> None:
    """Check PyPI for a newer version and notify the user.

    Non-blocking: runs in a background thread so it doesn't slow startup.
    Caches the result for 24 hours to avoid repeated network requests.
    """
    import threading

    cache_file = _efishent_home() / "update_check.json"
    thread = threading.Thread(target=_run_update_check, args=(cache_file,), daemon=True)
    thread.start()


_INSTALL_SCRIPT_URL = (
    "https://raw.githubusercontent.com/BBQuercus/eFISHent/main/install.sh"
)


def _upgrade_package(python: str) -> None:
    """Upgrade the eFISHent package in place, exiting on failure."""
    if shutil.which("uv"):
        cmd = ["uv", "pip", "install", "--python", python, "--upgrade", "efishent"]
    else:
        cmd = [python, "-m", "pip", "install", "--upgrade", "efishent"]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except Exception as e:
        print(f"Package update failed: {e}")
        sys.exit(1)
    if result.returncode != 0:
        print(f"Package update failed:\n{result.stderr.strip()}")
        sys.exit(1)


def _installed_version(python: str) -> str:
    """Return the eFISHent version importable by ``python``."""
    try:
        new_version_result = subprocess.run(
            [python, "-c", "from eFISHent import __version__; print(__version__)"],
            capture_output=True, text=True, timeout=10,
        )
        return new_version_result.stdout.strip()
    except Exception:
        return "unknown"


def _download_command(url: str) -> list[str] | None:
    """Return a curl/wget command that writes ``url`` to stdout, if available."""
    if shutil.which("curl"):
        return ["curl", "-fsSL", url]
    if shutil.which("wget"):
        return ["wget", "-qO-", url]
    return None


def _run_install_script(dl_cmd: list[str], prefix: Path) -> None:
    """Download install.sh and re-run it with --deps-only to update/verify deps."""
    dl_result = subprocess.run(dl_cmd, capture_output=True, text=True, timeout=30)
    if dl_result.returncode != 0:
        print("Warning: could not download install script.")
        return
    result = subprocess.run(
        ["sh", "-s", "--", "--deps-only", "--no-modify-rc",
         "--prefix", str(prefix)],
        input=dl_result.stdout, text=True, timeout=600,
    )
    if result.returncode != 0:
        print("Warning: dependency update had issues.")


def _update_dependencies(prefix: Path) -> None:
    """Update dependencies via install.sh if eFISHent was installed with it."""
    if not prefix.is_dir():
        print(f"No install directory at {prefix}, skipping dependency update.")
        print("Run the full installer to set up dependencies:")
        print(f"  curl -LsSf {_INSTALL_SCRIPT_URL} | sh")
        return

    dl_cmd = _download_command(_INSTALL_SCRIPT_URL)
    if not dl_cmd:
        print("Warning: neither curl nor wget found, skipping dependency update.")
        return

    try:
        _run_install_script(dl_cmd, prefix)
    except Exception as e:
        print(f"Warning: dependency update skipped ({e})")


def _print_dependency_summary() -> None:
    """Print which dependencies are found and how to install missing ones."""
    _ensure_deps_on_path()
    results = check_all_dependencies()
    all_found = all(r["found"] for r in results.values())
    for name, info in results.items():
        if info["found"]:
            version = info.get("version", "")
            print(f"  OK  {name}: {version}")
        else:
            print(f"  MISS {name}: not found ({info['needed_for']})")

    if not all_found:
        print("\nSome dependencies are missing. Run the full installer:")
        print(f"  curl -LsSf {_INSTALL_SCRIPT_URL} | sh")


def self_update() -> None:
    """Update eFISHent and its dependencies to the latest version."""
    python = sys.executable
    current_version = __version__

    print(f"Current version: {current_version}")

    # Step 1: Update the Python package
    print("\nUpdating eFISHent package...")
    _upgrade_package(python)

    new_version = _installed_version(python)
    if new_version == current_version:
        print(f"Already up to date (v{current_version}).")
    else:
        print(f"Updated: v{current_version} -> v{new_version}")

    # Step 2: Update dependencies via install.sh
    print("\nChecking dependencies...")
    _update_dependencies(_efishent_home())

    # Step 3: Verify dependencies
    print()
    _print_dependency_summary()
