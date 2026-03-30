"""eFISHent."""

try:
    from importlib.metadata import version

    __version__ = version("efishent")
except Exception:
    __version__ = "0.0.11"
