"""Locate, validate, and extract the BRENDA bulk JSON dump.

BRENDA requires a click-through licence acceptance before exposing
the bulk-download URL, so the file cannot be fetched non-interactively.
The user downloads the ``brenda_<release>.json.tar.gz`` tarball once
from <https://www.brenda-enzymes.org/download.php> and drops it into
the cache directory; this module then validates, extracts, and
returns the path to the unpacked JSON.
"""
from __future__ import annotations

import hashlib
import logging
import tarfile
from pathlib import Path

logger = logging.getLogger(__name__)


class BrendaDownloadError(RuntimeError):
    """Raised when the BRENDA bulk JSON cannot be located or validated."""


_INSTRUCTIONS = (
    "Download brenda_<release>.json.tar.gz from "
    "https://www.brenda-enzymes.org/download.php "
    "(accept the CC BY 4.0 licence) and place it in {cache_dir}, "
    "or drop the unpacked .json there directly."
)


def ensure_brenda_json(
    cache_dir: str | Path,
    *,
    expected_sha256: str | None = None,
) -> Path:
    """Return the path to the unpacked BRENDA JSON, extracting if needed.

    Look-up order in ``cache_dir``:
        1. ``*.json`` already unpacked: return as-is.
        2. ``*.json.tar.gz`` tarball: validate ``expected_sha256``
           (if given), extract the single JSON member, return its path.
        3. Neither present: raise ``BrendaDownloadError`` with
           instructions for manual download.

    Parameters
    ----------
    cache_dir
        Directory in which the user drops the BRENDA download.
        Created if missing.
    expected_sha256
        If given, the tarball must hash to this value or the function
        raises ``BrendaDownloadError``. Ignored when only an unpacked
        ``.json`` is present.

    Returns
    -------
    Path
        Absolute path to the unpacked BRENDA JSON.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    unpacked = sorted(cache_dir.glob("*.json"))
    if unpacked:
        return unpacked[0].resolve()

    tarballs = sorted(cache_dir.glob("*.json.tar.gz"))
    if not tarballs:
        raise BrendaDownloadError(_INSTRUCTIONS.format(cache_dir=cache_dir))

    tarball = tarballs[0]
    if expected_sha256 is not None:
        actual = sha256_of(tarball)
        if actual != expected_sha256:
            raise BrendaDownloadError(
                f"sha256 mismatch for {tarball.name}: "
                f"expected {expected_sha256}, got {actual}"
            )

    json_path = extract_brenda_json(tarball, cache_dir)
    return json_path.resolve()


def extract_brenda_json(tarball_path: str | Path, dest_dir: str | Path) -> Path:
    """Extract the single JSON member from a BRENDA tarball.

    Parameters
    ----------
    tarball_path
        Path to ``brenda_<release>.json.tar.gz``.
    dest_dir
        Directory to extract the JSON into. The file is written under
        its name as stored inside the tarball.

    Returns
    -------
    Path
        Path to the extracted JSON file.
    """
    tarball_path = Path(tarball_path)
    dest_dir = Path(dest_dir)
    with tarfile.open(tarball_path, "r:gz") as tar:
        json_members = [m for m in tar.getmembers() if m.name.endswith(".json")]
        if not json_members:
            raise BrendaDownloadError(
                f"no .json member found in {tarball_path.name}"
            )
        if len(json_members) > 1:
            raise BrendaDownloadError(
                f"expected exactly one .json member in {tarball_path.name}, "
                f"found {len(json_members)}: {[m.name for m in json_members]}"
            )
        member = json_members[0]
        # Flatten any directory prefix; the extracted file lives directly
        # in ``dest_dir`` under its basename.
        member.name = Path(member.name).name
        tar.extract(member, dest_dir, filter="data")
    extracted = dest_dir / member.name
    logger.info("extracted %s to %s", tarball_path.name, extracted)
    return extracted


def sha256_of(path: str | Path, *, chunk_size: int = 65536) -> str:
    """Return the hex sha256 digest of a file's contents."""
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()
