"""Invoke DLKcat (via Docker) to predict kcat values for a model.

Ported from GECKO MATLAB:
src/geckomat/gather_kcats/runDLKcat.m.
"""
from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


_DEFAULT_DOCKER_IMAGE = "ghcr.io/sysbiochalmers/dlkcat-gecko:0.1"
_INPUT_FILENAME = "DLKcat.tsv"
_OUTPUT_FILENAME = "DLKcatOutput.tsv"


def run_dlkcat(
    input_path: str | Path,
    *,
    output_path: Optional[str | Path] = None,
    docker_image: str = _DEFAULT_DOCKER_IMAGE,
    timeout: Optional[float] = None,
) -> Path:
    """Run DLKcat on a prepared input TSV and return the output path.

    Ported from GECKO MATLAB:
    src/geckomat/gather_kcats/runDLKcat.m.

    Requires Docker to be installed and the daemon to be running.
    The DLKcat invocation uses a published image and may take several
    minutes (longer on first run when the image is downloaded).

    Internally the function copies ``input_path`` into a temporary
    directory, mounts that directory into the Docker container at
    ``/data``, runs DLKcat, and copies the resulting predictions back
    to ``output_path``. The temporary directory is cleaned up on exit
    regardless of success or failure.

    MATLAB-COMPAT: GECKO MATLAB takes a ``modelAdapter`` and resolves
    the path as ``adapter.params.path/data/DLKcat.tsv``. geckopy takes
    the file path explicitly; the caller resolves it.

    MATLAB-COMPAT: GECKO MATLAB writes temporary files into the
    project's ``data/`` folder. geckopy uses a real temp directory
    (auto-cleaned by ``tempfile.TemporaryDirectory``), so a failed
    run cannot leave stray files in the project.

    MATLAB-COMPAT: The macOS PATH workaround
    (``setenv('PATH', '/usr/local/bin:' + PATH)``) needed when MATLAB
    is launched from Finder is omitted; Python users invoking from a
    shell already have a sensible PATH.

    Parameters
    ----------
    input_path
        Path to the DLKcat input TSV (typically the output of
        ``write_dlkcat_input``).
    output_path
        Where to write DLKcat's predictions. If ``None`` (default),
        ``input_path`` is overwritten, matching MATLAB.
    docker_image
        Image tag to invoke. Defaults to the published GECKO build.
    timeout
        Optional ``subprocess.run`` timeout in seconds for the
        Docker invocation. ``None`` (default) disables it.

    Returns
    -------
    pathlib.Path
        The path to the output file.

    Raises
    ------
    FileNotFoundError
        If ``input_path`` does not exist.
    RuntimeError
        If Docker is not installed, the DLKcat container exits
        non-zero, or no output file is produced.
    """
    input_path = Path(input_path)
    if not input_path.is_file():
        raise FileNotFoundError(f"DLKcat input file not found: {input_path}")

    out_path = Path(output_path) if output_path is not None else input_path

    _check_docker_available()

    with tempfile.TemporaryDirectory(prefix="geckopy_dlkcat_") as tmp_str:
        tmp_dir = Path(tmp_str)
        tmp_input = tmp_dir / _INPUT_FILENAME
        tmp_output = tmp_dir / _OUTPUT_FILENAME
        shutil.copyfile(input_path, tmp_input)

        cmd = [
            "docker", "run", "--rm",
            "-v", f"{tmp_dir}:/data",
            docker_image,
            "/bin/bash", "-c",
            f"python DLKcat.py /data/{_INPUT_FILENAME} "
            f"/data/{_OUTPUT_FILENAME}",
        ]

        logger.info(
            "run_dlkcat: invoking %s (this may take several minutes, "
            "longer on first run)", docker_image,
        )
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"DLKcat Docker run failed (exit code {result.returncode}). "
                f"stderr:\n{result.stderr.strip()}"
            )

        if not tmp_output.is_file():
            raise RuntimeError(
                f"DLKcat Docker run completed but produced no output file "
                f"at /data/{_OUTPUT_FILENAME}. stderr:\n"
                f"{result.stderr.strip()}"
            )

        out_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(tmp_output, out_path)

    logger.info("run_dlkcat: completed; output written to %s", out_path)
    return out_path


def _check_docker_available() -> None:
    """Raise RuntimeError if `docker --version` is not invokable."""
    try:
        result = subprocess.run(
            ["docker", "--version"],
            capture_output=True, text=True, timeout=10,
        )
    except FileNotFoundError as e:
        raise RuntimeError(
            "Cannot find Docker. Install Docker (https://www.docker.com) "
            "and ensure the daemon is running before invoking run_dlkcat."
        ) from e
    if result.returncode != 0:
        raise RuntimeError(
            f"Docker is installed but `docker --version` failed with exit "
            f"code {result.returncode}. stderr:\n{result.stderr.strip()}"
        )
