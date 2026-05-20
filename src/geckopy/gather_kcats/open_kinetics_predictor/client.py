"""HTTP client for the OpenKineticsPredictor (OKP) REST API.

https://predictor.openkinetics.org/api/v1 -- async job model:
submit -> poll status -> download result. Auth is a Bearer key.
"""
from __future__ import annotations

import json
import logging
from typing import Optional

import requests

logger = logging.getLogger(__name__)

_BASE_URL = "https://predictor.openkinetics.org/api/v1"
_TIMEOUT = 60
_RESULT_TIMEOUT = 120


class OKPError(RuntimeError):
    """Raised when the OpenKineticsPredictor API returns an error."""


def _bool_str(value: bool) -> str:
    return "true" if value else "false"


class OKPClient:
    """Thin wrapper over the OKP REST endpoints.

    The ``session`` argument exists so tests can inject a mocked
    ``requests.Session``; in production a fresh session is created.
    """

    def __init__(
        self,
        api_key: str,
        *,
        base_url: str = _BASE_URL,
        session: Optional[requests.Session] = None,
    ) -> None:
        self.api_key = api_key
        self.base_url = base_url.rstrip("/")
        self.session = session or requests.Session()

    @property
    def _auth(self) -> dict[str, str]:
        return {"Authorization": f"Bearer {self.api_key}"}

    def submit(
        self,
        csv_bytes: bytes,
        *,
        method: str,
        handle_long_sequences: str = "truncate",
        include_similarity_columns: bool = True,
        canonicalize_substrates: bool = True,
    ) -> dict:
        """POST the input CSV; return the parsed JSON (includes ``jobId``)."""
        resp = self.session.post(
            f"{self.base_url}/submit/",
            headers=self._auth,
            files={"file": ("OKP.csv", csv_bytes, "text/csv")},
            data={
                "targets": '["kcat"]',
                "methods": json.dumps({"kcat": method}),
                "handleLongSequences": handle_long_sequences,
                "includeSimilarityColumns": _bool_str(include_similarity_columns),
                "canonicalizeSubstrates": _bool_str(canonicalize_substrates),
            },
            timeout=_TIMEOUT,
        )
        return self._json_or_raise(resp, "submit")

    def status(self, job_id: str) -> dict:
        """GET the job status JSON."""
        resp = self.session.get(
            f"{self.base_url}/status/{job_id}/",
            headers=self._auth,
            timeout=_TIMEOUT,
        )
        return self._json_or_raise(resp, "status")

    def result_csv(self, job_id: str) -> str:
        """GET the result as CSV text (the JSON form contains bare NaN)."""
        resp = self.session.get(
            f"{self.base_url}/result/{job_id}/",
            headers=self._auth,
            timeout=_RESULT_TIMEOUT,
        )
        if not resp.ok:
            self._raise(resp, "result")
        return resp.text

    def methods(self) -> dict:
        """GET the available-methods listing (no auth required)."""
        resp = self.session.get(f"{self.base_url}/methods/", timeout=_TIMEOUT)
        return self._json_or_raise(resp, "methods")

    # ----------------------------------------------------------------- #
    def _json_or_raise(self, resp: requests.Response, what: str) -> dict:
        if not resp.ok:
            self._raise(resp, what)
        return resp.json()

    def _raise(self, resp: requests.Response, what: str) -> None:
        try:
            message = resp.json().get("error", resp.text)
        except (ValueError, AttributeError):
            message = resp.text
        raise OKPError(
            f"OKP {what} failed (HTTP {resp.status_code}): {message}"
        )
