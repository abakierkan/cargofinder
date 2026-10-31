"""Unit tests for profile metadata parsing and profile discovery."""

import gzip
import json
from unittest.mock import MagicMock, patch

import pytest

from cargofinder.config import DEFAULT_SHELL_PROFILES
from cargofinder.profiles import (
    _extract_hmm_header_metadata,
    fetch_default_profiles,
    fetch_pfam_profile,
    find_profile_files,
)

SAMPLE_HMM_HEADER = """HMMER3/f [3.3 | Nov 2019]
NAME  Linocin_M18
ACC   PF04454.19
DESC  Encapsulating protein for peroxidase
LENG  253
DATE  Wed May 28 07:46:41 2025
HMM          A        C        D
"""


def test_extract_hmm_header_metadata():
    meta = _extract_hmm_header_metadata(SAMPLE_HMM_HEADER)
    assert meta["name"] == "Linocin_M18"
    assert meta["acc"] == "PF04454.19"
    assert meta["desc"] == "Encapsulating protein for peroxidase"
    assert meta["leng"] == "253"
    assert "2025" in meta["date"]


def test_find_profile_files(tmp_path):
    with pytest.raises(FileNotFoundError, match="Profile directory not found"):
        find_profile_files(tmp_path / "nonexistent")

    with pytest.raises(FileNotFoundError, match="No .hmm or .hmm.gz files found"):
        find_profile_files(tmp_path)

    (tmp_path / "PF04454.hmm").write_text(SAMPLE_HMM_HEADER)
    (tmp_path / "other.hmm.gz").write_bytes(b"dummy")
    (tmp_path / "ignore.txt").write_text("ignore")

    found = find_profile_files(tmp_path)
    assert len(found) == 2
    assert any(f.name == "PF04454.hmm" for f in found)
    assert any(f.name == "other.hmm.gz" for f in found)


def test_fetch_pfam_profile_mocked(tmp_path):
    gzipped_bytes = gzip.compress(SAMPLE_HMM_HEADER.encode("utf-8"))

    mock_resp = MagicMock()
    mock_resp.read.return_value = gzipped_bytes
    mock_resp.__enter__.return_value = mock_resp
    mock_resp.__exit__.return_value = None

    with patch("urllib.request.urlopen", return_value=mock_resp):
        out_file = fetch_pfam_profile("PF04454", tmp_path)
        assert out_file.exists()
        assert out_file.read_text().startswith("HMMER3")

        meta_file = tmp_path / "PF04454.meta.json"
        assert meta_file.exists()
        meta = json.loads(meta_file.read_text())
        assert meta["acc"] == "PF04454.19"
        assert meta["name"] == "Linocin_M18"


def test_fetch_default_profiles_mocked(tmp_path):
    gzipped_bytes = gzip.compress(SAMPLE_HMM_HEADER.encode("utf-8"))

    mock_resp = MagicMock()
    mock_resp.read.return_value = gzipped_bytes
    mock_resp.__enter__.return_value = mock_resp
    mock_resp.__exit__.return_value = None

    with patch("urllib.request.urlopen", return_value=mock_resp):
        files = fetch_default_profiles(tmp_path, DEFAULT_SHELL_PROFILES)
        assert len(files) == 1
        assert files[0].name == "PF04454.hmm"

        manifest_file = tmp_path / "manifest.json"
        assert manifest_file.exists()
        manifest = json.loads(manifest_file.read_text())
        assert "PF04454" in manifest
        assert manifest["PF04454"]["family_label"] == "Family 1 (Linocin_M18)"
