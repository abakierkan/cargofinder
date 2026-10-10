"""HMM profile management and fetching from public repositories (InterPro/Pfam)."""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

from cargofinder.config import DEFAULT_SHELL_PROFILES, ShellProfileDef

logger = logging.getLogger(__name__)

INTERPRO_HMM_URL_TEMPLATE = (
    "https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}?annotation=hmm"
)


def _extract_hmm_header_metadata(hmm_text: str) -> Dict[str, str]:
    """Parse key header fields from an HMMER3 profile text."""
    meta: Dict[str, str] = {}
    for line in hmm_text.splitlines():
        if line.startswith("HMMER3"):
            meta["format"] = line.strip()
        elif line.startswith("NAME"):
            meta["name"] = line.split(maxsplit=1)[1].strip()
        elif line.startswith("ACC"):
            meta["acc"] = line.split(maxsplit=1)[1].strip()
        elif line.startswith("DESC"):
            meta["desc"] = line.split(maxsplit=1)[1].strip()
        elif line.startswith("LENG"):
            meta["leng"] = line.split(maxsplit=1)[1].strip()
        elif line.startswith("DATE"):
            meta["date"] = line.split(maxsplit=1)[1].strip()
        elif line.startswith("HMM "):
            # Reached model transitions table
            break
    return meta


def fetch_pfam_profile(
    pfam_id: str,
    out_dir: Path | str,
    timeout: int = 30,
) -> Path:
    """Download a single Pfam HMM profile from InterPro and decompress it.

    Parameters
    ----------
    pfam_id : str
        Pfam accession (e.g. 'PF04454').
    out_dir : Path | str
        Destination directory.
    timeout : int
        Network request timeout in seconds.

    Returns
    -------
    Path
        Path to the saved `.hmm` file.
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    target_file = out_path / f"{pfam_id}.hmm"

    url = INTERPRO_HMM_URL_TEMPLATE.format(pfam_id=pfam_id)
    logger.info("Downloading %s HMM from %s", pfam_id, url)

    req = urllib.request.Request(
        url,
        headers={"User-Agent": "cargofinder/0.1.0 (academic research)"},
    )

    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            raw_bytes = response.read()
    except urllib.error.URLError as exc:
        raise RuntimeError(
            f"Failed to download Pfam profile {pfam_id} from InterPro ({url}): {exc}"
        ) from exc

    # InterPro returns application/gzip
    if raw_bytes.startswith(b"\x1f\x8b"):
        try:
            decompressed_bytes = gzip.decompress(raw_bytes)
        except Exception as exc:
            raise ValueError(
                f"Failed to decompress gzip stream for {pfam_id}: {exc}"
            ) from exc
    else:
        decompressed_bytes = raw_bytes

    text = decompressed_bytes.decode("utf-8", errors="replace")
    if not text.startswith("HMMER3"):
        raise ValueError(
            f"Downloaded payload for {pfam_id} is not a valid HMMER3 profile. "
            f"Header: {text[:80]!r}"
        )

    target_file.write_text(text, encoding="utf-8")

    # Record checksum and header metadata
    sha256 = hashlib.sha256(decompressed_bytes).hexdigest()
    header_meta = _extract_hmm_header_metadata(text)
    header_meta["sha256"] = sha256
    header_meta["fetched_at"] = datetime.now(timezone.utc).isoformat()
    header_meta["source_url"] = url

    meta_file = out_path / f"{pfam_id}.meta.json"
    meta_file.write_text(json.dumps(header_meta, indent=2), encoding="utf-8")

    logger.info(
        "Saved %s to %s (Acc: %s)", pfam_id, target_file, header_meta.get("acc")
    )
    return target_file


def fetch_default_profiles(
    out_dir: Path | str,
    profiles: Optional[List[ShellProfileDef]] = None,
    timeout: int = 30,
) -> List[Path]:
    """Fetch default encapsulin shell profiles and write a consolidated manifest.

    Parameters
    ----------
    out_dir : Path | str
        Destination directory.
    profiles : list of ShellProfileDef, optional
        Profile definitions to fetch. Defaults to DEFAULT_SHELL_PROFILES (PF04454).
    timeout : int
        Request timeout.

    Returns
    -------
    list of Path
        Paths to downloaded `.hmm` files.
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    if profiles is None:
        profiles = DEFAULT_SHELL_PROFILES

    downloaded: List[Path] = []
    manifest: Dict[str, dict] = {}

    for prof in profiles:
        hmm_file = fetch_pfam_profile(prof.pfam_id, out_path, timeout=timeout)
        downloaded.append(hmm_file)

        meta_file = out_path / f"{prof.pfam_id}.meta.json"
        if meta_file.exists():
            meta = json.loads(meta_file.read_text(encoding="utf-8"))
        else:
            meta = {}
        meta["family_label"] = prof.family
        meta["family_description"] = prof.description
        manifest[prof.pfam_id] = meta

    manifest_file = out_path / "manifest.json"
    manifest_file.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    return downloaded


def find_profile_files(profile_dir: Path | str) -> List[Path]:
    """Locate all .hmm profiles in a directory.

    Parameters
    ----------
    profile_dir : Path | str
        Path to directory containing profile HMM files.

    Returns
    -------
    list of Path
        List of matching profile file paths sorted alphabetically.
    """
    p = Path(profile_dir)
    if not p.is_dir():
        raise FileNotFoundError(f"Profile directory not found: {profile_dir}")
    files = sorted(list(p.glob("*.hmm")) + list(p.glob("*.hmm.gz")))
    if not files:
        raise FileNotFoundError(
            f"No .hmm or .hmm.gz files found in profile directory: {profile_dir}"
        )
    return files
