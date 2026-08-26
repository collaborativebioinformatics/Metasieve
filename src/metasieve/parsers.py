"""Native Python parsing: FASTA I/O, Kraken2 filtering, SeqScreen ORF extraction.

These functions replace the Nextflow ``bin/*.py`` helpers. FASTA read/write uses
Biopython; ORF translation uses NCBI table 11 (bacterial / plastid).
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator, Optional, Sequence

from Bio import SeqIO
from Bio.Seq import Seq

LOGGER = logging.getLogger(__name__)

STOP = "*"
START_CODONS = {"ATG", "GTG", "TTG"}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
AA_FASTA_SUFFIXES = {".faa", ".fa", ".fasta", ".fna", ".fas"}
PROTEIN_NAME_HINTS = (
    "translated",
    "translation",
    "six_frame",
    "sixframe",
    "orf",
    "protein",
    "aa.fasta",
    "aa.fa",
    ".faa",
)
STRUCTURE_SUFFIXES = {".pdb", ".cif", ".mmcif", ".ent"}


@dataclass
class FastaRecord:
    record_id: str
    header: str
    sequence: str


@dataclass
class UnclassifiedResult:
    fasta: Path
    ids: Path
    manifest: Path
    stats: Path
    n_written: int


@dataclass
class OrfResult:
    fasta: Path
    manifest: Path
    split_dir: Path
    stats: Path
    n_orfs: int
    split_fastas: list[Path] = field(default_factory=list)


def read_fasta(path: Path) -> Iterator[FastaRecord]:
    """Yield FASTA records. ``header`` is the full description line (no ``>``)."""
    path = Path(path)
    for record in SeqIO.parse(str(path), "fasta"):
        yield FastaRecord(
            record_id=record.id,
            header=record.description,
            sequence=str(record.seq).replace(" ", "").replace("\n", ""),
        )


def wrap_fasta(seq: str, width: int = 80) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def write_fasta_record(handle, header: str, seq: str, width: int = 80) -> None:
    handle.write(f">{header}\n")
    handle.write("\n".join(wrap_fasta(seq, width)) + "\n")


def filter_fasta_min_length(
    source: Path,
    dest: Path,
    min_length: int = 0,
) -> int:
    """Copy FASTA records, dropping sequences shorter than *min_length*.

    Returns the number of records written.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    n_too_short = 0
    with dest.open("w", encoding="utf-8") as handle:
        for rec in read_fasta(source):
            if min_length and len(rec.sequence) < min_length:
                n_too_short += 1
                continue
            write_fasta_record(handle, rec.header, rec.sequence)
            n_written += 1
    LOGGER.info(
        "Unclassified FASTA %s: kept=%d dropped_short=%d min_length=%d",
        source,
        n_written,
        n_too_short,
        min_length,
    )
    return n_written


def sanitize_id(text: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    cleaned = re.sub(r"_+", "_", cleaned).strip("._")
    return cleaned or "unnamed"


def contig_id_from_header(header: str) -> str:
    rec_id = header.split()[0]
    if "|" in rec_id:
        return rec_id.split("|", 1)[1]
    return rec_id


def sha256_short(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii", errors="ignore")).hexdigest()[:12]


# ---------------------------------------------------------------------------
# Kraken2
# ---------------------------------------------------------------------------


def parse_kraken_assignments(path: Path) -> dict[str, dict[str, Any]]:
    """Parse Kraken2 per-read output.

    Columns: classified_flag, sequence_id, taxid, length, lca_mapping
    Unclassified when flag is ``U`` or taxid is ``0``.
    """
    assignments: dict[str, dict[str, Any]] = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                LOGGER.warning("Skipping malformed Kraken2 line %d in %s", lineno, path)
                continue
            flag = parts[0].strip()
            seq_id = parts[1].strip()
            taxid = parts[2].strip()
            length = parts[3].strip() if len(parts) > 3 else ""
            lca = parts[4] if len(parts) > 4 else ""
            unclassified = flag.upper().startswith("U") or taxid in {"0", "unclassified"}
            assignments[seq_id] = {
                "classified_flag": flag,
                "taxid": taxid,
                "length": length,
                "lca_mapping": lca,
                "unclassified": unclassified,
            }
    return assignments


def parse_kraken_report_unclassified(path: Path) -> dict[str, Any]:
    """Pull the unclassified (taxid 0) row from a standard Kraken2 report."""
    info: dict[str, Any] = {
        "report_unclassified_pct": None,
        "report_unclassified_reads": None,
    }
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            parts = raw.rstrip("\n\r").split("\t")
            if len(parts) < 6:
                continue
            taxid = parts[4].strip()
            name = parts[5].strip().lower()
            if taxid == "0" or name == "unclassified":
                try:
                    info["report_unclassified_pct"] = float(parts[0].strip())
                except ValueError:
                    pass
                try:
                    info["report_unclassified_reads"] = int(parts[1].strip())
                except ValueError:
                    pass
                break
    return info


def extract_unclassified_contigs(
    *,
    contigs: Path,
    kraken_out: Path,
    sample_id: str,
    output_fasta: Path,
    output_ids: Path,
    output_manifest: Path,
    output_stats: Path,
    kraken_report: Path | None = None,
    min_length: int = 0,
) -> UnclassifiedResult:
    """Keep Kraken2-unclassified contigs and write FASTA + tracking files."""
    contigs = Path(contigs)
    kraken_out = Path(kraken_out)
    if not contigs.is_file():
        raise FileNotFoundError(f"Contig FASTA not found: {contigs}")
    if not kraken_out.is_file():
        raise FileNotFoundError(f"Kraken2 assignment file not found: {kraken_out}")

    assignments = parse_kraken_assignments(kraken_out)
    unclassified_ids = {seq_id for seq_id, rec in assignments.items() if rec["unclassified"]}
    LOGGER.info(
        "Kraken2 assignments=%d unclassified=%d",
        len(assignments),
        len(unclassified_ids),
    )

    output_fasta = Path(output_fasta)
    output_ids = Path(output_ids)
    output_manifest = Path(output_manifest)
    output_stats = Path(output_stats)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_ids.parent.mkdir(parents=True, exist_ok=True)
    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    output_stats.parent.mkdir(parents=True, exist_ok=True)

    n_fasta = 0
    n_written = 0
    n_too_short = 0
    n_missing_assignment = 0
    manifest_rows: list[dict[str, Any]] = []

    with output_fasta.open("w", encoding="utf-8") as fasta_out:
        for rec in read_fasta(contigs):
            n_fasta += 1
            kraken_rec = assignments.get(rec.record_id)
            if kraken_rec is None:
                n_missing_assignment += 1
                LOGGER.debug("No Kraken2 record for contig %s", rec.record_id)
                continue
            if not kraken_rec["unclassified"]:
                continue
            if min_length and len(rec.sequence) < min_length:
                n_too_short += 1
                continue

            new_header = (
                f"{sample_id}|{rec.record_id} "
                f"sample={sample_id} contig_id={rec.record_id} "
                f"kraken_flag={kraken_rec['classified_flag']} "
                f"kraken_taxid={kraken_rec['taxid']} "
                f"length={len(rec.sequence)} orig_header={rec.header}"
            )
            write_fasta_record(fasta_out, new_header, rec.sequence)
            n_written += 1
            manifest_rows.append(
                {
                    "sample_id": sample_id,
                    "contig_id": rec.record_id,
                    "fasta_record_id": f"{sample_id}|{rec.record_id}",
                    "kraken_flag": kraken_rec["classified_flag"],
                    "kraken_taxid": kraken_rec["taxid"],
                    "nt_length": len(rec.sequence),
                    "orig_header": rec.header,
                }
            )

    with output_ids.open("w", encoding="utf-8") as handle:
        for row in manifest_rows:
            handle.write(row["contig_id"] + "\n")

    with output_manifest.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "sample_id",
            "contig_id",
            "fasta_record_id",
            "kraken_flag",
            "kraken_taxid",
            "nt_length",
            "orig_header",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)

    stats: dict[str, Any] = {
        "sample_id": sample_id,
        "contigs_in_fasta": n_fasta,
        "kraken_assignments": len(assignments),
        "kraken_unclassified_ids": len(unclassified_ids),
        "unclassified_written": n_written,
        "unclassified_too_short": n_too_short,
        "fasta_records_without_kraken": n_missing_assignment,
        "min_length": min_length,
    }
    if kraken_report and Path(kraken_report).is_file():
        stats.update(parse_kraken_report_unclassified(Path(kraken_report)))

    with output_stats.open("w", encoding="utf-8") as handle:
        json.dump(stats, handle, indent=2)
        handle.write("\n")

    LOGGER.info("Wrote %d unclassified contigs -> %s", n_written, output_fasta)
    if n_written == 0:
        LOGGER.warning("No unclassified contigs passed filters.")

    return UnclassifiedResult(
        fasta=output_fasta,
        ids=output_ids,
        manifest=output_manifest,
        stats=output_stats,
        n_written=n_written,
    )


# ---------------------------------------------------------------------------
# SeqScreen / ORF extraction
# ---------------------------------------------------------------------------


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate_frame(nt: str, frame: int, table: int = 11) -> str:
    """Translate one forward frame (0, 1, 2) with the given NCBI table."""
    seq = Seq(nt[frame:])
    trim = len(seq) - (len(seq) % 3)
    if trim <= 0:
        return ""
    return str(seq[:trim].translate(table=table))


def _maybe_add_orf(
    hits: list[dict[str, Any]],
    peptide: str,
    start_aa: int,
    end_aa: int,
    frame: int,
    strand: str,
    nt_len: int,
    min_aa: int,
    max_aa: int,
) -> None:
    peptide = peptide.replace(STOP, "")
    if len(peptide) < min_aa or len(peptide) > max_aa:
        return
    if strand == "+":
        nt_start = frame + start_aa * 3 + 1
        nt_end = frame + end_aa * 3
    else:
        nt_end = nt_len - frame - start_aa * 3
        nt_start = nt_len - frame - end_aa * 3 + 1
        nt_start = max(1, nt_start)
        nt_end = min(nt_len, nt_end)
    hits.append(
        {
            "aa_seq": peptide,
            "nt_start": int(nt_start),
            "nt_end": int(nt_end),
            "strand": strand,
            "frame": frame + 1 if strand == "+" else -(frame + 1),
            "aa_length": len(peptide),
            "source": "six_frame",
        }
    )


def orfs_from_translation(
    aa: str,
    frame: int,
    strand: str,
    nt_len: int,
    require_start: bool,
    min_aa: int,
    max_aa: int,
) -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    start_aa: int | None = None
    for idx, residue in enumerate(aa):
        if start_aa is None:
            if require_start and residue != "M":
                continue
            start_aa = idx
        if residue == STOP:
            if start_aa is not None:
                peptide = aa[start_aa:idx]
                _maybe_add_orf(
                    hits, peptide, start_aa, idx, frame, strand, nt_len, min_aa, max_aa
                )
            start_aa = None
    if start_aa is not None:
        peptide = aa[start_aa:]
        _maybe_add_orf(
            hits, peptide, start_aa, len(aa), frame, strand, nt_len, min_aa, max_aa
        )
    return hits


def find_orfs_six_frame(
    nt: str,
    min_aa: int,
    max_aa: int,
    require_start: bool,
    table: int = 11,
) -> list[dict[str, Any]]:
    seq = re.sub(r"[^ACGTNacgtn]", "N", nt).upper()
    rc = reverse_complement(seq)
    orfs: list[dict[str, Any]] = []
    for frame in range(3):
        orfs.extend(
            orfs_from_translation(
                translate_frame(seq, frame, table=table),
                frame,
                "+",
                len(seq),
                require_start,
                min_aa,
                max_aa,
            )
        )
        orfs.extend(
            orfs_from_translation(
                translate_frame(rc, frame, table=table),
                frame,
                "-",
                len(seq),
                require_start,
                min_aa,
                max_aa,
            )
        )
    if require_start:
        filtered = []
        for orf in orfs:
            if orf["strand"] == "+":
                codon = seq[orf["nt_start"] - 1 : orf["nt_start"] + 2]
            else:
                codon = reverse_complement(seq[orf["nt_end"] - 3 : orf["nt_end"]])
            if codon in START_CODONS:
                filtered.append(orf)
        return filtered
    return orfs


def find_seqscreen_report(seqscreen_dir: Path) -> Optional[Path]:
    seqscreen_dir = Path(seqscreen_dir)
    candidates = [
        seqscreen_dir / "report_generation" / "seqscreen_report.tsv",
        seqscreen_dir / "seqscreen_report.tsv",
    ]
    for path in candidates:
        if path.is_file():
            return path
    matches = sorted(seqscreen_dir.rglob("*seqscreen_report*.tsv"))
    return matches[0] if matches else None


def load_seqscreen_tsv(path: Path) -> dict[str, dict[str, Any]]:
    """Index SeqScreen rows by query / contig identifier."""
    by_query: dict[str, dict[str, Any]] = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return by_query
        query_keys = [
            key
            for key in reader.fieldnames
            if key
            and key.lower()
            in {
                "query",
                "query_id",
                "seqid",
                "seq_id",
                "sequence_id",
                "read",
                "read_id",
                "contig",
                "contig_id",
            }
        ]
        if not query_keys:
            query_keys = [reader.fieldnames[0]]
        for row in reader:
            raw = None
            for key in query_keys:
                if row.get(key):
                    raw = row[key].strip()
                    break
            if not raw:
                continue
            query_id = raw.split()[0]
            contig_id = query_id.split("|")[-1]
            annotation = {k: (v if v is not None else "") for k, v in row.items() if k}
            annotation["_query_id"] = query_id
            annotation["_contig_id"] = contig_id
            by_query[query_id] = annotation
            by_query[contig_id] = annotation
    return by_query


def annotation_blob(ann: Optional[dict[str, Any]]) -> str:
    if not ann:
        return ""
    skip = {"_query_id", "_contig_id"}
    parts = []
    for key, value in ann.items():
        if key in skip or value in (None, ""):
            continue
        if key.lower() in {"protein_sequence", "translation", "aa_seq", "sequence"}:
            continue
        parts.append(f"{key}={str(value).replace(chr(9), ' ')}")
    return ";".join(parts)


def seqscreen_protein_from_row(ann: Optional[dict[str, Any]]) -> Optional[str]:
    if not ann:
        return None
    wanted = {"protein_sequence", "translation", "aa_seq", "predicted_protein"}
    for candidate_key, value in ann.items():
        if candidate_key and candidate_key.lower() in wanted and value:
            seq = re.sub(r"[^A-Za-z*]", "", value)
            if len(seq) >= 10:
                return seq.replace("*", "")
    return None


def discover_protein_fastas(seqscreen_dir: Path) -> list[Path]:
    found: list[Path] = []
    seqscreen_dir = Path(seqscreen_dir)
    if not seqscreen_dir.is_dir():
        return found
    for path in seqscreen_dir.rglob("*"):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        name = path.name.lower()
        if suffix not in AA_FASTA_SUFFIXES and not name.endswith(".faa"):
            continue
        if any(hint in name or hint in str(path).lower() for hint in PROTEIN_NAME_HINTS):
            found.append(path)
    return sorted(set(found))


def index_seqscreen_protein_fastas(
    paths: Sequence[Path],
) -> dict[str, list[tuple[str, str]]]:
    index: dict[str, list[tuple[str, str]]] = {}
    for path in paths:
        try:
            for rec in read_fasta(path):
                peptide = re.sub(r"[^A-Za-z*]", "", rec.sequence)
                if not peptide:
                    continue
                keys = {rec.record_id, rec.record_id.split("|")[-1]}
                for key in keys:
                    index.setdefault(key, []).append((rec.header, peptide))
        except OSError as exc:
            LOGGER.warning("Could not read %s: %s", path, exc)
    return index


def lookup_annotation(
    seqscreen_ann: dict[str, dict[str, Any]], rec_id: str, contig_id: str
) -> Optional[dict[str, Any]]:
    for key in (rec_id, contig_id, rec_id.split("|")[-1]):
        if key in seqscreen_ann:
            return seqscreen_ann[key]
    return None


def pick_orfs_for_contig(
    contig_id: str,
    rec_id: str,
    nt_seq: str,
    seqscreen_ann: dict[str, dict[str, Any]],
    seqscreen_proteins: dict[str, list[tuple[str, str]]],
    min_aa: int,
    max_aa: int,
    require_start: bool,
    table: int = 11,
) -> list[dict[str, Any]]:
    orfs: list[dict[str, Any]] = []
    for key in (rec_id, contig_id, rec_id.split("|")[-1]):
        for header, aa in seqscreen_proteins.get(key, []):
            peptide = aa.replace("*", "")
            if min_aa <= len(peptide) <= max_aa:
                orfs.append(
                    {
                        "aa_seq": peptide,
                        "nt_start": "",
                        "nt_end": "",
                        "strand": "",
                        "frame": "",
                        "aa_length": len(peptide),
                        "source": "seqscreen_fasta",
                        "seqscreen_header": header,
                    }
                )
        ann = seqscreen_ann.get(key)
        prot = seqscreen_protein_from_row(ann)
        if prot and min_aa <= len(prot) <= max_aa:
            orfs.append(
                {
                    "aa_seq": prot,
                    "nt_start": "",
                    "nt_end": "",
                    "strand": "",
                    "frame": "",
                    "aa_length": len(prot),
                    "source": "seqscreen_tsv",
                    "seqscreen_header": "",
                }
            )

    if orfs:
        uniq: dict[str, dict[str, Any]] = {}
        for orf in orfs:
            uniq.setdefault(orf["aa_seq"], orf)
        return list(uniq.values())

    return find_orfs_six_frame(nt_seq, min_aa, max_aa, require_start, table=table)


def extract_seqscreen_orfs(
    *,
    sample_id: str,
    unclassified_fasta: Path,
    seqscreen_dir: Path,
    output_fasta: Path,
    output_manifest: Path,
    output_split_dir: Path,
    output_stats: Path,
    min_aa: int = 50,
    max_aa: int = 1024,
    require_start: bool = True,
    genetic_code: int = 11,
) -> OrfResult:
    """Call ORFs on unclassified contigs and merge SeqScreen annotations.

    Prefers SeqScreen-derived protein records when they match a contig ID,
    otherwise six-frame translates the nucleotide FASTA. FASTA record IDs are
    filesystem-safe: ``{sample}__{contig}__ORF_0001``.
    """
    unclassified_fasta = Path(unclassified_fasta)
    seqscreen_dir = Path(seqscreen_dir)
    if not unclassified_fasta.is_file():
        raise FileNotFoundError(f"Unclassified FASTA not found: {unclassified_fasta}")

    report_path = find_seqscreen_report(seqscreen_dir)
    seqscreen_ann: dict[str, dict[str, Any]] = {}
    if report_path:
        LOGGER.info("SeqScreen report: %s", report_path)
        seqscreen_ann = load_seqscreen_tsv(report_path)
        LOGGER.info("Indexed %d SeqScreen query keys", len(seqscreen_ann))
    else:
        LOGGER.warning("No seqscreen_report.tsv under %s", seqscreen_dir)

    protein_fastas = discover_protein_fastas(seqscreen_dir)
    LOGGER.info("SeqScreen protein FASTA candidates: %d", len(protein_fastas))
    seqscreen_proteins = index_seqscreen_protein_fastas(protein_fastas)

    output_split_dir = Path(output_split_dir)
    output_fasta = Path(output_fasta)
    output_manifest = Path(output_manifest)
    output_stats = Path(output_stats)
    output_split_dir.mkdir(parents=True, exist_ok=True)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_manifest.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "sample_id",
        "contig_id",
        "orf_id",
        "record_id",
        "fasta_header",
        "split_fasta",
        "source",
        "nt_start",
        "nt_end",
        "strand",
        "frame",
        "aa_length",
        "aa_sha256_12",
        "seqscreen_query",
        "seqscreen_annotation",
        "orig_nt_header",
    ]

    n_contigs = 0
    n_orfs = 0
    n_from_seqscreen = 0
    n_from_sixframe = 0
    split_fastas: list[Path] = []

    with output_fasta.open("w", encoding="utf-8") as fasta_out, output_manifest.open(
        "w", encoding="utf-8", newline=""
    ) as csv_out:
        writer = csv.DictWriter(csv_out, fieldnames=fieldnames)
        writer.writeheader()

        for rec in read_fasta(unclassified_fasta):
            n_contigs += 1
            contig_id = contig_id_from_header(rec.header)
            ann = lookup_annotation(seqscreen_ann, rec.record_id, contig_id)
            orfs = pick_orfs_for_contig(
                contig_id,
                rec.record_id,
                rec.sequence,
                seqscreen_ann,
                seqscreen_proteins,
                min_aa,
                max_aa,
                require_start,
                table=genetic_code,
            )
            safe_contig = sanitize_id(contig_id)
            for idx, orf in enumerate(orfs, start=1):
                orf_id = f"ORF_{idx:04d}"
                record_id = f"{sanitize_id(sample_id)}__{safe_contig}__{orf_id}"
                fasta_header = (
                    f"{record_id} sample={sample_id} contig_id={contig_id} "
                    f"orf_id={orf_id} start={orf.get('nt_start', '')} "
                    f"end={orf.get('nt_end', '')} strand={orf.get('strand', '')} "
                    f"frame={orf.get('frame', '')} aa_len={orf['aa_length']} "
                    f"source={orf['source']}"
                )
                write_fasta_record(fasta_out, fasta_header, orf["aa_seq"])

                split_name = f"{record_id}.faa"
                split_path = output_split_dir / split_name
                with split_path.open("w", encoding="utf-8") as split_handle:
                    write_fasta_record(split_handle, fasta_header, orf["aa_seq"])
                split_fastas.append(split_path)

                writer.writerow(
                    {
                        "sample_id": sample_id,
                        "contig_id": contig_id,
                        "orf_id": orf_id,
                        "record_id": record_id,
                        "fasta_header": fasta_header,
                        "split_fasta": split_name,
                        "source": orf["source"],
                        "nt_start": orf.get("nt_start", ""),
                        "nt_end": orf.get("nt_end", ""),
                        "strand": orf.get("strand", ""),
                        "frame": orf.get("frame", ""),
                        "aa_length": orf["aa_length"],
                        "aa_sha256_12": sha256_short(orf["aa_seq"]),
                        "seqscreen_query": (ann or {}).get("_query_id", ""),
                        "seqscreen_annotation": annotation_blob(ann),
                        "orig_nt_header": rec.header,
                    }
                )
                n_orfs += 1
                if str(orf["source"]).startswith("seqscreen"):
                    n_from_seqscreen += 1
                else:
                    n_from_sixframe += 1

    stats = {
        "sample_id": sample_id,
        "contigs": n_contigs,
        "orfs_written": n_orfs,
        "orfs_from_seqscreen": n_from_seqscreen,
        "orfs_from_six_frame": n_from_sixframe,
        "seqscreen_report": str(report_path) if report_path else "",
        "seqscreen_protein_fastas": [str(p) for p in protein_fastas],
        "min_aa": min_aa,
        "max_aa": max_aa,
    }
    with output_stats.open("w", encoding="utf-8") as handle:
        json.dump(stats, handle, indent=2)
        handle.write("\n")

    LOGGER.info("Wrote %d ORFs from %d contigs -> %s", n_orfs, n_contigs, output_fasta)
    if n_orfs == 0:
        LOGGER.warning("No ORFs passed length filters.")

    return OrfResult(
        fasta=output_fasta,
        manifest=output_manifest,
        split_dir=output_split_dir,
        stats=output_stats,
        n_orfs=n_orfs,
        split_fastas=sorted(split_fastas),
    )


# ---------------------------------------------------------------------------
# Structure tracking
# ---------------------------------------------------------------------------


def parse_fasta_header_fields(fasta: Path) -> dict[str, Any]:
    header = ""
    seq_chunks: list[str] = []
    with Path(fasta).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if header:
                    break
                header = line[1:].strip()
            else:
                seq_chunks.append(re.sub(r"\s+", "", line))
    rec_id = header.split()[0] if header else Path(fasta).stem
    fields: dict[str, Any] = {
        "record_id": rec_id,
        "fasta_header": header,
        "aa_length": len("".join(seq_chunks)),
    }
    for match in re.finditer(r"(\w+)=([^\s]+)", header):
        fields[match.group(1)] = match.group(2)
    return fields


def write_fold_sidecar(
    *,
    fasta: Path,
    structure: Path,
    method: str,
    sample_id: str,
    output: Path,
    model_name: str = "",
    plddt: str = "",
) -> Path:
    """Write a JSON sidecar linking a predicted structure to its ORF record."""
    meta = parse_fasta_header_fields(fasta)
    payload = {
        "sample_id": sample_id,
        "record_id": meta.get("record_id", Path(fasta).stem),
        "contig_id": meta.get("contig_id", ""),
        "orf_id": meta.get("orf_id", ""),
        "fasta_header": meta.get("fasta_header", ""),
        "aa_length": meta.get("aa_length", ""),
        "method": method,
        "structure_file": str(structure),
        "model_name": model_name,
        "plddt": plddt,
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    return output


def load_orf_manifest(path: Path) -> dict[str, dict[str, str]]:
    by_id: dict[str, dict[str, str]] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            record_id = (row.get("record_id") or "").strip()
            if record_id:
                by_id[record_id] = row
                by_id[record_id.lower()] = row
    return by_id


def _load_json(path: Path) -> Optional[dict[str, Any]]:
    try:
        with Path(path).open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        LOGGER.warning("Could not parse sidecar %s: %s", path, exc)
        return None


def guess_record_id(name: str, orf_index: dict[str, dict[str, str]]) -> Optional[str]:
    stem = Path(name).stem
    stem = re.sub(r"\.(esmfold|unrelaxed_rank.*|rank_.*)$", "", stem)
    if stem in orf_index:
        return orf_index[stem]["record_id"]
    for record_id, row in orf_index.items():
        if record_id.lower() == stem.lower():
            return row["record_id"]
        if stem.startswith(record_id) or record_id in stem:
            return row["record_id"]
    return None


def stamp_pdb(src: Path, dest: Path, remarks: Iterable[str]) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    body = Path(src).read_text(encoding="utf-8", errors="replace")
    remark_block = "".join(f"REMARK   0 {line}\n" for line in remarks)
    dest.write_text(remark_block + body, encoding="utf-8")


def copy_cif(src: Path, dest: Path, remarks: Iterable[str]) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    comment = "".join(f"# METASIEVE {line}\n" for line in remarks)
    body = Path(src).read_text(encoding="utf-8", errors="replace")
    dest.write_text(comment + body, encoding="utf-8")


def build_structure_manifest(
    *,
    orf_manifest: Path,
    output_manifest: Path,
    output_dir: Path,
    sidecar_paths: Sequence[Path] | None = None,
    structure_paths: Sequence[Path] | None = None,
) -> Path:
    """Rename/stamp PDB/CIF files and write contig → ORF → structure CSV."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    orf_index = load_orf_manifest(Path(orf_manifest))

    sidecar_paths = [Path(p) for p in (sidecar_paths or [])]
    structure_paths = [Path(p) for p in (structure_paths or [])]

    sidecars: list[dict[str, Any]] = []
    for path in sidecar_paths:
        if path.is_dir():
            for child in path.rglob("*.json"):
                payload = _load_json(child)
                if payload:
                    sidecars.append(payload)
        elif path.suffix.lower() == ".json":
            payload = _load_json(path)
            if payload:
                sidecars.append(payload)

    if not sidecars:
        for path in structure_paths:
            if path.suffix.lower() not in STRUCTURE_SUFFIXES:
                continue
            sidecars.append(
                {
                    "structure_file": str(path),
                    "record_id": guess_record_id(path.name, orf_index) or "",
                    "method": "unknown",
                    "sample_id": "",
                }
            )

    fieldnames = [
        "sample_id",
        "contig_id",
        "orf_id",
        "record_id",
        "folding_method",
        "structure_file",
        "published_structure",
        "aa_length",
        "nt_start",
        "nt_end",
        "strand",
        "frame",
        "aa_sha256_12",
        "fasta_header",
        "seqscreen_query",
        "source_orf",
        "model_name",
        "plddt",
    ]

    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    for payload in sidecars:
        src = Path(payload.get("structure_file") or payload.get("pdb") or "")
        if not src.is_file():
            matches = [
                p
                for p in structure_paths
                if p.name == src.name or p.name == Path(str(payload.get("structure_file", ""))).name
            ]
            if matches:
                src = matches[0]
            else:
                LOGGER.warning("Missing structure file for sidecar: %s", payload)
                continue

        record_id = payload.get("record_id") or guess_record_id(src.name, orf_index) or ""
        orf_row = orf_index.get(record_id) or orf_index.get(str(record_id).lower(), {})
        method = (payload.get("method") or payload.get("folding_method") or "unknown").lower()
        sample_id = payload.get("sample_id") or orf_row.get("sample_id", "")
        dest_name = f"{record_id or src.stem}.{method}{src.suffix.lower()}"
        dest = output_dir / method / dest_name
        key = (str(record_id), method, src.name)
        if key in seen:
            continue
        seen.add(key)

        remarks = [
            "METASIEVE TRACKING",
            f"sample_id={sample_id}",
            f"contig_id={orf_row.get('contig_id', '')}",
            f"orf_id={orf_row.get('orf_id', '')}",
            f"record_id={record_id}",
            f"folding_method={method}",
            f"fasta_header={orf_row.get('fasta_header', '')}",
            f"aa_sha256_12={orf_row.get('aa_sha256_12', '')}",
        ]
        if src.suffix.lower() == ".pdb":
            stamp_pdb(src, dest, remarks)
        else:
            copy_cif(src, dest, remarks)

        rows.append(
            {
                "sample_id": sample_id,
                "contig_id": orf_row.get("contig_id", ""),
                "orf_id": orf_row.get("orf_id", ""),
                "record_id": record_id,
                "folding_method": method,
                "structure_file": src.name,
                "published_structure": str(Path(method) / dest_name),
                "aa_length": orf_row.get("aa_length", payload.get("aa_length", "")),
                "nt_start": orf_row.get("nt_start", ""),
                "nt_end": orf_row.get("nt_end", ""),
                "strand": orf_row.get("strand", ""),
                "frame": orf_row.get("frame", ""),
                "aa_sha256_12": orf_row.get("aa_sha256_12", ""),
                "fasta_header": orf_row.get("fasta_header", ""),
                "seqscreen_query": orf_row.get("seqscreen_query", ""),
                "source_orf": orf_row.get("source", ""),
                "model_name": payload.get("model_name", ""),
                "plddt": payload.get("plddt", ""),
            }
        )

    output_manifest = Path(output_manifest)
    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    with output_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    LOGGER.info("Indexed %d structures -> %s", len(rows), output_manifest)
    return output_manifest
