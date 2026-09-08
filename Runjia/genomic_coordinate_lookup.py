# ============================================================
# GENOMIC COORDINATE LOOKUP for E. coli sgRNAID -> target locus sequence
#
# Implements literature-review Section 5, Step 1 ("do first" — enabler
# for families D, E, F, G, H and for Section A of
# proposed_extensions_flanking_folding_readcount.py):
#
#   1. Parse sgRNAID (e.g. "aaeAb3241_107_Cas9") into gene name, b-number,
#      and offset.
#   2. Map the b-number through the MG1655 GenBank annotation to a genomic
#      coordinate and strand.
#   3. Pull the +/-200 nt window around that coordinate from NC_000913.3.
#   4. Sanity-check by reconstructing the protospacer from the genome and
#      confirming it's plausible (exact match against the matrix's own
#      positional encoding needs the p1..p20 column semantics, which
#      aren't confirmed here -- see the CAVEAT below).
#
# STATUS / WHAT IS AND ISN'T VERIFIED:
#   - The sgRNAID parser is built from the ONE example format given in
#     beyond_QCT_literature_review.md ("aaeAb3241_107_Cas9" -> gene aaeA,
#     b-number b3241, offset 107). This script does not have access to
#     ecoli_feature_matrix.csv (it lives on your machine, not in this
#     sandbox), so the regex has not been validated against your real
#     sgRNAID values. Run validate_sgRNAID_parsing() against your actual
#     CSV first -- it reports every ID that fails to parse instead of
#     silently skipping rows.
#   - The offset-to-coordinate convention (does "offset 107" count from
#     the gene's start codon on the SENSE strand, or is it an absolute
#     left-to-right genome offset regardless of strand?) is NOT specified
#     in the literature review and is not guessable from one example.
#     Both conventions are implemented (`offset_mode="sense_relative"`
#     default, or `"absolute"`) -- pick the one that reconstructs a
#     protospacer matching what's already encoded in your matrix, and
#     treat this as something to confirm, not assume.
#   - Needs a local MG1655 GenBank record (NC_000913.3) and matching
#     FASTA. This sandbox has no network access to fetch them; download_
#     reference_genome() will attempt an NCBI Entrez fetch when you run
#     this on a machine that does have internet access, or point
#     GENBANK_PATH at a file you've already downloaded.
# ============================================================

# %% CELL 0 — Install missing packages
import sys
import subprocess


def install_if_missing(import_name, pip_name=None):
    pip_name = pip_name or import_name
    try:
        __import__(import_name)
    except ImportError:
        print(f"Installing {pip_name}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name])


for _mod, _pip in [("Bio", "biopython"), ("pandas", None)]:
    install_if_missing(_mod, _pip)

# %% CELL 1 — Imports and config
import os
import re
import pandas as pd
from Bio import SeqIO

REFERENCE_ACCESSION = "NC_000913.3"  # E. coli K-12 MG1655
GENBANK_PATH = "NC_000913.3.gb"  # local path once downloaded
FLANK_WINDOW_NT = 200  # +/- window pulled around the target coordinate

# sgRNAID format from the literature review's own example:
#   "aaeAb3241_107_Cas9" -> gene="aaeA", bnum="b3241", offset=107, rest="Cas9"
# Non-greedy gene group stops right before the b-number pattern.
SGRNAID_PATTERN = re.compile(
    r"^(?P<gene>[A-Za-z0-9]+?)(?P<bnum>b\d{1,5})_(?P<offset>\d+)_(?P<rest>.*)$"
)


# %% CELL 2 — sgRNAID parsing
def parse_sgRNAID(sgRNAID: str) -> dict:
    """Parse one sgRNAID into {gene, bnum, offset, rest}. Raises ValueError
    (not a silent None) if the ID doesn't match the expected format, since
    a silently-skipped guide would quietly shrink your dataset."""
    m = SGRNAID_PATTERN.match(sgRNAID)
    if not m:
        raise ValueError(
            f"sgRNAID {sgRNAID!r} did not match the expected "
            f"<gene><bnumber>_<offset>_<rest> pattern (e.g. 'aaeAb3241_107_Cas9')."
        )
    return {
        "gene": m.group("gene"),
        "bnum": m.group("bnum"),
        "offset": int(m.group("offset")),
        "rest": m.group("rest"),
    }


def validate_sgRNAID_parsing(sgRNAID_series) -> pd.DataFrame:
    """Run parse_sgRNAID over every value in a pandas Series (e.g.
    df['sgRNAID']) and report failures instead of assuming the format
    holds project-wide. Run this against your real CSV BEFORE trusting
    any coordinate this script produces."""
    records = []
    for sid in sgRNAID_series:
        try:
            parsed = parse_sgRNAID(sid)
            parsed["sgRNAID"] = sid
            parsed["parsed_ok"] = True
        except ValueError:
            parsed = {"sgRNAID": sid, "parsed_ok": False}
        records.append(parsed)
    result = pd.DataFrame(records)
    n_fail = (~result["parsed_ok"]).sum()
    print(f"Parsed {len(result) - n_fail}/{len(result)} sgRNAIDs successfully.")
    if n_fail:
        print(f"{n_fail} failed to parse — inspect result[~result.parsed_ok] "
              "before trusting the format assumption.")
    return result


# %% CELL 3 — Reference genome loading
def download_reference_genome(email: str, out_path: str = GENBANK_PATH):
    """Fetch NC_000913.3 from NCBI via Biopython's Entrez. Requires network
    access and a valid email (NCBI's usage policy). This sandbox cannot
    reach NCBI -- run this cell on your local Colab runtime instead, which
    does have internet access."""
    from Bio import Entrez
    Entrez.email = email
    print(f"Fetching {REFERENCE_ACCESSION} from NCBI (this can take a minute)...")
    handle = Entrez.efetch(
        db="nucleotide", id=REFERENCE_ACCESSION, rettype="gbwithparts", retmode="text"
    )
    with open(out_path, "w") as f:
        f.write(handle.read())
    print(f"Saved to {out_path}")
    return out_path


def load_reference_genome(genbank_path: str = GENBANK_PATH):
    """Load a local GenBank record. Raises FileNotFoundError with a clear
    pointer rather than silently returning None."""
    if not os.path.exists(genbank_path):
        raise FileNotFoundError(
            f"{genbank_path!r} not found. Either call download_reference_genome() "
            "on a machine with internet access, or download NC_000913.3 manually "
            "from https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3 (GenBank "
            "flat-file format) and point GENBANK_PATH at it."
        )
    print(f"Loading {genbank_path}...")
    record = SeqIO.read(genbank_path, "genbank")
    print(f"Loaded {record.id}: {len(record.seq)} bp, "
          f"{sum(1 for f in record.features if f.type == 'gene')} gene features.")
    return record


def build_bnumber_index(record) -> dict:
    """Map b-number (locus_tag, e.g. 'b3241') -> {start, end, strand, gene}
    (0-based, end-exclusive, matching Biopython's FeatureLocation
    convention) from a loaded GenBank record."""
    index = {}
    for feature in record.features:
        if feature.type != "gene":
            continue
        locus_tags = feature.qualifiers.get("locus_tag", [])
        if not locus_tags:
            continue
        bnum = locus_tags[0]
        index[bnum] = {
            "start": int(feature.location.start),
            "end": int(feature.location.end),
            "strand": feature.location.strand,  # +1 or -1
            "gene": feature.qualifiers.get("gene", [bnum])[0],
        }
    print(f"Indexed {len(index)} genes by b-number.")
    return index


# %% CELL 4 — Coordinate resolution and window extraction
def resolve_genomic_coordinate(bnum: str, offset: int, bnumber_index: dict,
                                offset_mode: str = "sense_relative") -> dict:
    """Convert (b-number, offset) into an absolute genome coordinate + strand.

    offset_mode="sense_relative" (default): offset counts from the gene's
      start codon in the direction of transcription -- i.e. genome_coord =
      gene_start + offset on the + strand, or gene_end - offset on the -
      strand. This is the more common convention in guide-library papers
      but is NOT confirmed against this project's data.
    offset_mode="absolute": offset is a plain left-to-right genome offset
      from gene_start regardless of strand.

    Try both against a few known guides and see which one reconstructs a
    protospacer that matches your matrix's existing positional encoding --
    that's the actual validation this script can't do without the CSV.
    """
    if bnum not in bnumber_index:
        raise KeyError(f"b-number {bnum!r} not found in the genome annotation index.")
    gene = bnumber_index[bnum]
    strand = gene["strand"]

    if offset_mode == "absolute":
        coord = gene["start"] + offset
    elif offset_mode == "sense_relative":
        coord = gene["start"] + offset if strand == 1 else gene["end"] - offset
    else:
        raise ValueError(f"Unknown offset_mode: {offset_mode!r}")

    return {"bnum": bnum, "gene": gene["gene"], "coord": coord, "strand": strand}


def extract_locus_window(record, coord: int, strand: int,
                          protospacer_len: int = 20,
                          flank_window_nt: int = FLANK_WINDOW_NT) -> dict:
    """Pull the protospacer + PAM + +/-flank_window_nt around a resolved
    coordinate, oriented to the guide's strand (reverse-complemented if
    strand == -1 so position 1 of the returned sequence is always the
    5' end of the protospacer as Cas9 sees it)."""
    genome_len = len(record.seq)
    lo = max(0, coord - flank_window_nt)
    hi = min(genome_len, coord + protospacer_len + 3 + flank_window_nt)  # +3 for PAM
    window = record.seq[lo:hi]
    if strand == -1:
        window = window.reverse_complement()
    return {
        "coord": coord,
        "strand": strand,
        "window_start": lo,
        "window_end": hi,
        "sequence": str(window),
    }


def build_locus_table(sgRNAID_series, bnumber_index, record,
                       offset_mode: str = "sense_relative") -> pd.DataFrame:
    """End-to-end: sgRNAID -> parsed fields -> genome coordinate -> locus
    window, for a whole Series of sgRNAIDs. Rows that fail at any step are
    kept with an `error` column rather than dropped silently, so you can
    see exactly how much of the dataset the coordinate lookup actually
    covers before building features on top of it."""
    rows = []
    for sid in sgRNAID_series:
        row = {"sgRNAID": sid}
        try:
            parsed = parse_sgRNAID(sid)
            row.update(parsed)
            resolved = resolve_genomic_coordinate(
                parsed["bnum"], parsed["offset"], bnumber_index, offset_mode
            )
            row.update(resolved)
            locus = extract_locus_window(record, resolved["coord"], resolved["strand"])
            row.update(locus)
            row["error"] = None
        except (ValueError, KeyError) as e:
            row["error"] = str(e)
        rows.append(row)
    df = pd.DataFrame(rows)
    n_ok = df["error"].isna().sum()
    print(f"Resolved {n_ok}/{len(df)} sgRNAIDs to a genome locus.")
    return df


# ============================================================
# DEMO
# ============================================================
if __name__ == "__main__":
    print("Demo: parsing the literature review's own example sgRNAID.")
    print(parse_sgRNAID("aaeAb3241_107_Cas9"))
    print(
        "\nCoordinate resolution and window extraction need a local "
        "NC_000913.3 GenBank file (see load_reference_genome's docstring) "
        "and your real sgRNAID values (to pick offset_mode and confirm the "
        "parser) -- neither is available in this sandbox. Run this on your "
        "local Colab runtime: \n"
        "  1. record = download_reference_genome(email='you@example.com')  "
        "# or load_reference_genome('NC_000913.3.gb') if already downloaded\n"
        "  2. bnumber_index = build_bnumber_index(record)\n"
        "  3. validate_sgRNAID_parsing(df['sgRNAID'])   # check the format assumption\n"
        "  4. locus_df = build_locus_table(df['sgRNAID'], bnumber_index, record)\n"
        "  5. Compare locus_df['sequence'][:20] (or its reverse complement) "
        "against whatever the existing p1..p20 matrix columns encode, for a "
        "handful of guides, under both offset_mode values -- whichever "
        "reconstructs the known protospacer is the right convention."
    )
