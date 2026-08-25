from __future__ import annotations

import argparse
import hashlib
import sqlite3
from dataclasses import dataclass
from pathlib import Path


SCHEMA_VERSION = "0.2.0"
EXPECTED_MIN_CYTOBAND_ROWS = 800


@dataclass(frozen=True, slots=True)
class BuildInput:
    genome_build: str
    common_name: str
    path: Path
    source_url: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_chromosome(chromosome: str) -> str | None:
    value = chromosome.removeprefix("chr")
    if value in {str(item) for item in range(1, 23)} | {"X", "Y"}:
        return f"chr{value}"
    return None


def read_cytobands(build_input: BuildInput) -> list[tuple[str, int, int, str, str]]:
    rows: list[tuple[str, int, int, str, str]] = []
    for line_number, raw_line in enumerate(
        build_input.path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            fields = line.split()
        if len(fields) < 5:
            raise ValueError(
                f"{build_input.path}:{line_number}: expected at least five columns"
            )
        chromosome = normalize_chromosome(fields[0])
        if chromosome is None:
            continue
        start = int(fields[1])
        end = int(fields[2])
        if start < 0 or end <= start:
            raise ValueError(f"{build_input.path}:{line_number}: invalid interval")
        rows.append((chromosome, start, end, fields[3], fields[4]))
    if len(rows) < EXPECTED_MIN_CYTOBAND_ROWS:
        raise ValueError(
            f"{build_input.path}: only {len(rows)} supported cytoband rows; "
            f"expected at least {EXPECTED_MIN_CYTOBAND_ROWS}"
        )
    return rows


def create_schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;

        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );

        CREATE TABLE sources (
            source_id INTEGER PRIMARY KEY,
            source_key TEXT NOT NULL UNIQUE,
            title TEXT NOT NULL,
            source_type TEXT NOT NULL,
            url TEXT,
            license_note TEXT,
            sha256 TEXT,
            notes TEXT
        );

        CREATE TABLE genome_builds (
            build_id INTEGER PRIMARY KEY,
            build_name TEXT NOT NULL UNIQUE,
            common_name TEXT NOT NULL UNIQUE,
            source_id INTEGER,
            record_count INTEGER NOT NULL DEFAULT 0,
            FOREIGN KEY(source_id) REFERENCES sources(source_id)
        );

        CREATE TABLE cytobands (
            cytoband_id INTEGER PRIMARY KEY,
            build_id INTEGER NOT NULL,
            chromosome TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            band TEXT NOT NULL,
            stain TEXT NOT NULL,
            arm TEXT NOT NULL CHECK(arm IN ('p','q')),
            FOREIGN KEY(build_id) REFERENCES genome_builds(build_id),
            UNIQUE(build_id, chromosome, start, end, band)
        );

        CREATE TABLE centromeres (
            centromere_id INTEGER PRIMARY KEY,
            build_id INTEGER NOT NULL,
            chromosome TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            source_method TEXT NOT NULL,
            FOREIGN KEY(build_id) REFERENCES genome_builds(build_id),
            UNIQUE(build_id, chromosome)
        );

        CREATE TABLE event_types (
            event_type_id INTEGER PRIMARY KEY,
            code TEXT NOT NULL UNIQUE,
            name TEXT NOT NULL,
            template TEXT,
            auto_render_supported INTEGER NOT NULL CHECK(auto_render_supported IN (0,1)),
            notes TEXT
        );

        CREATE TABLE rules (
            rule_id INTEGER PRIMARY KEY,
            rule_key TEXT NOT NULL UNIQUE,
            rule_group TEXT NOT NULL,
            rule_summary TEXT NOT NULL,
            implementation_note TEXT,
            source_id INTEGER,
            source_locator TEXT,
            FOREIGN KEY(source_id) REFERENCES sources(source_id)
        );

        CREATE TABLE examples (
            example_id INTEGER PRIMARY KEY,
            iscn_string TEXT NOT NULL,
            interpretation TEXT NOT NULL,
            expected_valid INTEGER NOT NULL CHECK(expected_valid IN (0,1)),
            source_id INTEGER,
            source_locator TEXT,
            FOREIGN KEY(source_id) REFERENCES sources(source_id)
        );

        CREATE TABLE external_validators (
            validator_id INTEGER PRIMARY KEY,
            package_name TEXT NOT NULL UNIQUE,
            pinned_version TEXT NOT NULL,
            license TEXT NOT NULL,
            role TEXT NOT NULL,
            enabled_by_default INTEGER NOT NULL CHECK(enabled_by_default IN (0,1)),
            notes TEXT
        );

        CREATE INDEX idx_cytobands_lookup
            ON cytobands(build_id, chromosome, start, end);
        CREATE INDEX idx_cytobands_band
            ON cytobands(build_id, chromosome, band);
        """
    )


def seed_static_data(connection: sqlite3.Connection) -> None:
    connection.executemany(
        "INSERT INTO metadata(key,value) VALUES (?,?)",
        [
            ("database_name", "ONTSeq Cytogenomic Knowledge Base"),
            ("schema_version", SCHEMA_VERSION),
            ("target_iscn_edition", "ISCN 2024"),
            (
                "scope",
                "Open/reference coordinate data plus reformulated implementation rules; "
                "does not reproduce the copyrighted ISCN 2024 publication.",
            ),
            ("clinical_status", "Research Use Only; expert review required"),
        ],
    )

    connection.execute(
        """
        INSERT INTO sources(source_key,title,source_type,url,license_note,notes)
        VALUES (?,?,?,?,?,?)
        """,
        (
            "lea_evers_thesis_2026",
            "Lea Evers MSc thesis: Establishing Oxford Nanopore Sequencing technology as a rapid method for karyotyping in clinical practice",
            "project thesis",
            None,
            "User-provided project source; only short reformulated rules/examples are stored.",
            "Relevant sections: ISCN fundamentals and implementation/cytoband merging.",
        ),
    )
    thesis_source_id = connection.execute(
        "SELECT source_id FROM sources WHERE source_key='lea_evers_thesis_2026'"
    ).fetchone()[0]

    connection.executemany(
        "INSERT INTO event_types(code,name,template,auto_render_supported,notes) VALUES (?,?,?,?,?)",
        [
            ("+", "whole chromosome gain", "+chr", 1, "Numerical event."),
            ("-", "whole chromosome loss", "-chr", 1, "Numerical event."),
            ("del", "deletion", "del(chr)(band_or_interval)", 1, "Subset renderer."),
            ("dup", "duplication", "dup(chr)(band_or_interval)", 1, "Subset renderer."),
            ("inv", "inversion", "inv(chr)(breakpoint_interval)", 1, "Subset renderer."),
            ("t", "translocation", "t(chr1;chr2)(bp1;bp2)", 1, "Subset renderer."),
            ("ins", "insertion", "ins(...)", 0, "Stored as event; not safely auto-rendered yet."),
            ("der", "derivative chromosome", "der(...)", 0, "Future renderer extension."),
            ("add", "additional material", "add(chr)(band)", 0, "Future renderer extension."),
        ],
    )

    rules = [
        (
            "general.header_order",
            "general",
            "Render chromosome count first, then sex chromosome complement, then abnormalities.",
            "Renderer constructs a fixed header before event fragments.",
            "p. 18",
        ),
        (
            "general.no_comma_spaces",
            "formatting",
            "Do not emit spaces around comma separators.",
            "Renderer joins fragments with literal commas.",
            "p. 18",
        ),
        (
            "event.numeric_sign",
            "numerical",
            "Represent whole-chromosome gain/loss with plus/minus signs.",
            "CHROMOSOME_GAIN/LOSS map to signed chromosome fragments.",
            "p. 18",
        ),
        (
            "ordering.numeric_before_structural",
            "ordering",
            "For the same chromosome, numerical events precede structural events.",
            "Event sort key assigns numerical events a lower category rank.",
            "p. 18",
        ),
        (
            "cytoband.coordinate_order",
            "cytoband",
            "Determine band adjacency from genomic coordinates rather than alphanumeric sorting.",
            "CytobandIndex is ordered by genomic start coordinate.",
            "pp. 40-41",
        ),
        (
            "cytoband.no_pq_merge",
            "cytoband",
            "Do not collapse simple deletion/duplication intervals across the centromere.",
            "Cross-arm del/dup events are omitted from the automatic subset.",
            "p. 41",
        ),
        (
            "report.current_edition",
            "reporting",
            "State the target ISCN edition in report metadata.",
            "ISCNProposal.standard_edition remains locked to ISCN 2024.",
            "p. 27",
        ),
    ]
    connection.executemany(
        """
        INSERT INTO rules(rule_key,rule_group,rule_summary,implementation_note,source_id,source_locator)
        VALUES (?,?,?,?,?,?)
        """,
        [(*row[:4], thesis_source_id, row[4]) for row in rules],
    )

    connection.executemany(
        """
        INSERT INTO examples(iscn_string,interpretation,expected_valid,source_id,source_locator)
        VALUES (?,?,?,?,?)
        """,
        [
            ("47,XY,+8", "Male karyotype with an additional chromosome 8.", 1, thesis_source_id, "p. 18"),
            ("44,X,-X,-3", "Female karyotype with loss of X and chromosome 3.", 1, thesis_source_id, "p. 18"),
            ("46,XX,inv(2)(p23p13)", "Inversion of chromosome 2.", 1, thesis_source_id, "p. 19"),
            ("46,XX,t(9;22)(q34;q11.2)", "Translocation between chromosomes 9 and 22.", 1, thesis_source_id, "p. 20"),
        ],
    )

    connection.execute(
        """
        INSERT INTO external_validators(package_name,pinned_version,license,role,enabled_by_default,notes)
        VALUES (?,?,?,?,?,?)
        """,
        (
            "iscn-authenticator",
            "0.2.1",
            "MIT",
            "Independent parser/rule-engine cross-check for generated ISCN strings",
            0,
            "Optional until separately validated against the ONTSeq fixture corpus.",
        ),
    )


def add_build(connection: sqlite3.Connection, build_input: BuildInput) -> None:
    rows = read_cytobands(build_input)
    checksum = sha256_file(build_input.path)
    source_key = f"ucsc_{build_input.common_name}_cytoband"
    connection.execute(
        """
        INSERT INTO sources(source_key,title,source_type,url,license_note,sha256,notes)
        VALUES (?,?,?,?,?,?,?)
        """,
        (
            source_key,
            f"UCSC {build_input.common_name} cytoband table",
            "open genomic reference",
            build_input.source_url,
            "Reference data provenance retained; review upstream terms before redistribution.",
            checksum,
            "Imported from a five-column UCSC-style cytoband table.",
        ),
    )
    source_id = connection.execute(
        "SELECT source_id FROM sources WHERE source_key=?", (source_key,)
    ).fetchone()[0]
    connection.execute(
        """
        INSERT INTO genome_builds(build_name,common_name,source_id,record_count)
        VALUES (?,?,?,?)
        """,
        (build_input.genome_build, build_input.common_name, source_id, len(rows)),
    )
    build_id = connection.execute(
        "SELECT build_id FROM genome_builds WHERE build_name=?", (build_input.genome_build,)
    ).fetchone()[0]
    connection.executemany(
        """
        INSERT INTO cytobands(build_id,chromosome,start,end,band,stain,arm)
        VALUES (?,?,?,?,?,?,?)
        """,
        [
            (build_id, chromosome, start, end, band, stain, band[0])
            for chromosome, start, end, band, stain in rows
        ],
    )

    chromosomes = sorted({row[0] for row in rows})
    for chromosome in chromosomes:
        acen = [row for row in rows if row[0] == chromosome and row[4] == "acen"]
        if not acen:
            continue
        connection.execute(
            """
            INSERT INTO centromeres(build_id,chromosome,start,end,source_method)
            VALUES (?,?,?,?,?)
            """,
            (
                build_id,
                chromosome,
                min(row[1] for row in acen),
                max(row[2] for row in acen),
                "span of UCSC acen cytobands",
            ),
        )


def build_database(output: Path, inputs: list[BuildInput], bootstrap_only: bool) -> None:
    if output.exists():
        output.unlink()
    connection = sqlite3.connect(output)
    try:
        create_schema(connection)
        seed_static_data(connection)
        if not bootstrap_only:
            if len(inputs) != 2:
                raise ValueError("provide both hg19 and hg38 cytoband files")
            for build_input in inputs:
                add_build(connection, build_input)
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the ONTSeq ISCN knowledge-base SQLite file")
    parser.add_argument("--hg19", type=Path)
    parser.add_argument("--hg38", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--bootstrap-only",
        action="store_true",
        help="Build schema/rules/examples without cytoband rows (development only).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inputs: list[BuildInput] = []
    if args.hg19:
        inputs.append(
            BuildInput(
                genome_build="GRCh37",
                common_name="hg19",
                path=args.hg19,
                source_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz",
            )
        )
    if args.hg38:
        inputs.append(
            BuildInput(
                genome_build="GRCh38",
                common_name="hg38",
                path=args.hg38,
                source_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
            )
        )
    build_database(args.output, inputs, args.bootstrap_only)
    print(args.output)


if __name__ == "__main__":
    main()
