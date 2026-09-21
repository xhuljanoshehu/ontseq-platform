from __future__ import annotations

import html
from pathlib import Path

from .cnv.qdnaseq import QDNAseqCallReport
from .report_interactive import read_plot
from .report_plots import CnvChromosomeBar, cnv_genome_svg


def cnv_html_section(root: Path, cnv: QDNAseqCallReport) -> str:
    fit_rows = "".join(
        "<tr>"
        f"<td>{fit.bin_size_kbp}</td><td>{fit.cellularity:.3f}</td>"
        f"<td>{fit.ploidy:.3f}</td><td>{fit.fit_error:.6g}</td>"
        f"<td>{fit.segment_count}</td></tr>"
        for fit in sorted(cnv.fits, key=lambda value: value.bin_size_kbp)
    )
    chromosome_rows = "".join(
        "<tr>"
        f"<td>{html.escape(chromosome.chromosome)}</td>"
        f"<td>{chromosome.median_copy_number:.3f}</td>"
        f"<td>{chromosome.rounded_copy_number}</td>"
        f"<td>{chromosome.agreeing_bins}/{chromosome.contributing_bins}</td>"
        f"<td>{chromosome.min_copy_number:.3f}–{chromosome.max_copy_number:.3f}</td>"
        "</tr>"
        for chromosome in cnv.chromosome_consensus
    )
    bars = [
        CnvChromosomeBar(
            chromosome=chromosome.chromosome,
            median_cn=chromosome.median_copy_number,
            min_cn=chromosome.min_copy_number,
            max_cn=chromosome.max_copy_number,
            rounded_cn=chromosome.rounded_copy_number,
        )
        for chromosome in cnv.chromosome_consensus
    ]
    genome_svg = cnv_genome_svg(
        bars,
        title="Genome-wide copy-number overview",
        baseline=cnv.primary_fit.ploidy,
    )
    genome_figure = (
        "<h3>Genome-wide copy-number overview</h3><figure class='plot' "
        "style='margin:14px 0'>"
        + genome_svg
        + "<figcaption style='color:#5e687a;font-size:12px;margin-top:6px'>Median copy "
        "number per chromosome from the multi-bin consensus; whiskers show each "
        "chromosome's min–max range across bin sizes; the dashed line marks the "
        "fitted ACE ploidy. Autosomes only — no X/Y statement is made or implied. "
        "Descriptive technical evidence, not an adequacy or reportability assessment."
        "</figcaption></figure>"
        if genome_svg
        else ""
    )
    images: dict[str, str] = {}
    for label, name in (
        ("Kopienzahlprofil · QDNAseq + ACE", cnv.primary_fit.copy_number_plot),
        ("ACE · Tumoranteil und Ploidie", cnv.primary_fit.fit_plot),
    ):
        image = read_plot(root, name)
        if image is None:
            continue
        images[name] = (
            f"<figure class='plot'><h3>{html.escape(label)}</h3>"
            f"<img alt='{html.escape(label)}' style='max-width:100%;height:auto' "
            f"src='{image}'>"
            f"<figcaption>{html.escape(name)} · Originalabbildung aus dem Lauf"
            "</figcaption></figure>"
        )
    profile_image = images.get(cnv.primary_fit.copy_number_plot, "")
    fit_images = images.get(cnv.primary_fit.fit_plot, "")
    return (
        "<section id='cnv'><h2>Kopienzahlprofil · QDNAseq + ACE</h2>"
        "<p class='muted'>Copy-number analysis — QDNAseq + ACE. "
        "Modellwerte aus dem Lauf; keine eigenständige klinische Freigabe.</p>"
        "<div class='cnv-facts'>"
        f"<div><span>Primäre Bin-Größe</span><strong>{cnv.primary_fit.bin_size_kbp} kbp"
        "</strong></div>"
        f"<div><span>Zellularität</span><strong>{cnv.primary_fit.cellularity:.3f}"
        "</strong></div>"
        f"<div><span>Ploidie</span><strong>{cnv.primary_fit.ploidy:.3f}</strong></div>"
        f"<div><span>Anpassungsfehler</span><strong>{cnv.primary_fit.fit_error:.6g}"
        "</strong></div></div>"
        + profile_image
        + genome_figure
        + "</section><section id='cnv-fit'><h2>Modellanpassung über die Bin-Größen</h2>"
        + fit_images
        + "<h3>Multi-resolution fits</h3><div class='table-wrap'><table><thead><tr>"
        "<th>Bin (kbp)</th><th>Cellularity</th><th>Ploidy</th><th>Fit error</th>"
        f"<th>Segments</th></tr></thead><tbody>{fit_rows}</tbody></table></div>"
        "</section><section id='cnv-chromosomes'><h2>Kopienzahl je Chromosom</h2>"
        "<h3>Chromosome-level consensus</h3><div class='table-wrap'><table><thead><tr>"
        "<th>Chromosome</th><th>Median CN</th><th>Rounded CN</th><th>Agreement</th>"
        f"<th>Range</th></tr></thead><tbody>{chromosome_rows}</tbody></table></div></section>"
    )
