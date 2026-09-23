// Presentation-only model calculations. Never change normalized events or release gates.
const finite = value => typeof value === 'number' && Number.isFinite(value);
export const validFit = fit => fit && finite(fit.cellularity) && fit.cellularity > 0 && fit.cellularity <= 1 && finite(fit.ploidy) && fit.ploidy > 0 && finite(fit.fit_error) && fit.fit_error >= 0;
export function ratioAt(copies, fraction, ploidy) {
  if (!finite(copies) || copies < 0 || !finite(fraction) || fraction <= 0 || fraction > 1 || !finite(ploidy) || ploidy <= 0) return null;
  return (fraction * copies + (1 - fraction) * 2) / (fraction * ploidy + (1 - fraction) * 2);
}
export function copyAt(ratio, fraction, ploidy) {
  if (!finite(ratio) || ratio < 0 || !finite(fraction) || fraction <= 0 || fraction > 1 || !finite(ploidy) || ploidy <= 0) return null;
  const copies = (ratio * (fraction * ploidy + (1 - fraction) * 2) - (1 - fraction) * 2) / fraction;
  return copies < -1e-9 ? null : Math.max(0, copies);
}
export function dilutedRatio(copies, ploidy, fraction) {
  if (!finite(copies) || copies < 0 || !finite(ploidy) || ploidy <= 0 || !finite(fraction) || fraction < 0 || fraction > 1) return null;
  return (fraction * copies + (1 - fraction) * 2) / (fraction * ploidy + (1 - fraction) * 2);
}
export function solutions(fit) {
  if (!fit) return [];
  const list = [{...fit, pipeline: true}, ...(fit.alternatives || []).map(row => ({...row, pipeline:false}))];
  return list.filter(validFit).filter((row,index,rows) => rows.findIndex(other => other.cellularity === row.cellularity && other.ploidy === row.ploidy && other.fit_error === row.fit_error) === index);
}
export function moduleState(name, modules, requested) {
  const recorded = modules.find(row => row.name === name);
  return recorded || {name, status:requested.includes(name) ? 'NOT_RECORDED' : 'NOT_REQUESTED', reason:requested.includes(name) ? 'Angefordert, aber kein Ausführungszustand im Ergebnis vorhanden.' : 'Für diesen Lauf nicht angefordert.'};
}
export const moduleLabels = {qc:'Qualitätskontrolle',cnv:'Kopienzahl',sv:'Strukturvarianten',fusion:'Fusionsprüfung',iscn:'ISCN',report:'Bericht',small_variants:'Kleine Varianten',methylation:'Methylierung'};
export const statusLabels = {COMPLETED:'gelaufen',NO_CALL:'ohne belastbaren Call',FAILED:'fehlgeschlagen',NOT_RUN:'nicht gelaufen',NOT_REQUESTED:'nicht angefordert',NOT_RECORDED:'nicht dokumentiert'};
export const number = (value, digits = 2) => value === null || value === undefined || (typeof value === 'number' && !Number.isFinite(value)) ? '—' : typeof value === 'number' ? value.toLocaleString('de-DE', {maximumFractionDigits:digits}) : String(value);
export const percent = value => value === null || value === undefined ? '—' : `${number(value * 100, 1)} %`;
export function modelSummary(rows, ploidy) {
  const measured = rows.filter(row => finite(row.copies));
  if (!measured.length) return 'Keine auswertbaren Chromosomenwerte für diese Modellansicht.';
  const changed = measured.filter(row => Math.abs(row.copies - ploidy) > .5);
  const detail = changed.map(row => `${row.chromosome}: ${number(row.copies)} Kopien`).join(' · ') || 'Keine Abweichung über 0,5 Kopien von der Modellploidie in den verfügbaren Chromosomenmitteln.';
  return `${detail} ${measured.length} Chromosomen auswertbar${rows.length > measured.length ? ` · ${rows.length - measured.length} nicht bestimmbar` : ''}. Deskriptive Modellansicht, kein normaler Karyotyp abgeleitet.`;
}

export function marlinAssessmentFacts(report) {
  const fraction = report.feature_summary?.observed_fraction;
  return [
    ['Anteil beobachteter Modell-CpGs', fraction == null ? 'nicht verfügbar' : `${number(fraction * 100,6)} %`],
    ['Rohe Modellscore-Schwelle', number(report.model_score_threshold,2)],
    ['Modellscore-Schwelle erreicht', report.model_score_threshold_met == null ? 'nicht verfügbar' : report.model_score_threshold_met ? 'ja · nur Modellscore' : 'nein · nur Modellscore'],
    ['Assay-Beurteilbarkeit', `${report.assay_assessability ?? 'NOT_ESTABLISHED'} · Bewertung nicht validiert`],
  ];
}
