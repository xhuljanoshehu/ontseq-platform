import React, { useEffect, useState } from 'react';
import { createRoot } from 'react-dom/client';
import { copyAt, ratioAt, dilutedRatio, solutions, moduleState, modelSummary, moduleLabels as baseModuleLabels, statusLabels, number as n, percent } from './model.js';
import './report.css';

const data = JSON.parse(document.getElementById('ontseq-befund-data').textContent);
const original = document.getElementById('ontseq-static-report');
// Reuse escaped, server-rendered evidence sections. Do not interpret arbitrary HTML input.
const evidence = {};
for (const id of ['coverage','methylation','iscn','provenance','events','qc','sv-review']) {
  const node = original.querySelector(`[id="${id}"]`);
  if (!node) continue;
  const clone = node.cloneNode(true);
  clone.querySelectorAll('script').forEach(script => script.remove());
  const ids = new Map([...clone.querySelectorAll('[id]')].map(element => [element.id, `bf-evidence-${id}-${element.id}`]));
  clone.querySelectorAll('[id]').forEach(element => {element.id=ids.get(element.id);});
  clone.querySelectorAll('*').forEach(element => {
    for (const attr of [...element.attributes]) {
      if (attr.name.startsWith('on')) {element.removeAttribute(attr.name);continue;}
      const value=attr.value.replace(/url\(#([^)]*)\)/g,(match,key)=>ids.has(key)?`url(#${ids.get(key)})`:match);
      if ((attr.name==='href'||attr.name==='xlink:href') && value.startsWith('#')) {
        if(ids.has(value.slice(1))) element.setAttribute(attr.name,`#${ids.get(value.slice(1))}`);
        else element.removeAttribute(attr.name);
      } else element.setAttribute(attr.name,value);
    }
  });
  evidence[id] = clone.innerHTML;
}
const moduleLabels = {...baseModuleLabels, marlin: "MARLIN-Klassifikation"};
const names = Object.keys(moduleLabels);
const labelType = type => ({insertion:'Insertion',deletion:'Deletion',duplication:'Duplikation',inversion:'Inversion',translocation:'Translokation',chromosome_gain:'Chromosomenzugewinn',chromosome_loss:'Chromosomenverlust'}[type] || type);
function Section({id,title,note,tone='',children}) {
  return <section className={`bf-section ${tone}`} id={`befund-${id}`}><h2>{title}</h2>{note && <p className="bf-note">{note}</p>}{children}</section>;
}
function Empty({children}) { return <div className="bf-empty">{children || 'Für diesen Lauf liegen keine entsprechenden Messwerte vor.'}</div>; }
function Native({id}) { return evidence[id] ? <div className="bf-native" dangerouslySetInnerHTML={{__html:evidence[id]}}/> : <Empty/>; }
function Facts({rows}) { return <div className="bf-facts">{rows.map(([label,value]) => <div key={label}><span>{label}</span><strong>{value}</strong></div>)}</div>; }
function Status({value}) { return <span className={`bf-status s-${value}`}>{statusLabels[value] || value}</span>; }
function Marlin({report}) {
  if (!report) return <Empty>Für diesen Lauf kein MARLIN-Ausführungsnachweis vorhanden.</Empty>;
  const feature = report.feature_summary;
  return <><p><Status value={report.status}/> · {report.reason}</p>
    <p className="bf-alert">{report.validation_status} · Nur für Forschungszwecke; nicht klinisch validiert. Modellscores sind keine Erkrankungswahrscheinlichkeiten.</p>
    {report.decision === 'UNKNOWN' && <p className="bf-alert">UNKNOWN · Keine hinreichend sichere Klassifikation. Eine führende Modellklasse ist keine bestätigte Klasse.</p>}
    <Facts rows={[["Entscheidung",report.decision ?? 'nicht verfügbar'],["Beobachtete Modell-CpGs",n(feature?.observed_model_feature_count,0)],["Erwartete Modell-CpGs",n(feature?.expected_feature_count,0)],["Explizit fehlende CpGs",n(feature?.explicit_na_feature_count,0)],["Nicht beobachtete Modell-CpGs",n(feature?.absent_feature_count,0)],["Führende Modellklasse",report.top_class ?? 'nicht verfügbar'],["Führender Modellscore",n(report.top_class_score,6)],["Konfidenzschwelle",n(report.confidence_threshold,2)]]}/>
    {[["Klassen",report.class_scores],["Familien",report.family_scores],["Linien",report.lineage_scores]].map(([title,rows])=>rows?.length>0 && <div key={title} className="bf-table"><h3>{title} · Modellscores</h3><table><thead><tr><th>Rang</th><th>Bezeichnung</th><th>Modellscore</th></tr></thead><tbody>{rows.map((row,index)=><tr key={row.label}><td>{index+1}</td><td>{row.label}</td><td>{n(row.score,6)}</td></tr>)}</tbody></table></div>)}
    {[...(report.warnings || []),...(report.limitations || [])].map((text,index)=><p key={index} className="bf-alert">{text}</p>)}
    <details><summary>MARLIN-Werkzeuge, Parameter und Prüfsummen</summary><pre>{JSON.stringify({schema_version:report.schema_version,adapter_version:report.adapter_version,run_id:report.run_id,tools:report.tools,parameters:report.parameters,input_fingerprints:report.input_fingerprints,installation_signature:report.installation_signature,feature_summary:feature},null,2)}</pre></details>
  </>;
}
function Chromosomes({rows,baseline=2,model=false}) {
  return <div className="bf-chromosomes">{Array.from({length:22},(_,index)=>`chr${index+1}`).map(chr => {
    const row = rows.find(item=>item.chromosome===chr); const cn = row?.copies;
    const tone = cn == null ? 'missing' : cn > baseline + .35 ? 'gain' : cn < baseline - .35 ? 'loss' : '';
    return <div className={`bf-chromosome ${tone}`} key={chr}><span>{chr.replace('chr','')}</span><strong data-chromosome={chr}>{n(cn)}</strong><small>{cn == null ? 'nicht bestimmbar' : model ? 'Modellkopien' : 'Kopien'}</small>{row?.agreement && <small>{row.agreement}</small>}</div>;
  })}</div>;
}
function EventTable({events}) {
  const [filter,setFilter] = useState('');
  const rows = events.filter(event=>JSON.stringify(event).toLocaleLowerCase().includes(filter.toLocaleLowerCase()));
  return <><label className="bf-search">Ereignisse durchsuchen<input value={filter} onChange={e=>setFilter(e.target.value)} placeholder="Chromosom, Gen oder Ereigniskennung"/></label><p>{rows.length} von {events.length} Ereignissen</p><div className="bf-table"><table><thead><tr><th>Kennung / Typ</th><th>Position / Gene</th><th>Nachweis</th><th>Bewertungsgrenze</th></tr></thead><tbody>{rows.map(event=><tr key={event.event_id}><td><strong>{event.event_id}</strong><br/>{labelType(event.event_type)}</td><td className="bf-mono">{event.primary_locus}{event.secondary_locus && <><br/>{event.secondary_locus}</>}<br/>{event.genes.join(' · ') || 'Keine Genannotation'}</td><td>{event.evidence.map((ev,i)=><div key={i}>{ev.caller} {ev.caller_version} · {n(ev.support_reads,0)} Reads · VAF {percent(ev.variant_allele_fraction)}</div>)}</td><td><span>{event.confidence}</span><p>{event.reportability_text}</p><details><summary>Details und Annotationen</summary>{event.notes.map((note,i)=><p key={i}>{note}</p>)}{event.annotations.map((annotation,i)=><p key={i}>{annotation.source_id} · {annotation.assertion} · {annotation.scope_note}</p>)}</details></td></tr>)}</tbody></table>{!rows.length && <Empty>Keine passenden Ereignisse. Dies ist keine Aussage über biologische Abwesenheit.</Empty>}</div></>;
}
function App() {
  const view=data.view, cnv=data.cnv;
  const [bin,setBin]=useState(cnv?.primary_fit.bin_size_kbp ?? null);
  const [selected,setSelected]=useState(0);
  const [fraction,setFraction]=useState(null);
  const [vaf,setVaf]=useState(.25);
  const [anchor,setAnchor]=useState('');
  const [showAudit,setShowAudit]=useState(false);
  const fit=cnv?.fits.find(row=>row.bin_size_kbp===bin) || cnv?.primary_fit;
  const alternatives=solutions(fit), solution=alternatives[selected] || alternatives[0];
  const primary=cnv?.primary_fit;
  const chromosomes=(cnv?.chromosomes || []).filter(row=>/^chr([1-9]|1[0-9]|2[0-2])$/.test(row.chromosome));
  // Consensus copy numbers are reconstructed under the recorded primary fit. They are
  // not raw bin ratios; this distinction is printed next to every exploratory view.
  const ratios=chromosomes.map(row=>({...row,ratio:ratioAt(row.median_copy_number,primary?.cellularity,primary?.ploidy)}));
  const modeled=ratios.map(row=>({...row,copies:solution ? copyAt(row.ratio,solution.cellularity,solution.ploidy) : null}));
  const dilution=fraction ?? solution?.cellularity ?? 0;
  const anchorRatio=ratios.find(row=>row.chromosome===anchor)?.ratio;
  const anchorPloidy=anchorRatio > 0 && vaf > 0 ? (2/anchorRatio - 2*(1-2*vaf))/(2*vaf) : null;
  const anchorRows=ratios.map(row=>({...row,copies:copyAt(row.ratio,2*vaf,anchorPloidy)}));
  const modules=names.map(name=>moduleState(name,view.modules,data.requested_modules));
  const cnvEvents=view.events.filter(event=>event.copy_number !== null);
  const svEvents=view.events.filter(event=>event.copy_number === null);
  const plot=data.plots[String(bin)];
  const methyl=modules.find(row=>row.name==='methylation');
  const qcEntries=view.qc_metrics;
  useEffect(()=>{document.body.classList.add('befund-ready');},[]);
  useEffect(()=>{
    let closed=[];
    const before=()=>{closed=[...document.querySelectorAll('#ontseq-befund-root details:not([open])')];closed.forEach(node=>{node.open=true;});};
    const after=()=>{closed.forEach(node=>{node.open=false;});closed=[];};
    window.addEventListener('beforeprint',before);window.addEventListener('afterprint',after);
    return ()=>{window.removeEventListener('beforeprint',before);window.removeEventListener('afterprint',after);};
  },[]);
  useEffect(()=>{original.classList.toggle('bf-audit-visible',showAudit);},[showAudit]);
  return <div className="bf-app">
    <div className="bf-ruo">RESEARCH USE ONLY · NICHT KLINISCH VALIDIERT</div>
    <header className="bf-header"><div className="bf-header-inner"><div className="bf-brand">ONTSeq <span>BEFUND · {view.pipeline_version}</span></div><div className="bf-heading"><div><p className="bf-kicker">Genomische Auswertung · Einzelprobe</p><h1>Auswertung Oxford-Nanopore-Sequenzierung</h1><p>Analytische Ergebnisse und ihre Nachweise. Ohne fachliche Vidierung nicht freigegeben.</p></div><div className="bf-release"><span>Freigabestatus</span><strong>{view.release_status}</strong><small>Keine elektronische Signatur</small></div></div><Facts rows={[["Probe",view.sample_id],["Lauf",view.run_id],["Genom",view.genome_build],["Assay",view.assay_mode],["Profil",view.analysis_profile],["Erstellt",view.created_at]]}/></div></header>
    <div className="bf-shell"><nav className="bf-nav" aria-label="Befundabschnitte">{[['key','Kernbefunde'],['execution','Ausführung'],['profile','Kopienzahl'],['solutions','Lösungsraum'],['chromosomes','Chromosomen'],['sv','Strukturvarianten'],['qc','QC'],['methylation','Methylierung'],['iscn','ISCN'],['review','Vidierung']].map(([id,title])=><a key={id} href={`#befund-${id}`}>{title}</a>)}<button onClick={()=>window.print()}>Drucken / PDF</button></nav>
    <main className="bf-main">
    <Section id="key" title="Kernbefunde" note="Normalisierte Ereignisse aus diesem Lauf. Eine alternative Modellansicht ändert die gespeicherten Befunde nicht.">
      <div className="bf-keycards"><div><span>Kopienzahlereignisse</span><strong>{n(cnvEvents.length,0)}</strong><p><Status value={modules.find(row=>row.name==='cnv').status}/></p>{solution && <small data-testid="selected-solution">Modellansicht: {percent(solution.cellularity)} · Ploidie {n(solution.ploidy)}</small>}</div><div><span>Strukturvarianten</span><strong>{n(svEvents.length,0)}</strong><p><Status value={modules.find(row=>row.name==='sv').status}/></p><small>Kandidaten zur fachlichen Prüfung</small></div><div><span>Qualitätskontrolle</span><strong>{view.qc_verdict}</strong><p>{view.qc_failed_gates.length} nicht erfüllte Kriterien</p><small>{view.analysis_intent} · {view.genome_build}</small></div></div>
      {solution && <div className="bf-model-strip">Ausgewählte Modellansicht · {modelSummary(modeled,solution.ploidy)}</div>}
      {view.alerts.map((alert,i)=><p key={i} className={`bf-alert ${alert.level}`}><strong>{alert.title}</strong> · {alert.detail}</p>)}
    </Section>
    <Section id="execution" title="Ausführungszustand · warum eine Stufe kein Ergebnis hat" note="Nicht angefordert, nicht gelaufen, fehlgeschlagen und NO_CALL sind unterschiedliche Zustände. Keiner davon bedeutet einen negativen biologischen Befund."><div className="bf-table"><table><thead><tr><th>Modul</th><th>Zustand</th><th>Begründung aus dem Lauf</th></tr></thead><tbody>{modules.map(row=><tr key={row.name}><td>{moduleLabels[row.name]}</td><td><Status value={row.status}/><small className="bf-code">{row.status}</small></td><td>{row.reason}</td></tr>)}</tbody></table></div></Section>
    <Section id="profile" title="Kopienzahlprofil · QDNAseq + ACE" note="Originalabbildungen des gewählten Laufs. X und Y werden im autosomalen CNV-Modell nicht bewertet.">{fit ? <><div className="bf-controls"><span>Bin-Größe</span>{cnv.fits.map(row=><button key={row.bin_size_kbp} aria-pressed={bin===row.bin_size_kbp} onClick={()=>{setBin(row.bin_size_kbp);setSelected(0);setFraction(null);}}>{row.bin_size_kbp} kbp</button>)}</div><Facts rows={[["Pipeline-Zellularität",percent(fit.cellularity)],["Pipeline-Ploidie",n(fit.ploidy)],["Anpassungsfehler",n(fit.fit_error,6)],["Segmente",n(fit.segment_count,0)]]}/>{plot?.copy_number ? <figure><img alt={`Original-Kopienzahlprofil ${bin} kbp`} src={plot.copy_number}/><figcaption>Unveränderte Originalgrafik · {bin} kbp · Pipeline-Anpassung</figcaption></figure> : <Empty>Die Originalgrafik wurde nicht mit diesem Bericht bereitgestellt.</Empty>}{plot?.fit && <details><summary>Originalgrafik der ACE-Anpassung</summary><img alt="ACE-Modellanpassung" src={plot.fit}/></details>}</> : <Empty>Kein QDNAseq-/ACE-Ergebnis für diesen Lauf verfügbar. Siehe Ausführungszustand.</Empty>}</Section>
    <Section id="solutions" title="Dieselben Daten, mehrere Lesarten" tone="bf-caution" note="Explorative Modellrechnung. Die Tabelle enthält die gespeicherte Pipeline-Lösung und verfügbare alternative Minima. Ein kleinerer Anpassungsfehler ist kein unabhängiger Nachweis für Tumoranteil oder Ploidie.">{alternatives.length ? <><div className="bf-table"><table><thead><tr><th>Modellansicht wählen</th><th>Tumoranteil</th><th>Ploidie</th><th>Fehler</th><th>Relativ zur Pipeline</th></tr></thead><tbody>{alternatives.map((row,i)=><tr key={i} className={i===selected?'bf-selected':''}><td><button data-solution={i} aria-pressed={i===selected} onClick={()=>{setSelected(i);setFraction(null);}}>{row.pipeline?'Pipeline-Auswahl':`Alternative ${i}`}</button></td><td>{percent(row.cellularity)}</td><td>{n(row.ploidy)}</td><td>{n(row.fit_error,6)}</td><td>{fit.fit_error>0 ? `${n(row.fit_error/fit.fit_error)} ×` : 'nicht bestimmbar'}</td></tr>)}</tbody></table></div><p className="bf-footnote">Die interaktiven Kopienzahlen werden aus dem Chromosomenkonsens unter der primären Anpassung rekonstruiert. Sie sind keine erneute Segmentierung und keine bin-spezifischen Messwerte. Nicht mögliche negative Kopienzahlen bleiben unbestimmbar. Die normalisierten Ereignisse, Originalgrafiken und ISCN bleiben unverändert.</p></> : <Empty>Keine gültigen Reinheits-/Ploidie-Lösungen verfügbar.</Empty>}</Section>
    <Section id="anchor" title="Baseline auf 2 · Anker aus der Allelfrequenz" tone="bf-model" note="Hypothetische Modellrechnung, keine gemessene VAF. Vorausgesetzt werden eine klonale heterozygote somatische Variante und eine unabhängig bestätigte kopienzahlneutrale Ankerregion.">{ratios.length ? <><div className="bf-controls"><label>Angenommener neutraler Anker<select value={anchor} onChange={e=>setAnchor(e.target.value)}><option value="">Bitte bewusst auswählen</option>{ratios.map(row=><option key={row.chromosome} value={row.chromosome}>{row.chromosome} · Verhältnis {n(row.ratio,3)}</option>)}</select></label><label>Hypothetische VAF · {percent(vaf)}<input aria-label="Hypothetische VAF" type="range" min=".01" max=".5" step=".01" value={vaf} onChange={e=>setVaf(Number(e.target.value))}/></label></div>{anchor ? <><Facts rows={[["Angenommener Tumoranteil f = 2 × VAF",percent(vaf*2)],["Abgeleitete Modellploidie",anchorPloidy>0?n(anchorPloidy):'nicht möglich']]}/><Chromosomes rows={anchorRows} baseline={2} model/></> : <Empty>Erst die ausdrückliche Auswahl eines hypothetischen neutralen Ankers aktiviert diese Rechnung.</Empty>}<p className="bf-footnote">Die Neutralität des ausgewählten Chromosoms wird hier nicht geprüft. Diese Modellansicht bestimmt weder eine Nachweisgrenze noch eine Risikogruppe.</p></> : <Empty>Ohne Chromosomenwerte lässt sich kein Anker-Modell darstellen.</Empty>}</Section>
    <Section id="dilution" title="In-silico-Verdünnung" note="Hypothetische Beimischung diploider Normalzellen bei festgehaltener ausgewählter Tumorploidie und Kopienzahl. Keine experimentelle Verdünnungsreihe.">{solution && modeled.length ? <><label className="bf-slider">Tumoranteil · <strong data-testid="dilution-value">{percent(dilution)}</strong><input aria-label="Tumoranteil der Verdünnung" type="range" min="0" max="1" step=".01" value={dilution} onChange={e=>setFraction(Number(e.target.value))}/></label><div className="bf-ratios">{modeled.map(row=>{const value=dilutedRatio(row.copies,solution.ploidy,dilution);return <div key={row.chromosome}><span>{row.chromosome}</span><div className="bf-ratio-track"><i style={{width:`${Math.max(0,Math.min(100,(value??0)/3*100))}%`}}/></div><strong>{n(value,3)}</strong></div>;})}</div><p className="bf-footnote">r = (f × CN + (1 − f) × 2) / (f × P + (1 − f) × 2). Zahlen zeigen das modellierte relative Signal; die Balken reichen bis Verhältnis 3. Bei 0 % Tumoranteil ergibt sich Verhältnis 1.</p></> : <Empty>Keine Modellwerte für eine Verdünnung verfügbar.</Empty>}</Section>
    <Section id="series" title="Verdünnungsreihe und technische Nachweisgrenze" note="Eine interaktive Simulation ersetzt keine technische Validierung."><Empty>In diesem Ergebnisvertrag ist keine experimentelle Verdünnungsreihe hinterlegt. Keine Nachweisgrenze ableitbar.</Empty></Section>
    <Section id="chromosomes" title="Kopienzahl je Chromosom" note={solution?`Aktuelle Modellansicht: ${percent(solution.cellularity)} Tumoranteil · Ploidie ${n(solution.ploidy)}. Nur Autosomen.`:'Kein vollständiger Chromosomenkonsens verfügbar.'}><Chromosomes rows={modeled} baseline={solution?.ploidy ?? 2} model/>{chromosomes.length>0 && <details><summary>Unveränderter Konsens über die Bin-Größen</summary><div className="bf-table"><table><thead><tr><th>Chromosom</th><th>Median CN</th><th>Spanne</th><th>Übereinstimmung</th></tr></thead><tbody>{chromosomes.map(row=><tr key={row.chromosome}><td>{row.chromosome}</td><td>{n(row.median_copy_number)}</td><td>{n(row.min_copy_number)}–{n(row.max_copy_number)}</td><td>{row.agreeing_bins}/{row.contributing_bins}</td></tr>)}</tbody></table></div></details>}</Section>
    <Section id="bands" title="Betroffene Zytobanden" note="Annotationen der normalisierten Ereignisse; keine zusätzliche Klassifikationsschwelle aus der Gestaltungsvorlage.">{view.events.some(row=>row.cytobands) ? <div className="bf-table"><table><thead><tr><th>Ereignis</th><th>Zytobanden</th><th>Typ</th></tr></thead><tbody>{view.events.filter(row=>row.cytobands).map(row=><tr key={row.event_id}><td>{row.event_id}</td><td>{row.cytobands}</td><td>{labelType(row.event_type)}</td></tr>)}</tbody></table></div>:<Empty>Keine Zytobandenannotation in den erfassten Ereignissen.</Empty>}</Section>
    <Section id="sv" title="Strukturvarianten" note="Alle vorhandenen Caller-Nachweise bleiben sichtbar. Übereinstimmung mehrerer Caller ist technische Evidenz, keine unabhängige Validierung."><EventTable events={svEvents}/><details><summary>Vollständige SV-Nachweise und Bewertungsgrenzen</summary><Native id="sv-review"/></details></Section>
    <Section id="events" title="Normalisierte Kopienzahlereignisse"><EventTable events={cnvEvents}/><details><summary>Vollständige Ereignisse im statischen Nachweisformat</summary><Native id="events"/></details></Section>
    <Section id="qc" title="Qualitätskontrolle und Readlängen"><Facts rows={qcEntries.map(([key,value])=>[key.replaceAll('_',' '),n(value)])}/><Native id="qc"/></Section>
    <Section id="coverage" title="Abdeckung und Assay-Kontext"><Native id="coverage"/></Section>
    <Section id="methylation" title="Methylierung"><p><Status value={methyl.status}/> · {methyl.reason}</p>{evidence.methylation?<Native id="methylation"/>:<Empty>Keine auswertbaren Methylierungswerte in diesem Bericht. Nicht gemessen ist nicht null.</Empty>}</Section>
    <Section id="marlin" title="Methylierungsklassifikation · MARLIN"><Marlin report={data.marlin}/></Section>
    <Section id="iscn" title="ISCN · Vorschlag zur fachlichen Prüfung" tone="bf-caution"><Native id="iscn"/></Section>
    <Section id="provenance" title="Technische Nachweise"><Facts rows={[["Pipeline",view.pipeline_version],["Git-Stand",view.git_commit],["Darstellung",data.contract],["Referenz",view.reference_id]]}/><details><summary>Werkzeuge, Parameter, Referenzbündel und Prüfsummen</summary><Native id="provenance"/></details></Section>
    <Section id="scope" title="Plattformumfang"><div className="bf-scope">{['COMPLETED','NO_CALL','FAILED','NOT_RUN','NOT_REQUESTED','NOT_RECORDED'].map(status=>{const rows=modules.filter(row=>row.status===status);return rows.length?<div key={status}><Status value={status}/><p>{rows.map(row=>moduleLabels[row.name]).join(' · ')}</p></div>:null;})}</div></Section>
    <Section id="limits" title="Grenzen und Warnungen">{[...new Set([...view.warnings,...(cnv?.warnings||[]),...(cnv?.limitations||[])])].map((text,index)=><p key={index} className="bf-alert">{text}</p>)}<p>Die Modellansichten liefern keine neue Messung und keine klinische Freigabe. Fehlende Daten und fehlgeschlagene Module dürfen nicht als unauffälliger Befund interpretiert werden.</p></Section>
    <Section id="risk" title="Risikogruppe nach Leitlinie" tone="bf-caution"><Empty>Keine automatische Risikozuordnung durch diese Oberfläche. Der Befundvertrag enthält keine freigegebene Leitlinienklassifikation.</Empty></Section>
    <Section id="review" title="Vidierung" tone="bf-blocked"><Facts rows={[["Gespeicherter Prüfstatus",view.release_status],["ISCN-Prüfstatus",data.iscn.review_status],["QC",view.qc_verdict]]}/><p>Dieser Bericht ist ein Forschungsartefakt. Regler und Modellwahl haben keine Wirkung auf Ergebnisse, Prüfsummen oder Freigaben. Eine klinische Freigabe kann hier nicht erteilt werden.</p><button onClick={()=>setShowAudit(!showAudit)} aria-expanded={showAudit}>{showAudit?'Statischen Nachweis schließen':'Statischen Gesamtnachweis anzeigen'}</button></Section>
    <footer className="bf-footer">ONTSeq · {data.contract} · Offline-Bericht · Forschungszwecke<br/>Modellansichten werden nicht gespeichert. Datenbasis und Werkzeugversionen siehe technische Nachweise.</footer>
    </main></div>
  </div>;
}
class ReportBoundary extends React.Component {
  state={failed:false};
  static getDerivedStateFromError(){return {failed:true};}
  componentDidCatch(){document.body.classList.remove('befund-ready');}
  render(){return this.state.failed?<p role="alert">Die interaktive Ansicht konnte nicht geladen werden. Der vollständige statische Befund folgt.</p>:this.props.children;}
}
if(data.contract !== 'befund-interactive-v1') throw new Error('Unknown report contract');
createRoot(document.getElementById('ontseq-befund-root')).render(<ReportBoundary><App/></ReportBoundary>);
