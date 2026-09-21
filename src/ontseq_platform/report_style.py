"""Offline presentation styles for the portable analytical report."""

REPORT_CSS = """
    :root {
      color-scheme: light;
      --ink:#101c26; --muted:#5b6b7a; --line:#dae3ea; --panel:#ffffff;
      --canvas:#e6ecf2; --accent:#0b5580; --accent-soft:#eaf2f7;
      --critical:#8f1d1d; --critical-soft:#fff0f0; --warning:#8a4b08;
      --warning-soft:#fff7e8; --info:#36566f; --info-soft:#eef5f9;
      --ok:#245c45; --ok-soft:#edf7f1; --neutral:#586174;
      --neutral-soft:#f0f2f5;
    }
    * { box-sizing:border-box; }
    body {
      margin:0; font:14px/1.5 Inter,Segoe UI,system-ui,sans-serif; color:var(--ink);
      background:var(--canvas);
    }
    .ruo {
      position:sticky; top:0; z-index:20; background:#731c1c; color:white;
      padding:9px 20px; text-align:center; font-weight:800; letter-spacing:.05em;
    }
    .shell { max-width:1460px; margin:0 auto; padding:24px; }
    .masthead {
      background:var(--panel); border:1px solid var(--line); border-radius:14px;
      padding:24px;
    }
    .masthead h1 { margin:0 0 6px; font-size:28px; }
    .eyebrow {
      color:var(--muted); font-size:11px; font-weight:800; letter-spacing:.08em;
      text-transform:uppercase;
    }
    .identity {
      display:grid; grid-template-columns:repeat(6,minmax(120px,1fr)); gap:10px;
      margin-top:20px;
    }
    .identity div { border-top:2px solid var(--line); padding-top:8px; min-width:0; }
    .identity span {
      display:block; color:var(--muted); font-size:11px; text-transform:uppercase;
    }
    .identity strong { display:block; margin-top:3px; overflow-wrap:anywhere; }
    .layout {
      display:grid; grid-template-columns:220px minmax(0,1fr); gap:18px;
      margin-top:18px;
    }
    nav {
      align-self:start; position:sticky; top:58px; background:var(--panel);
      border:1px solid var(--line); border-radius:12px; padding:10px;
    }
    nav a {
      display:block; padding:9px 10px; border-radius:8px; color:var(--ink);
      text-decoration:none;
    }
    nav a:hover, nav a:focus-visible { background:var(--accent-soft); outline:none; }
    main { min-width:0; }
    section {
      background:var(--panel); border:1px solid var(--line); border-radius:12px;
      padding:20px; margin-bottom:16px;
    }
    h2 { margin:0 0 14px; font-size:20px; }
    h3 { margin:2px 0 0; font-size:17px; }
    h4 { margin:18px 0 8px; font-size:14px; }
    p { margin:6px 0; }
    .muted { color:var(--muted); }
    .module-strip {
      display:grid; grid-template-columns:repeat(auto-fit,minmax(125px,1fr));
      gap:8px; margin-top:12px;
    }
    .module-state {
      border:1px solid var(--line); border-left-width:5px; border-radius:9px;
      padding:10px;
    }
    .module-state span {
      display:block; color:var(--muted); font-size:11px; text-transform:uppercase;
    }
    .module-state strong { display:block; margin-top:3px; }
    .state-completed { border-left-color:var(--ok); background:var(--ok-soft); }
    .state-no-call { border-left-color:var(--warning); background:var(--warning-soft); }
    .state-failed { border-left-color:var(--critical); background:var(--critical-soft); }
    .state-not-run { border-left-color:var(--neutral); background:var(--neutral-soft); }
    .state-label {
      display:inline-block; padding:3px 7px; border-radius:999px;
      border:1px solid currentColor; font-size:11px; font-weight:800;
    }
    .alert {
      border-left:5px solid; padding:12px 14px; margin:10px 0; border-radius:8px;
    }
    .alert-critical { color:var(--critical); background:var(--critical-soft); }
    .alert-warning { color:var(--warning); background:var(--warning-soft); }
    .alert-info { color:var(--info); background:var(--info-soft); }
    .table-wrap { overflow-x:auto; margin-top:10px; }
    table { width:100%; border-collapse:collapse; min-width:620px; }
    caption { text-align:left; font-weight:800; margin:0 0 8px; }
    th,td {
      padding:9px 10px; border-bottom:1px solid var(--line); text-align:left;
      vertical-align:top;
    }
    th { background:#f7f8fa; font-size:12px; }
    code {
      font:12px/1.45 ui-monospace,SFMono-Regular,Consolas,monospace;
      overflow-wrap:anywhere;
    }
    .event-card {
      border:1px solid var(--line); border-radius:11px; padding:16px;
      margin:14px 0; background:#fcfcfd;
    }
    .event-heading {
      display:flex; justify-content:space-between; gap:14px; align-items:flex-start;
    }
    .reportability {
      max-width:360px; border:1px solid var(--line); border-radius:8px;
      padding:7px 9px; font-size:12px; font-weight:700; background:white;
    }
    .event-grid {
      display:grid; grid-template-columns:repeat(4,minmax(120px,1fr)); gap:10px;
      margin:14px 0;
    }
    .event-grid div { border-top:1px solid var(--line); padding-top:7px; min-width:0; }
    dt { color:var(--muted); font-size:11px; text-transform:uppercase; }
    dd { margin:2px 0 0; overflow-wrap:anywhere; }
    .boundary {
      background:var(--info-soft); border-left:4px solid var(--info); padding:10px 12px;
      border-radius:7px;
    }
    .gate-failure {
      background:var(--critical-soft); color:var(--critical); border-radius:8px;
      padding:10px 12px; margin-top:10px;
    }
    .iscn {
      font:700 19px/1.5 ui-monospace,SFMono-Regular,Consolas,monospace;
      color:var(--accent); overflow-wrap:anywhere;
    }
    .empty-state {
      background:var(--info-soft); border:1px solid #bfd1dd; border-radius:8px;
      padding:14px;
    }
    tr.critical { background:var(--critical-soft); }
    input[type="search"] { padding:8px; max-width:100%; }
    footer { color:var(--muted); font-size:12px; padding:4px 2px 24px; }
    @media (max-width:1000px) {
      .identity { grid-template-columns:repeat(3,minmax(120px,1fr)); }
      .layout { grid-template-columns:1fr; }
      nav { position:static; display:flex; overflow-x:auto; gap:4px; }
      nav a { white-space:nowrap; }
      .event-grid { grid-template-columns:repeat(2,minmax(120px,1fr)); }
    }
    @media (max-width:620px) {
      .shell { padding:12px; }
      .masthead { padding:17px; }
      .identity { grid-template-columns:1fr 1fr; }
      .event-heading { display:block; }
      .reportability { margin-top:9px; max-width:none; }
      .event-grid { grid-template-columns:1fr; }
      section { padding:15px; }
    }
    @media print {
      body { background:white; }
      .ruo { position:static; }
      nav { display:none; }
      .layout { display:block; }
      section,.masthead,.event-card { break-inside:avoid; box-shadow:none; }
    }

    body { font-family:"Segoe UI",system-ui,sans-serif; font-size:14px; }
    .ruo { position:static; background:#04222f; color:#f0c674;
      font-size:11px; padding:7px 20px; letter-spacing:.12em; }
    .masthead { border:0; border-radius:0; padding:32px 40px;
      background:linear-gradient(135deg,#04222f 0%,#062a45 45%,#0b5580 100%);
      color:#fff; }
    .masthead-inner { max-width:1200px; margin:auto; display:grid;
      grid-template-columns:minmax(0,1.65fr) minmax(250px,1fr); gap:40px; align-items:end; }
    .masthead .eyebrow { color:#9dc6dd; letter-spacing:.22em; font-weight:500; }
    .masthead h1 { margin:12px 0; font-size:34px; font-weight:600; line-height:1.15; }
    .masthead p { color:#cfe2ee; max-width:66ch; }
    .masthead .identity { display:grid; grid-template-columns:1fr; gap:5px; margin:0;
      font-family:Consolas,ui-monospace,monospace; font-size:12px; }
    .masthead .identity div { display:grid; grid-template-columns:85px minmax(0,1fr);
      gap:14px; padding:0; border:0; }
    .masthead .identity span { color:#9dc6dd; text-transform:none; }
    .masthead .identity strong { font-weight:400; margin:0; }
    .masthead .release { color:#f0c674; }
    .shell { max-width:1280px; padding:24px 40px 40px; }
    .layout { display:block; margin:0; }
    nav { display:flex; flex-wrap:wrap; gap:4px; position:static;
      background:transparent; border:0; padding:0 0 18px; }
    nav a { padding:5px 9px; font-size:12px; color:var(--accent); }
    section { padding:26px 28px 22px; margin-bottom:20px; scroll-margin-top:16px; }
    h2 { font-size:12px; letter-spacing:.15em; text-transform:uppercase;
      color:var(--muted); margin:0 0 14px; }
    h3 { font-size:16px; margin:18px 0 10px; font-weight:600; }
    h4 { color:var(--muted); }
    #overview { border-left:5px solid var(--accent); }
    .summary-grid { display:grid; grid-template-columns:repeat(3,minmax(0,1fr));
      gap:14px; margin:18px 0; }
    .summary-card { border:1px solid #e4ebf1; border-top:3px solid var(--accent);
      border-radius:9px; padding:15px 17px; min-width:0; }
    .summary-card strong { display:block; font:600 18px/1.4 Consolas,monospace;
      overflow-wrap:anywhere; }
    .summary-value { overflow-wrap:anywhere; font-size:25px; font-weight:600;
      color:var(--accent); margin:7px 0; }
    .summary-card p { font-size:12px; color:var(--muted); overflow-wrap:anywhere; }
    .summary-card a { color:var(--accent); font-size:12px; }
    .summary-context { font-size:12px; color:var(--muted); margin:16px 0; }
    .status-grid { display:grid; grid-template-columns:repeat(4,minmax(0,1fr));
      gap:14px; margin:0 0 20px; }
    .status-card { background:#fff; border:1px solid var(--line); border-left:4px solid;
      border-radius:10px; padding:16px 18px; min-width:0; }
    .status-card span { font-size:11px; letter-spacing:.12em; text-transform:uppercase;
      color:var(--muted); }
    .status-card strong { display:block; margin:5px 0; font-size:20px; overflow-wrap:anywhere; }
    .status-card p { color:var(--muted); font-size:12px; }
    .tone-ok { border-left-color:var(--ok); color:var(--ok); }
    .tone-warning { border-left-color:var(--warning); color:var(--warning); }
    .tone-critical { border-left-color:var(--critical); color:var(--critical); }
    .tone-neutral { border-left-color:var(--neutral); color:var(--neutral); }
    .identity { grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); }
    .identity div { border-top:1px solid #e4ebf1; }
    th { background:#f7fafc; color:var(--muted); font-weight:600; }
    td { font-size:12px; }
    table { min-width:620px; }
    section > .table-wrap { margin:14px 0; }
    .boundary { font-size:12px; }
    #warnings { border-left:5px solid var(--warning); }
    svg,img { max-width:100%; }
    .plot svg { height:auto; }
    .cnv-facts { display:grid; grid-template-columns:repeat(4,minmax(0,1fr));
      gap:14px; padding:16px 18px; margin:12px 0 20px; background:#f7fafc;
      border:1px solid #e7eef3; border-radius:10px; }
    .cnv-facts span { color:var(--muted); font-size:12px; }
    .cnv-facts strong { display:block; font:600 18px/1.5 Consolas,monospace; }
    .plot { margin:18px 0; }
    .plot figcaption { color:var(--muted); font-size:12px; margin-top:8px; }
    .event-heading > div { min-width:0; }
    .event-heading h3 { overflow-wrap:anywhere; }
    @media (max-width:850px) {
      .masthead-inner { grid-template-columns:1fr; gap:24px; }
      .status-grid,.summary-grid,.cnv-facts { grid-template-columns:repeat(2,minmax(0,1fr)); }
      .shell { padding:20px; }
      .masthead { padding:28px 24px; }
      .masthead h1 { font-size:29px; }
      nav { overflow:visible; }
    }
    @media (max-width:520px) {
      .shell { padding:16px 12px; }
      section { padding:20px 16px; }
      .masthead { padding:24px 20px; }
      .masthead h1 { font-size:26px; }
      .summary-grid,.cnv-facts { grid-template-columns:1fr; }
      .status-card { padding:12px; }
      .status-card strong { font-size:16px; }
      .event-grid,.identity { grid-template-columns:1fr; }
      .module-strip { grid-template-columns:1fr 1fr; }
      nav a { white-space:normal; }
    }
    @media print {
      @page { margin:14mm; }
      body { font-size:11px; }
      .shell { max-width:none; padding:12px 0; }
      .masthead { background:#fff; color:#101c26; border-bottom:2px solid #0b5580; padding:12px 0; }
      .masthead h1 { font-size:24px; }
      .masthead p,.masthead .eyebrow,.masthead .identity span,.masthead .release { color:#334b5c; }
      .ruo { background:white; color:#8a4b08; }
      nav { display:none; }
      section { break-inside:auto; padding:16px; }
      h2,h3 { break-after:avoid; }
      .summary-card,.status-card,figure,tr { break-inside:avoid; }
      .table-wrap { overflow:visible; }
      table { min-width:0; table-layout:fixed; }
      td,th { overflow-wrap:anywhere; padding:5px; }
      input[type=search] { display:none; }
    }
"""
