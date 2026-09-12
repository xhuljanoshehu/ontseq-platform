#!/usr/bin/env node
/** Build the local, authenticated workspace without a CDN or any demo dependency.
 *
 * Normal build: cd web && npm ci && npm run build
 * Offline bootstrap: set ONTSEQ_WEB_NODE_MODULES to an existing node_modules directory
 * containing the exact lockfile versions, then run this script with Node.
 */
import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { mkdir, readFile, realpath, writeFile } from "node:fs/promises";
import { dirname, join, relative, resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const repository = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const webRoot = join(repository, "web");
const sourceRoot = join(webRoot, "src");
const outputPath = join(repository, "src", "ontseq_platform", "service", "workspace.html");
const modulesRoot = process.env.ONTSEQ_WEB_NODE_MODULES
  ? resolve(process.env.ONTSEQ_WEB_NODE_MODULES)
  : join(webRoot, "node_modules");
const packageJson = JSON.parse(await readFile(join(webRoot, "package.json"), "utf8"));
const lockfile = JSON.parse(await readFile(join(webRoot, "package-lock.json"), "utf8"));
const project = await readFile(join(repository, "pyproject.toml"), "utf8");
const productVersion = /^version\s*=\s*"([^"]+)"/m.exec(project)?.[1];
assert.equal(packageJson.version, productVersion, "Web and core product versions must match");
assert.match(productVersion, /^\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?$/);
assert.equal(lockfile.packages[""].version, productVersion);

// A pre-existing dependency directory is not a version bypass. It must satisfy the lock.
const shippedPackages = ["react", "react-dom", "scheduler"];
const runtimeRoots = await Promise.all(shippedPackages.map((name) => realpath(join(modulesRoot, name))));
for (const name of [...shippedPackages, "esbuild"]) {
  const installed = JSON.parse(await readFile(join(modulesRoot, name, "package.json"), "utf8"));
  assert.equal(installed.version, lockfile.packages[`node_modules/${name}`].version, `${name} differs from package-lock.json`);
}
const { build } = await import(pathToFileURL(join(modulesRoot, "esbuild", "lib", "main.js")).href);
const output = await build({
  absWorkingDir: webRoot,
  entryPoints: [join(sourceRoot, "main.jsx")],
  bundle: true,
  write: false,
  outdir: join(webRoot, ".build-unused"),
  format: "iife",
  platform: "browser",
  target: ["chrome100", "edge100", "firefox100", "safari15.4"],
  jsx: "automatic",
  tsconfigRaw: { compilerOptions: {} },
  define: { "process.env.NODE_ENV": '"production"' },
  minify: true,
  sourcemap: false,
  legalComments: "inline",
  charset: "utf8",
  alias: Object.fromEntries(shippedPackages.map((name) => [name, join(modulesRoot, name)])),
  nodePaths: [modulesRoot],
  metafile: true,
  logLevel: "warning",
});

// Restrict the final module graph to this source tree and the three runtime libraries.
const within = (path, root) => { const part = relative(root, path); return part === "" || (!part.startsWith("..") && !/^[A-Za-z]:/.test(part)); };
for (const input of Object.keys(output.metafile.inputs)) {
  const absolute = resolve(webRoot, input);
  assert.ok(within(absolute, sourceRoot) || runtimeRoots.some((root) => within(absolute, root)), `Unexpected live module: ${input}`);
}
const script = output.outputFiles.find((file) => file.path.endsWith(".js"))?.text;
const css = output.outputFiles.find((file) => file.path.endsWith(".css"))?.text;
assert.ok(script && css, "Expected one JS bundle and one CSS bundle");
assert.equal(output.outputFiles.length, 2, "Live runtime must have no external assets");
assert.equal(script.includes("__ONTSEQ_TOKEN__"), false, "Only HTML metadata may carry the replaceable token marker");
assert.doesNotMatch(css, /@import|url\s*\(/i, "Live CSS must not reference network assets");

const licenses = await Promise.all(shippedPackages.map(async (name) => {
  const text = await readFile(join(modulesRoot, name, "LICENSE"), "utf8");
  return `${name}@${lockfile.packages[`node_modules/${name}`].version}\n${"=".repeat(64)}\n${text.trim()}\n`;
}));
const escapeScript = (value) => value.replace(/<\/script/gi, "<\\/script");
const escapeStyle = (value) => value.replace(/<\/style/gi, "<\\/style");
const html = `<!doctype html>
<html lang="de">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <meta name="ontseq-session-token" content="__ONTSEQ_TOKEN__" />
    <meta name="ontseq-ui-version" content="${productVersion}" />
    <title>ONTSeq ${productVersion} · Lokaler Live-Arbeitsplatz</title>
    <style>${escapeStyle(css)}</style>
  </head>
  <body>
    <div id="root"></div>
    <noscript>JavaScript ist für den lokalen Arbeitsplatz erforderlich. Ohne JavaScript wird keine Analyse gestartet.</noscript>
    <script>${escapeScript(script)}</script>
    <script type="text/plain" id="ontseq-third-party-notices">${escapeScript(licenses.join("\n"))}</script>
  </body>
</html>
`;
assert.equal(html.split("__ONTSEQ_TOKEN__").length - 1, 1);
assert.doesNotMatch(html, /<(?:script|link)\b[^>]*(?:src|href)\s*=/i, "HTML must be self-contained");

const digest = createHash("sha256").update(html).digest("hex");
if (process.argv.includes("--check")) {
  assert.equal(await readFile(outputPath, "utf8"), html, "workspace.html is stale; run npm run build in web/");
  process.stdout.write(`Live workspace is reproducible: SHA256 ${digest}\n`);
} else {
  await mkdir(dirname(outputPath), { recursive: true });
  await writeFile(outputPath, html, "utf8");
  process.stdout.write(`Built ${relative(repository, outputPath)} (${Buffer.byteLength(html)} bytes)\nSHA256 ${digest}\n`);
}
