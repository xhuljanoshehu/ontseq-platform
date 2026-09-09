import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

const html = await readFile(new URL("../../src/ontseq_platform/service/workspace.html", import.meta.url), "utf8").catch(() => null);
const packageJson = JSON.parse(await readFile(new URL("../package.json", import.meta.url), "utf8"));
const lock = JSON.parse(await readFile(new URL("../package-lock.json", import.meta.url), "utf8"));

test("web dependencies are exact and the integrity-pinned graph contains no demo libraries", () => {
  assert.equal(packageJson.private, true);
  assert.deepEqual(packageJson.dependencies, { react: "19.2.0", "react-dom": "19.2.0" });
  assert.deepEqual(packageJson.devDependencies, { esbuild: "0.25.12" });
  assert.equal(lock.lockfileVersion, 3);
  assert.equal(lock.packages[""].version, packageJson.version);
  for (const [name, value] of Object.entries(lock.packages)) {
    if (!name) continue;
    assert.match(name, /^node_modules\/(?:react|react-dom|scheduler|esbuild|@esbuild\/[^/]+)$/);
    assert.match(value.integrity, /^sha512-[A-Za-z0-9+/]+=*$/);
  }
  for (const [name, version] of Object.entries(lock.packages["node_modules/esbuild"].optionalDependencies)) {
    assert.equal(lock.packages[`node_modules/${name}`].version, version);
    assert.equal(lock.packages[`node_modules/${name}`].optional, true);
  }
});

test("packaged workspace has one unissued token marker and no external assets", () => {
  assert.ok(html, "Run npm run build before the build-contract test");
  assert.equal(html.startsWith("<!doctype html>"), true);
  assert.equal(html.split("__ONTSEQ_TOKEN__").length - 1, 1);
  assert.match(html, /<meta name="ontseq-session-token" content="__ONTSEQ_TOKEN__"/);
  assert.doesNotMatch(html, /<(?:script|link)\b[^>]*(?:src|href)\s*=/i);
  assert.doesNotMatch(html, /<script\b[^>]*type=["']module/);
  assert.match(html, /Nur für Forschungszwecke/);
  assert.match(html, /Keine elektronische Signatur/);
});

test("runtime license notices and the product version are bundled", () => {
  assert.ok(html);
  assert.ok(html.includes(`<meta name="ontseq-ui-version" content="${packageJson.version}"`));
  assert.match(html, /id="ontseq-third-party-notices"/);
  for (const name of ["react", "react-dom", "scheduler"]) {
    assert.ok(html.includes(`${name}@${lock.packages[`node_modules/${name}`].version}`));
  }
  assert.match(html, /Permission is hereby granted, free of charge/);
});
