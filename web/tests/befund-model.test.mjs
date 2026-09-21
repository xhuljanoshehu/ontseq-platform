import assert from 'node:assert/strict';
import test from 'node:test';
import { copyAt, ratioAt, solutions, moduleState, dilutedRatio, modelSummary } from '../src/befund/model.js';

test('alternative solutions use the same ratio and never mutate measured values', () => {
  const ratio = ratioAt(9, .5, 4);
  assert.equal(ratio, 11/6);
  assert.ok(Math.abs(copyAt(ratio, .8, 2) - 49/12) < 1e-12);
  assert.equal(copyAt(ratio, 0, 2), null);
  assert.equal(copyAt(null, .8, 2), null);
  assert.equal(copyAt(.1, .1, 2), null); // impossible negative copy number
  assert.equal(ratioAt(0, 1, 2), 0); // zero is measured, not absent
  assert.equal(dilutedRatio(4, 4, 0), 1);
});

test('fit alternatives retain pipeline selection and reject non-finite or invalid inputs', () => {
  const fit = {cellularity:.5, ploidy:4, fit_error:0,
    alternatives:[{cellularity:.8, ploidy:2, fit_error:.02},
      {cellularity:0, ploidy:2, fit_error:0}, {cellularity:.5, ploidy:NaN, fit_error:0}]};
  const copy = JSON.stringify(fit);
  assert.equal(solutions(fit).length, 2);
  assert.equal(solutions(fit)[0].pipeline, true);
  assert.equal(JSON.stringify(fit), copy);
});

test('missing outcomes distinguish not requested from requested but unrecorded', () => {
  assert.equal(moduleState('methylation', [], []).status, 'NOT_REQUESTED');
  assert.equal(moduleState('methylation', [], ['methylation']).status, 'NOT_RECORDED');
  assert.equal(moduleState('methylation', [{name:'methylation',status:'FAILED',reason:'guard'}], ['methylation']).status, 'FAILED');
});

test('unavailable and partial chromosome models cannot imply absence of deviations', () => {
  assert.match(modelSummary([], 2), /Keine auswertbaren/);
  assert.match(modelSummary([{chromosome:'chr1',copies:null}], 2), /Keine auswertbaren/);
  assert.match(modelSummary([{chromosome:'chr1',copies:2},{chromosome:'chr2',copies:null}], 2), /1.*nicht bestimmbar/);
});
