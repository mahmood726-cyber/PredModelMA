/**
 * Truth-recovery assertions for PredModelMA.
 * Run: node truth-recovery/test-truth-recovery.mjs   (exit 0 pass, 1 fail)
 */
import { poolCstat, knappHartungC, poolOE, logit, expit, __setTransform, __setAlpha } from './engine.mjs';
import { generateCstatStudies, generateOEStudies } from './dgp-cstat.mjs';
import { runCell, runOECell } from './harness.mjs';

__setAlpha(0.05); __setTransform('logit');
let failures=0;
function ok(name,cond,detail=''){console.log(`${cond?'PASS':'FAIL'}  ${name}${detail?'  ['+detail+']':''}`);if(!cond)failures++;}

// 1. logit/expit round-trip (engine's own transform pair).
{ const C=0.73; ok('logit/expit round-trip', Math.abs(expit(logit(C))-C)<1e-12); }

// 2. Point estimate essentially unbiased at moderate heterogeneity.
{ const cell=runCell({Ctrue:0.75,tau:0.15,k:10,nReps:2000,baseSeed:42});
  ok('C-stat bias < 0.01', Math.abs(cell.biasC)<0.01, `bias=${cell.biasC.toFixed(5)}`); }

// 3. Pooled C in (0.5,1) and z-CI brackets the point estimate.
{ const {studies}=generateCstatStudies({Ctrue:0.78,tau:0.2,k:8,seed:7});
  const cp=poolCstat(studies,'reml');
  ok('pooled C in (0.5,1) & CI brackets it',
     cp.c>0.5&&cp.c<1&&cp.ci[0]<cp.c&&cp.c<cp.ci[1],
     `C=${cp.c.toFixed(4)} CI=[${cp.ci[0].toFixed(4)},${cp.ci[1].toFixed(4)}]`); }

// 4. Documented: default z-interval under-covers at k=5 high tau; HKSJ recovers.
{ const cell=runCell({Ctrue:0.75,tau:0.35,k:5,nReps:2000,baseSeed:99});
  ok('z-interval under-covers at k=5 high tau (documented)',
     cell.coverageZ<0.93, `covZ=${(cell.coverageZ*100).toFixed(1)}%`);
  ok('HKSJ (knappHartungC) coverage > z at k=5 high tau',
     cell.coverageHKSJ>cell.coverageZ, `covH=${(cell.coverageHKSJ*100).toFixed(1)}%`); }

// 5. O:E pooled on log scale recovers true O:E with small bias.
{ const cell=runOECell({OEtrue:1.0,tau:0.2,k:10,nReps:1500,baseSeed:3000});
  ok('O:E pooled (log scale) |bias| < 0.05', Math.abs(cell.bias)<0.05, `bias=${cell.bias.toFixed(5)}`); }

console.log(`\n${failures===0?'ALL TESTS PASSED':failures+' TEST(S) FAILED'}`);
process.exit(failures===0?0:1);
