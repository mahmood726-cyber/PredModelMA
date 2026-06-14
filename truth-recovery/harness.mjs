/**
 * Truth-recovery harness for PredModelMA.
 * Pools k studies around a known true C with the engine's OWN poolCstat (z-interval,
 * logit transform) and knappHartungC (HKSJ, t_{k-1}). Measures coverage of the TRUE C,
 * bias, and CI width across heterogeneity. Also checks O:E recovery via poolOE.
 */
import { poolCstat, knappHartungC, poolOE, __setTransform, __setAlpha } from './engine.mjs';
import { generateCstatStudies, generateOEStudies } from './dgp-cstat.mjs';

__setAlpha(0.05);
__setTransform('logit');

export function runCell({Ctrue,tau,k,nReps,baseSeed}){
  let covZ=0,covH=0,biasSum=0,wZ=0,wH=0,ok=0;
  for(let r=0;r<nReps;r++){
    const {studies}=generateCstatStudies({Ctrue,tau,k,nEvents:200,nControls:400,seed:baseSeed+r});
    const cp=poolCstat(studies,'reml');
    if(!cp||!isFinite(cp.c)||!isFinite(cp.ci[0])||!isFinite(cp.ci[1]))continue;
    const kh=knappHartungC(cp);
    ok++;
    if(Ctrue>=cp.ci[0]&&Ctrue<=cp.ci[1])covZ++;
    if(kh&&Ctrue>=kh.ci[0]&&Ctrue<=kh.ci[1])covH++;
    biasSum+=(cp.c-Ctrue);
    wZ+=(cp.ci[1]-cp.ci[0]);
    if(kh)wH+=(kh.ci[1]-kh.ci[0]);
  }
  return {tau,k,n:ok,coverageZ:covZ/ok,coverageHKSJ:covH/ok,biasC:biasSum/ok,widthZ:wZ/ok,widthHKSJ:wH/ok};
}

export function runOECell({OEtrue,tau,k,nReps,baseSeed}){
  let cov=0,bias=0,ok=0;
  for(let r=0;r<nReps;r++){
    const {studies}=generateOEStudies({OEtrue,tau,k,expectedEvents:100,seed:baseSeed+r});
    const op=poolOE(studies,'reml');
    if(!op||!isFinite(op.oe))continue;
    ok++;
    if(OEtrue>=op.ci[0]&&OEtrue<=op.ci[1])cov++;
    bias+=(op.oe-OEtrue);
  }
  return {tau,k,n:ok,coverage:cov/ok,bias:bias/ok};
}

export function runSweep({Ctrue=0.75,nReps=2000,baseSeed=2000}={}){
  const cells=[];
  for(const k of [5,10,20]) for(const tau of [0.0,0.15,0.35])
    cells.push(runCell({Ctrue,tau,k,nReps,baseSeed:baseSeed+k*100000+Math.round(tau*1000)*137}));
  return {Ctrue,nReps,cells};
}

if(process.argv[1]?.endsWith('harness.mjs')){
  const out=runSweep();
  console.log(`PredModelMA known-truth C recovery | C_true=${out.Ctrue} | reps=${out.nReps}/cell`);
  console.log('tau    k   covZ   covHKSJ  biasC      widthZ   widthHKSJ');
  for(const c of out.cells)
    console.log(`${c.tau.toFixed(2)}  ${String(c.k).padStart(2)}  ${(c.coverageZ*100).toFixed(1)}%  ${(c.coverageHKSJ*100).toFixed(1)}%   ${c.biasC>=0?'+':''}${c.biasC.toFixed(5)}  ${c.widthZ.toFixed(4)}  ${c.widthHKSJ.toFixed(4)}`);
  console.log('\nO:E recovery (O:E_true=1.0, k=10):');
  for(const tau of [0.0,0.2]){const o=runOECell({OEtrue:1.0,tau,k:10,nReps:2000,baseSeed:9000});console.log(`tau=${tau} coverage=${(o.coverage*100).toFixed(1)}% bias=${o.bias>=0?'+':''}${o.bias.toFixed(5)}`);}
}
