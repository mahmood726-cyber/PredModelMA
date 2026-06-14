/**
 * Standalone seeded known-truth DGP for C-statistic meta-analysis (PredModelMA).
 *
 * Self-contained seeded PRNG (splitmix32 + Box-Muller) so the DGP does NOT depend
 * on the engine internals (the inline engine has no seedable rnorm).
 *
 * Model (Debray 2017 / Riley 2016):
 *   - true population C-statistic C_true in (0.5,1).
 *   - study i true logit-C: theta_i = logit(C_true) + N(0, tau).
 *   - observed logit-C = theta_i + N(0, se_logit_i^2), se_logit from a finite
 *     validation sample (Hanley-McNeil natural SE -> logit via delta method).
 * Studies report C (natural) and its CI [c_lo,c_hi] so the engine's own parseData
 * path (se_c = (chi-clo)/(2*1.96)) is exercised.
 */
import { logit, expit } from './engine.mjs';

function splitmix32(a){return function(){a|=0;a=a+0x9e3779b9|0;let t=a^a>>>16;t=Math.imul(t,0x21f0aaad);t=t^t>>>15;t=Math.imul(t,0x735a2d97);return (t^t>>>15)>>>0;};}
export function makeRng(seed){
  const sm=splitmix32(seed);let s=[sm(),sm(),sm(),sm()];
  function u(){let r=(s[0]+s[3])>>>0|0;s[0]=r;const t=s[1]<<9;s[2]^=s[0];s[3]^=s[1];s[1]^=s[2];s[0]^=s[3];s[2]^=t;s[3]=(s[3]<<11|s[3]>>>21)>>>0;return (r>>>0)/4294967296;}
  return {u, rnorm(mu=0,sd=1){let u1=u();while(u1===0)u1=u();const u2=u();return mu+sd*Math.sqrt(-2*Math.log(u1))*Math.cos(2*Math.PI*u2);}};
}

/** Hanley-McNeil (1982) SE of the C-statistic on the natural scale. */
export function seCstatNatural(C,n1,n0){
  const Q1=C/(2-C), Q2=(2*C*C)/(1+C);
  const v=(C*(1-C)+(n1-1)*(Q1-C*C)+(n0-1)*(Q2-C*C))/(n1*n0);
  return Math.sqrt(Math.max(v,1e-12));
}

export function generateCstatStudies({Ctrue,tau,k,nEvents=200,nControls=400,seed}){
  const rng=makeRng(seed);
  const logitTrue=logit(Ctrue);
  const studies=[];
  for(let i=0;i<k;i++){
    const theta_i=logitTrue+(tau>0?rng.rnorm(0,tau):0);
    const Ci=expit(theta_i);
    const seC=seCstatNatural(Ci,nEvents,nControls);
    const seLogit=seC/(Ci*(1-Ci));
    let Cobs=expit(theta_i+rng.rnorm(0,seLogit));
    Cobs=Math.min(0.9999,Math.max(0.5001,Cobs));
    const seObs=seCstatNatural(Cobs,nEvents,nControls);
    studies.push({study:'S'+i, c:Cobs, se_c:seObs,
      clo:Cobs-1.96*seObs, chi:Cobs+1.96*seObs, n:nEvents+nControls, events:nEvents});
  }
  return {studies, Ctrue, tau};
}

export function generateOEStudies({OEtrue,tau,k,expectedEvents=100,seed}){
  const rng=makeRng(seed);
  const logTrue=Math.log(OEtrue);
  const studies=[];
  for(let i=0;i<k;i++){
    const theta_i=logTrue+(tau>0?rng.rnorm(0,tau):0);
    const OE_i=Math.exp(theta_i);
    const E=expectedEvents, meanO=E*OE_i;
    let O=Math.max(1,Math.round(meanO+rng.rnorm(0,Math.sqrt(meanO))));
    const oe=O/E, se_oe=oe/Math.sqrt(O); // natural-scale SE; engine uses log-CI path
    studies.push({study:'S'+i, oe, se_oe,
      oelo:Math.exp(Math.log(oe)-1.96/Math.sqrt(O)), oehi:Math.exp(Math.log(oe)+1.96/Math.sqrt(O))});
  }
  return {studies, OEtrue, tau};
}
