/**
 * Truth-recovery engine adapter for PredModelMA (ADDITIVE, non-invasive).
 *
 * The real engine lives inline in ../index.html as plain <script> functions.
 * We extract the PURE pooling/transform/HKSJ functions VERBATIM (no edits to
 * their bodies) below, stub the two `document.getElementById(...)` reads they
 * rely on (only 'alpha' and 'cTransform' are referenced by these functions),
 * and re-export. Math identical to what the browser runs.
 *
 * --- BEGIN verbatim extract from index.html lines 398-491, 494-654, 960-978 ---
 */

// Minimal DOM stub: poolCstatTransform/poolOE/knappHartungC read #alpha and #cTransform.
const __els = { alpha: { value: '0.05' }, cTransform: { value: 'logit' } };
globalThis.document = globalThis.document || {
  getElementById(id) { return __els[id] || { value: '' }; }
};
export function __setAlpha(a){ __els.alpha.value = String(a); }
export function __setTransform(t){ __els.cTransform.value = String(t); }

function normalCDF(z){const a1=0.254829592,a2=-0.284496736,a3=1.421413741,a4=-1.453152027,a5=1.061405429,p=0.3275911;const s=z<0?-1:1;const x=Math.abs(z)/Math.sqrt(2);const t=1/(1+p*x);const y=1-(((((a5*t+a4)*t)+a3)*t+a2)*t+a1)*t*Math.exp(-x*x);return 0.5*(1+s*y)}
function normalQuantile(p){if(p<=0)return-Infinity;if(p>=1)return Infinity;if(p===0.5)return 0;const a=p<0.5?p:1-p;const t=Math.sqrt(-2*Math.log(a));const c0=2.515517,c1=0.802853,c2=0.010328,d1=1.432788,d2=0.189269,d3=0.001308;let z=t-(c0+c1*t+c2*t*t)/(1+d1*t+d2*t*t+d3*t*t*t);if(p<0.5)z=-z;return z}
function chiSqP(x,df){if(x<=0)return 1;return 1-regGammaP(df/2,x/2)}
function regGammaP(a,x){if(x<=0)return 0;const lga=lgamma(a);if(x<a+1){let s=1/a,t=1/a;for(let n=1;n<200;n++){t*=x/(a+n);s+=t;if(Math.abs(t)<1e-12*Math.abs(s))break}return s*Math.exp(-x+a*Math.log(x)-lga)}else{let f=1,b2=x+1-a,c2=1/1e-30,d2=1/b2;f=d2;for(let i=1;i<=200;i++){const an=-i*(i-a),bn=x+2*i+1-a;d2=bn+an*d2;if(Math.abs(d2)<1e-30)d2=1e-30;d2=1/d2;c2=bn+an/c2;if(Math.abs(c2)<1e-30)c2=1e-30;f*=d2*c2;if(Math.abs(d2*c2-1)<1e-10)break}return 1-f*Math.exp(-x+a*Math.log(x)-lga)}}
function lgamma(x){const g=7,c=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];if(x<0.5)return Math.log(Math.PI/Math.sin(Math.PI*x))-lgamma(1-x);x-=1;let a=c[0];const t=x+g+0.5;for(let i=1;i<g+2;i++)a+=c[i]/(x+i);return 0.5*Math.log(2*Math.PI)+(x+0.5)*Math.log(t)-t+Math.log(a)}

function logit(p){return Math.log(p/(1-p))}
function expit(x){if(x>500)return 1;if(x<-500)return 0;return 1/(1+Math.exp(-x))}

// t-distribution quantile via Hill's algorithm (Abramowitz & Stegun approx for large df, iterative for small)
function tQuantile(p,df){
  if(df<=0)return NaN;
  if(p<=0)return-Infinity;if(p>=1)return Infinity;if(p===0.5)return 0;
  // For df>=30 use normal approx with correction
  if(df>=30){
    const z=normalQuantile(p);
    const g1=(z*z*z+z)/(4*df);
    const g2=(5*z*z*z*z*z+16*z*z*z+3*z)/(96*df*df);
    return z+g1+g2;
  }
  // Newton-Raphson on Student t CDF
  let t=normalQuantile(p);
  for(let iter=0;iter<50;iter++){
    const cdfVal=tCDF(t,df);
    const pdfVal=tPDF(t,df);
    if(pdfVal<1e-15)break;
    const delta=(cdfVal-p)/pdfVal;
    t-=delta;
    if(Math.abs(delta)<1e-10)break;
  }
  return t;
}
function tPDF(t,df){
  return Math.exp(lgamma((df+1)/2)-lgamma(df/2)-0.5*Math.log(df*Math.PI)-(df+1)/2*Math.log(1+t*t/df));
}
function tCDF(t,df){
  // Use regularized incomplete beta function
  const x=df/(df+t*t);
  const ib=regBetaI(x,df/2,0.5);
  return t>=0?(1-0.5*ib):(0.5*ib);
}
function regBetaI(x,a,b){
  // Regularized incomplete beta via continued fraction (Lentz)
  if(x<=0)return 0;if(x>=1)return 1;
  const lbeta=lgamma(a)+lgamma(b)-lgamma(a+b);
  const front=Math.exp(a*Math.log(x)+b*Math.log(1-x)-lbeta)/a;
  if(x<(a+1)/(a+b+2)){return front*betaCF(x,a,b)}
  else{return 1-Math.exp(a*Math.log(x)+b*Math.log(1-x)-lbeta)/b*betaCF(1-x,b,a)}
}
function betaCF(x,a,b){
  const maxIter=200;const eps=1e-12;
  let qab=a+b,qap=a+1,qam=a-1;
  let c=1,d=1-qab*x/qap;
  if(Math.abs(d)<1e-30)d=1e-30;d=1/d;let h=d;
  for(let m=1;m<=maxIter;m++){
    let m2=2*m;
    let aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1+aa*d;if(Math.abs(d)<1e-30)d=1e-30;d=1/d;
    c=1+aa/c;if(Math.abs(c)<1e-30)c=1e-30;
    h*=d*c;
    aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1+aa*d;if(Math.abs(d)<1e-30)d=1e-30;d=1/d;
    c=1+aa/c;if(Math.abs(c)<1e-30)c=1e-30;
    const del=d*c;h*=del;
    if(Math.abs(del-1)<eps)break;
  }
  return h;
}

// C-statistic transformations (Debray et al. 2017, Stat Med)
function cTransformLogit(c,se_c){
  const lc=logit(c);
  const se_lc=se_c/(c*(1-c));// delta method
  return{theta:lc,se:se_lc};
}
function cTransformLoglog(c,se_c){
  // g(c) = log(-log(c)); g'(c) = -1/(c*log(c))
  const clamp=Math.min(Math.max(c,0.501),0.9999);
  const theta=Math.log(-Math.log(clamp));
  const deriv=-1/(clamp*Math.log(clamp));
  const se=Math.abs(deriv)*se_c;
  return{theta,se};
}
function cTransformArcsine(c,se_c){
  // g(c) = asin(sqrt(c)); g'(c) = 1/(2*sqrt(c*(1-c)))
  const clamp=Math.min(Math.max(c,0.501),0.9999);
  const theta=Math.asin(Math.sqrt(clamp));
  const deriv=1/(2*Math.sqrt(clamp*(1-clamp)));
  const se=Math.abs(deriv)*se_c;
  return{theta,se};
}
function cBackTransformLogit(theta){return expit(theta)}
function cBackTransformLoglog(theta){return Math.exp(-Math.exp(theta))}
function cBackTransformArcsine(theta){const s=Math.sin(theta);return s*s}

function poolCstatTransform(data,method,transform){
  const valid=data.filter(d=>d.c>0.5&&d.c<1&&d.se_c>0);
  if(valid.length<2)return null;
  const tfn=transform==='loglog'?cTransformLoglog:transform==='arcsine'?cTransformArcsine:cTransformLogit;
  const backfn=transform==='loglog'?cBackTransformLoglog:transform==='arcsine'?cBackTransformArcsine:cBackTransformLogit;
  const transformed=valid.map(d=>{
    const{theta,se}=tfn(d.c,d.se_c);
    return{...d,lc:theta,se_lc:se,v:se*se};
  });
  const pooled=rePool(transformed.map(t=>t.lc),transformed.map(t=>t.v),method);
  const alpha=parseFloat(document.getElementById('alpha').value)||0.05;
  const zc=normalQuantile(1-alpha/2);
  const c_pooled=backfn(pooled.beta);
  const c_lo=backfn(pooled.beta-zc*pooled.se);
  const c_hi=backfn(pooled.beta+zc*pooled.se);
  // For loglog, back-transform reverses direction so swap if needed
  const ci=[Math.min(c_lo,c_hi),Math.max(c_lo,c_hi)];
  // Prediction interval (IntHout et al. 2016): t_{k-2} * sqrt(se^2 + tau^2)
  let pi=null;
  if(pooled.k>=3){
    const tVal=tQuantile(1-alpha/2,pooled.k-2);
    const pi_se=Math.sqrt(pooled.se*pooled.se+pooled.tau2);
    const pi_lo=backfn(pooled.beta-tVal*pi_se);
    const pi_hi=backfn(pooled.beta+tVal*pi_se);
    pi=[Math.min(pi_lo,pi_hi),Math.max(pi_lo,pi_hi)];
  }
  return{...pooled,c:c_pooled,ci,pi,transform,studies:transformed};
}

function parseData(){
  const raw=document.getElementById('dataInput').value.trim();
  if(!raw)return null;
  const lines=raw.split('\n').filter(l=>l.trim());
  const data=[];
  for(const line of lines){
    const p=line.trim().split(/[\t,]+/);
    if(p.length<6)continue;
    const study=p[0].trim();
    const n=parseInt(p[1]),events=parseInt(p[2]);
    const c=parseFloat(p[3]),clo=parseFloat(p[4]),chi=parseFloat(p[5]);
    let oe=null,oelo=null,oehi=null,calslope=null,probast='U';
    if(p.length>=9){oe=parseFloat(p[6]);oelo=parseFloat(p[7]);oehi=parseFloat(p[8])}
    if(p.length>=10)calslope=parseFloat(p[9]);
    if(p.length>=11)probast=p[10].trim().toUpperCase().charAt(0);
    // Optional extended columns: CalSlope_SE, NRI_events, NRI_nonevents, NRI_event_SE, NRI_nonevent_SE, Applicability
    let calslope_se=null, nri_event=null, nri_nonevent=null, nri_event_se=null, nri_nonevent_se=null, applicability='U';
    if(p.length>=12) calslope_se=parseFloat(p[11]);
    if(p.length>=13) nri_event=parseFloat(p[12]);
    if(p.length>=14) nri_nonevent=parseFloat(p[13]);
    if(p.length>=15) nri_event_se=parseFloat(p[14]);
    if(p.length>=16) nri_nonevent_se=parseFloat(p[15]);
    if(p.length>=17) applicability=p[16].trim().toUpperCase().charAt(0);
    // Extended columns: IDI, IDI_SE, CITL, CITL_SE, TP, FP, Threshold(%)
    let idi=null,idi_se=null,citl=null,citl_se=null,tp=null,fp=null,nb_threshold=null;
    if(p.length>=18) idi=parseFloat(p[17]);
    if(p.length>=19) idi_se=parseFloat(p[18]);
    if(p.length>=20) citl=parseFloat(p[19]);
    if(p.length>=21) citl_se=parseFloat(p[20]);
    if(p.length>=22) tp=parseFloat(p[21]);
    if(p.length>=23) fp=parseFloat(p[22]);
    if(p.length>=24) nb_threshold=parseFloat(p[23]);
    if(isNaN(n)||isNaN(c))continue;
    // Derive SE from CI if available
    let se_c=null;
    if(isFinite(clo)&&isFinite(chi)){se_c=(chi-clo)/(2*1.96)}
    // Newcombe variance approximation when CI is not available (Newcombe 2006)
    if((se_c===null||!isFinite(se_c)||se_c<=0)&&isFinite(c)&&c>0.5&&c<1&&events>0&&n>events){
      const nEvents=events;
      const nNon=n-events;
      const varC=c*(1-c)*(1+(nEvents/2-1)*(1-c)/(2-c)+(nNon/2-1)*c/(1+c))/(nEvents*nNon);
      se_c=Math.sqrt(varC);
    }
    let se_oe=null;
    if(isFinite(oelo)&&isFinite(oehi)&&oe>0){
      se_oe=(Math.log(oehi)-Math.log(oelo))/(2*1.96);
    }
    // Derive calslope SE if not provided: assume SE ~ 0.1 * |slope| as rough default
    if(isFinite(calslope)&&!isFinite(calslope_se)){calslope_se=null}
    data.push({study,n,events,c,clo,chi,se_c,oe,oelo,oehi,se_oe,calslope,calslope_se,probast,
      probast_d1:probast,probast_d2:probast,probast_d3:probast,probast_d4:probast,
      nri_event,nri_nonevent,nri_event_se,nri_nonevent_se,applicability,
      idi,idi_se,citl,citl_se,tp,fp,nb_threshold});
  }
  return data.length>=2?data:null;
}

// Pool c-statistics — uses selected transformation
function poolCstat(data,method){
  const transform=document.getElementById('cTransform').value||'logit';
  return poolCstatTransform(data,method,transform);
}

// Pool O:E ratios on log scale
function poolOE(data,method){
  const valid=data.filter(d=>d.oe>0&&d.se_oe>0);
  if(valid.length<2)return null;
  const transformed=valid.map(d=>{
    const loe=Math.log(d.oe);
    return{...d,loe,v:d.se_oe*d.se_oe};
  });
  const pooled=rePool(transformed.map(t=>t.loe),transformed.map(t=>t.v),method);
  const alpha=parseFloat(document.getElementById('alpha').value)||0.05;
  const zc=normalQuantile(1-alpha/2);
  const oe_pooled=Math.exp(pooled.beta);
  const oe_lo=Math.exp(pooled.beta-zc*pooled.se);
  const oe_hi=Math.exp(pooled.beta+zc*pooled.se);
  return{...pooled,oe:oe_pooled,ci:[oe_lo,oe_hi],studies:transformed};
}

// Generic random-effects pooling
function rePool(yi,vi,method){
  const k=yi.length;
  // Fixed-effect first
  let sumW=0,sumWY=0;
  for(let i=0;i<k;i++){const w=1/vi[i];sumW+=w;sumWY+=w*yi[i]}
  const fe=sumWY/sumW;
  let Q=0;
  for(let i=0;i<k;i++){const w=1/vi[i];Q+=w*(yi[i]-fe)*(yi[i]-fe)}
  const df=k-1;
  const Qp=chiSqP(Q,df);
  let sumW2=0;
  for(let i=0;i<k;i++)sumW2+=(1/vi[i])*(1/vi[i]);
  const C=sumW-sumW2/sumW;

  let tau2;
  if(method==='reml'){
    // REML via iterative Fisher scoring
    tau2=Math.max(0,(Q-df)/C);
    for(let iter=0;iter<50;iter++){
      let sumW_re=0,sumW2_re=0,sumWY_re=0;
      for(let i=0;i<k;i++){const w=1/(vi[i]+tau2);sumW_re+=w;sumW2_re+=w*w;sumWY_re+=w*yi[i]}
      const mu=sumWY_re/sumW_re;
      let dL=0,d2L=0;
      for(let i=0;i<k;i++){
        const w=1/(vi[i]+tau2);
        dL+=-w+w*w*(yi[i]-mu)*(yi[i]-mu);
        d2L+=w*w-2*w*w*w*(yi[i]-mu)*(yi[i]-mu);
      }
      dL*=0.5;d2L*=0.5;
      // REML adjustment
      dL+=0.5*sumW2_re/sumW_re;
      const step=-dL/d2L;
      tau2=Math.max(0,tau2+step);
      if(Math.abs(step)<1e-8)break;
    }
  }else{
    tau2=Math.max(0,(Q-df)/C);
  }

  let sumWre=0,sumWreY=0;
  const weights=[];
  for(let i=0;i<k;i++){const w=1/(vi[i]+tau2);sumWre+=w;sumWreY+=w*yi[i];weights.push(w)}
  const beta=sumWreY/sumWre;
  const se=1/Math.sqrt(sumWre);
  const z=beta/se;
  const p=2*(1-normalCDF(Math.abs(z)));
  const I2=Math.max(0,(Q-df)/Q)*100;
  const totalW=weights.reduce((s,w)=>s+w,0);
  const wPct=weights.map(w=>(w/totalW*100));
  return{beta,se,z,p,Q,Qp,I2,tau2,k,weights:wPct};
}

function knappHartungC(cPool) {
  if (!cPool || cPool.k < 2) return null;
  const k = cPool.k;
  const alpha = parseFloat(document.getElementById('alpha').value) || 0.05;
  // HKSJ: q_KH = max(1, Q/(k-1)) — floor at 1 per advanced-stats rules
  const qKH = Math.max(1, cPool.Q / (k - 1));
  // Use t_{k-1} NOT z
  const tCrit = tQuantile(1 - alpha / 2, k - 1);
  // KH-adjusted SE on transformed scale
  const seKH = cPool.se * Math.sqrt(qKH);
  // Back-transform CI
  const transform = cPool.transform || 'logit';
  const backfn = transform === 'loglog' ? cBackTransformLoglog : transform === 'arcsine' ? cBackTransformArcsine : cBackTransformLogit;
  const lo = backfn(cPool.beta - tCrit * seKH);
  const hi = backfn(cPool.beta + tCrit * seKH);
  const ci = [Math.min(lo, hi), Math.max(lo, hi)];
  const ciWidth = ci[1] - ci[0];
  return { ci, qKH, tCrit, seKH, ciWidth };
}

/* --- END verbatim extract --- */

export {
  logit, expit, normalCDF, normalQuantile, tQuantile,
  cTransformLogit, cBackTransformLogit,
  poolCstat, poolCstatTransform, poolOE, rePool, knappHartungC,
};
