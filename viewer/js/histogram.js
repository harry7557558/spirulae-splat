// Tiny canvas bar-chart for parameter distributions. Dependency-free.

export function drawHistogram(canvas, hist, label) {
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth, h = canvas.clientHeight;
  canvas.width = w*dpr; canvas.height = h*dpr;
  const ctx = canvas.getContext('2d');
  ctx.scale(dpr, dpr);
  ctx.clearRect(0,0,w,h);

  const { bins, min, max } = hist;
  const n = bins.length;
  let peak = 0; for (const b of bins) if (b>peak) peak = b;
  if (peak <= 0) peak = 1;

  const padL = 4, padR = 4, padT = 14, padB = 16;
  const gw = w - padL - padR, gh = h - padT - padB;

  ctx.fillStyle = '#5a5e72';
  ctx.font = '9px monospace';
  ctx.textBaseline = 'top';
  ctx.fillText(label, padL, 1);

  // bars (log-count scaling for dynamic range)
  const bw = gw / n;
  ctx.fillStyle = '#4f8ef7';
  const logPeak = Math.log10(peak+1);
  for (let i=0;i<n;i++) {
    const v = Math.log10(bins[i]+1)/logPeak;
    const bh = v*gh;
    ctx.fillRect(padL + i*bw, padT + gh - bh, Math.max(bw-0.5,0.5), bh);
  }

  // axis labels (min / max)
  ctx.fillStyle = '#5a5e72';
  ctx.textBaseline = 'bottom';
  const fmt = (x) => (Math.abs(x)>=1000||Math.abs(x)<0.01&&x!==0) ? x.toExponential(1) : x.toFixed(2);
  ctx.fillText(fmt(min), padL, h);
  ctx.textAlign = 'right';
  ctx.fillText(fmt(max), w-padR, h);
  ctx.textAlign = 'left';
}
