// Hover read-outs for the plots. The server sends, per plot output, the marks
// it drew as fractions of the image (origin bottom left) with a text label;
// see src/plot_hover.R. Nothing is computed here: labels are shown as sent.
(function () {
  var store = {}, tip = null, RADIUS = 12;
  function register() {
    Shiny.addCustomMessageHandler('regshiny-hover', function (m) {
      store[m.id] = m.groups || [];
      var box = document.getElementById(m.id);
      if (box) box.classList.toggle('regshiny-hoverable', store[m.id].length > 0);
    });
  }
  if (window.Shiny && Shiny.addCustomMessageHandler) register();
  else document.addEventListener('DOMContentLoaded', register);

  function hide() { if (tip) tip.style.display = 'none'; }
  function find(groups, fx, fy, w, h) {
    var best = null, bestD = RADIUS * RADIUS, i, k, g, d, dx, dy;
    for (i = 0; i < groups.length; i++) {
      g = groups[i]; if (g.kind !== 'xy') continue;
      for (k = 0; k < g.x.length; k++) {
        dx = (g.x[k] - fx) * w; dy = (g.y[k] - fy) * h; d = dx * dx + dy * dy;
        if (d <= bestD) { bestD = d; best = g.label[k]; }
      }
    }
    if (best !== null) return best;
    var area = Infinity, a;
    for (i = 0; i < groups.length; i++) {
      g = groups[i]; if (g.kind !== 'rect') continue;
      for (k = 0; k < g.x0.length; k++) {
        if (fx < g.x0[k] || fx > g.x1[k] || fy < g.y0[k] || fy > g.y1[k]) continue;
        a = (g.x1[k] - g.x0[k]) * (g.y1[k] - g.y0[k]);
        if (a < area) { area = a; best = g.label[k]; }
      }
    }
    if (best !== null) return best;
    for (i = 0; i < groups.length; i++) {
      g = groups[i]; if (g.kind !== 'x') continue;
      if (fx < g.box[0] || fx > g.box[1] || fy < g.box[2] || fy > g.box[3]) continue;
      d = Infinity;
      for (k = 0; k < g.x.length; k++) {
        dx = Math.abs(g.x[k] - fx);
        if (dx < d) { d = dx; best = g.label[k]; }
      }
      if (best !== null) return best;
    }
    return null;
  }
  document.addEventListener('mousemove', function (ev) {
    var img = ev.target;
    if (!img || img.tagName !== 'IMG') return hide();
    var box = img.closest('.shiny-plot-output');
    if (!box || !store[box.id] || !store[box.id].length) return hide();
    var r = img.getBoundingClientRect();
    if (!r.width || !r.height) return hide();
    var label = find(store[box.id], (ev.clientX - r.left) / r.width, 1 - (ev.clientY - r.top) / r.height, r.width, r.height);
    if (label === null) return hide();
    if (!tip) {
      tip = document.createElement('div'); tip.className = 'regshiny-tip'; tip.setAttribute('role', 'tooltip');
      document.body.appendChild(tip);
    }
    tip.textContent = label; tip.style.display = 'block';
    var tw = tip.offsetWidth, th = tip.offsetHeight;
    var x = ev.clientX + 14, y = ev.clientY + 14;
    if (x + tw > window.innerWidth - 4) x = ev.clientX - tw - 14;
    if (y + th > window.innerHeight - 4) y = ev.clientY - th - 14;
    tip.style.left = Math.max(4, x) + 'px'; tip.style.top = Math.max(4, y) + 'px';
  });
  document.addEventListener('scroll', hide, true);
})();
