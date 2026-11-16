"""Synthetic 1728 x 895 summary contract, real Chromium layout and screenshots."""
import json, os
from pathlib import Path
from urllib.parse import urlsplit
from playwright.sync_api import sync_playwright
root = Path(os.environ.get('SUMMARY_EVIDENCE_DIR', 'test-artifacts/evidence-summary'))
root.mkdir(parents=True, exist_ok=True)
visible_js = r'''root => {
  let count = 0;
  const walk = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
  while (walk.nextNode()) {
    const node = walk.currentNode;
    let el = node.parentElement, clip = root.getBoundingClientRect(), hidden = false;
    for (let a = el; a && a !== root.parentElement; a = a.parentElement) {
      const s = getComputedStyle(a);
      if (s.display === 'none' || s.visibility === 'hidden' || (a.tagName === 'DETAILS' && !a.open && !a.querySelector('summary').contains(el))) {hidden = true; break;}
      if (s.overflowX === 'hidden' || s.overflowY === 'hidden' || s.webkitLineClamp !== 'none') {
        const r = a.getBoundingClientRect();
        clip = {left:Math.max(clip.left,r.left),right:Math.min(clip.right,r.right),top:Math.max(clip.top,r.top),bottom:Math.min(clip.bottom,r.bottom)};
      }
    }
    if (hidden) continue;
    for (let i=0; i<node.length; i++) {
      const range=document.createRange();range.setStart(node,i);range.setEnd(node,i+1);
      const r=range.getBoundingClientRect();
      if (r.width>0 && r.height>0 && r.left>=clip.left-.5 && r.right<=clip.right+.5 && r.top>=clip.top-.5 && r.bottom<=clip.bottom+.5) count++;
    }
  }
  return count;
}'''
with sync_playwright() as p:
    browser = p.chromium.launch(headless=True)
    page = browser.new_page(viewport={'width':1728, 'height':895}, device_scale_factor=1)
    external = []
    page.on('request', lambda r: external.append(r.url) if urlsplit(r.url).hostname not in ('127.0.0.1','localhost',None) else None)
    page.goto(os.environ.get('SUMMARY_TEST_URL','http://127.0.0.1:8765'), wait_until='networkidle')
    page.locator('.regshiny-evidence-row').nth(6).wait_for()
    expected = json.loads(page.locator('#fixture-expected').text_content())
    report=[]
    for scenario in ['before','after']:
        page.locator('input[name="scenario"][value="'+scenario+'"]').check()
        page.wait_for_function('(s) => document.querySelectorAll(".regshiny-evidence-row .regshiny-step-"+s).length === 12', arg=scenario)
        table = page.locator('.regshiny-evidence')
        box=table.bounding_box();count=table.evaluate(visible_js)
        entries=[s.strip() for s in table.locator('li[data-role]').all_text_contents()]
        exp=next(x for x in expected if x['state']==scenario)
        assert entries and sorted(entries)==sorted(exp['text'])
        assert all(n<=2 for n in table.locator('.regshiny-evidence-row').evaluate_all('(rows)=>rows.map(r=>r.querySelectorAll(".regshiny-evidence-preview li").length)'))
        assert table.locator('.regshiny-evidence-preview li:not([data-role="evidence"])').count()==0
        assert all(a.get_attribute('target')=='_blank' and a.get_attribute('rel')=='noopener noreferrer' for a in table.locator('a').all())
        page.screenshot(path=str(root/(scenario+'.png')))
        page.screenshot(path=str(root/(scenario+'_full.png')), full_page=True)
        report.append(dict(state=scenario,height=round(box['height'],2),visible_characters=count,fixture_characters=exp['chars'],finding_range=[min(exp['finding_lengths']),max(exp['finding_lengths'])],evidence_range=[min(exp['evidence_lengths']),max(exp['evidence_lengths'])],evidence_rows=len(entries),screenshot=str(root/(scenario+'.png'))))
        print(json.dumps(report[-1]),flush=True)
        assert box['height']<=900, report[-1]
        assert count<=2600, report[-1]
        # Expanding the one disclosure preserves every line and reveals the full finding.
        compact_row_chars=table.locator('.regshiny-evidence-row').first.evaluate(visible_js)
        table.locator('details').first.locator('summary').click()
        assert table.locator('.regshiny-evidence-row').first.evaluate(visible_js)>compact_row_chars
        assert table.locator('details').first.get_attribute('open') is not None
        table.locator('details').first.locator('summary').click()
    page.locator('input[name="view"][value="Map"]').check()
    page.locator('svg.regshiny-evidence-map').wait_for()
    assert page.locator('.regshiny-map-question').count()==7
    assert page.locator('.regshiny-map-edge.pending').count()==0
    overflow=page.locator('svg.regshiny-evidence-map').evaluate('svg=>{const v=svg.viewBox.baseVal;return [...svg.querySelectorAll("text")].filter(t=>{const r=t.getBBox();return r.x<0||r.y<0||r.x+r.width>v.width||r.y+r.height>v.height}).map(t=>t.textContent)}')
    assert not overflow, overflow
    page.screenshot(path=str(root/'map.png'))
    page.locator('input[name="scenario"][value="before"]').check()
    page.wait_for_function('document.querySelectorAll(".regshiny-map-edge.pending").length === 12')
    assert not external, external
    (root/'metrics.json').write_text(json.dumps(dict(viewport=[1728,895],scenarios=report,map_text_overflow=overflow,render_external_requests=external),indent=2))
    print('EVIDENCE_SUMMARY_BROWSER_OK')
    browser.close()
