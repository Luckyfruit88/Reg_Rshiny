"""Synthetic 1728 x 895 summary contract, real Chromium layout and screenshots."""
import json, os, time
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import urlopen
from urllib.error import URLError
from playwright.sync_api import sync_playwright
root = Path(os.environ.get('SUMMARY_EVIDENCE_DIR', 'test-artifacts/evidence-summary-1260'))
root.mkdir(parents=True, exist_ok=True)
test_url = os.environ.get('SUMMARY_TEST_URL', 'http://127.0.0.1:8765')
assert urlsplit(test_url).hostname in ('127.0.0.1', 'localhost')
# CI starts R and this process together. Wait for readiness, not a fixed delay;
# startup failures still fail within 30 seconds and expose the Shiny log.
deadline = time.monotonic() + 30
while True:
    try:
        with urlopen(test_url, timeout=1) as response:
            if response.status == 200:
                break
    except (URLError, TimeoutError, ConnectionError):
        if time.monotonic() >= deadline:
            raise RuntimeError('Synthetic Shiny fixture did not become ready within 30 seconds')
    time.sleep(.1)
print('SUMMARY_FIXTURE_READY', flush=True)
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
preview_js = r'''el => {
  const box=el.getBoundingClientRect(), walk=document.createTreeWalker(el,NodeFilter.SHOW_TEXT);
  let visible=0; const tops=new Set();
  while(walk.nextNode()) {
    const node=walk.currentNode;
    for(let i=0;i<node.length;i++) {
      const range=document.createRange();range.setStart(node,i);range.setEnd(node,i+1);
      const r=range.getBoundingClientRect();
      // Count wrapping whitespace too: it represents a visible word boundary.
      if(r.height>0 && r.left>=box.left-.5 && r.right<=box.right+.5 && r.top>=box.top-.5 && r.bottom<=box.bottom+.5) {
        visible++; if(node.textContent[i].trim()) tops.add(Math.round(r.top*10)/10);
      }
    }
  }
  return {total:el.textContent.length, visible, lines:tops.size, width:box.width, lineClamp:getComputedStyle(el).webkitLineClamp};
}'''
with sync_playwright() as p:
    browser = p.chromium.launch(headless=True)
    page = browser.new_page(viewport={'width':1728, 'height':895}, device_scale_factor=1)
    external = []
    page.on('request', lambda r: external.append(r.url) if urlsplit(r.url).hostname not in ('127.0.0.1','localhost',None) else None)
    page.goto(test_url, wait_until='networkidle')
    page.locator('.regshiny-evidence-row').nth(6).wait_for()
    expected = json.loads(page.locator('#fixture-expected').text_content())
    report=[]
    for scenario in ['before','after']:
        page.locator('input[name="scenario"][value="'+scenario+'"]').check()
        page.wait_for_function('(s) => document.querySelectorAll(".regshiny-evidence-row .regshiny-step-"+s).length === 12', arg=scenario)
        table = page.locator('.regshiny-evidence')
        box=table.bounding_box();count=table.evaluate(visible_js)
        assert abs(box['width']-1260)<1, box
        entries=[s.strip() for s in table.locator('li[data-role]').all_text_contents()]
        exp=next(x for x in expected if x['state']==scenario)
        assert entries and sorted(entries)==sorted(exp['text'])
        for row in table.locator('.regshiny-evidence-row').all():
            previews=row.locator('.regshiny-evidence-preview li')
            assert previews.count()<=2
            role='evidence' if row.locator('li[data-role="evidence"]').count() else 'context'
            assert all(li.get_attribute('data-role')==role for li in previews.all())
        assert table.locator('[data-question="Pathway"] .regshiny-evidence-preview li[data-role="context"]').count()==2
        assert table.locator('li[data-role="notice"] a').count()==0
        assert all(a.get_attribute('target')=='_blank' and a.get_attribute('rel')=='noopener noreferrer' for a in table.locator('a').all())
        previews=[el.evaluate(preview_js) for el in table.locator('.regshiny-evidence-preview .regshiny-evidence-line').all()]
        for item in previews:
            assert item['visible']>=min(150,item['total']), item
            assert item['lines']<=2 and item['lineClamp']=='2', item
        fonts=table.locator('.regshiny-step,.regshiny-level,.regshiny-evidence-source,.regshiny-evidence-source small,.regshiny-prior-level').evaluate_all('(els)=>els.map(e=>parseFloat(getComputedStyle(e).fontSize))')
        assert min(fonts)>=12, fonts
        page.screenshot(path=str(root/(scenario+'.png')))
        page.screenshot(path=str(root/(scenario+'_full.png')), full_page=True)
        result=dict(state=scenario,container_width=round(box['width'],2),height=round(box['height'],2),visible_characters=count,
            preview_min_visible=min(x['visible'] for x in previews),long_preview_min_visible=min(x['visible'] for x in previews if x['total']>=150),
            preview_measurements=previews,min_badge_font_px=min(fonts),fixture_characters=exp['chars'],evidence_rows=len(entries),screenshot=str(root/(scenario+'.png')))
        report.append(result)
        print(json.dumps(result),flush=True)
        assert box['height']<=1050, result
        for row in table.locator('.regshiny-evidence-row').all():
            row.locator('details > summary').click()
            for line in row.locator('.regshiny-evidence-line,.regshiny-finding').all():
                m=line.evaluate(preview_js)
                assert m['visible']==m['total'], m
            row.locator('details > summary').click()
    page.locator('input[name="view"][value="Map"]').check()
    svg=page.locator('svg.regshiny-evidence-map');svg.wait_for()
    assert abs(svg.bounding_box()['width']-1260)<1
    assert page.locator('.regshiny-map-question').count()==7
    assert page.locator('.regshiny-map-edge.pending').count()==0
    overflow=svg.evaluate('svg=>{const v=svg.viewBox.baseVal;return [...svg.querySelectorAll("text")].filter(t=>{const r=t.getBBox();return r.x<0||r.y<0||r.x+r.width>v.width||r.y+r.height>v.height}).map(t=>t.textContent)}')
    font_sizes=svg.locator('text').evaluate_all('(els)=>els.map(e=>parseFloat(getComputedStyle(e).fontSize)*e.getScreenCTM().a)')
    assert not overflow, overflow
    assert min(font_sizes)>=12, font_sizes
    hover_results=[]
    for kind in ['question','source']:
        node=svg.locator('.regshiny-map-'+kind).first
        attr='data-'+kind+'-node'; node_id=node.get_attribute(attr)
        node.hover();page.wait_for_timeout(180)
        connected=svg.locator('.regshiny-map-edge['+attr+'="'+node_id+'"]')
        other=svg.locator('.regshiny-map-edge:not(['+attr+'="'+node_id+'"])')
        active=connected.evaluate_all('(es)=>es.map(e=>({opacity:+getComputedStyle(e).opacity,width:parseFloat(getComputedStyle(e).strokeWidth),color:getComputedStyle(e).stroke}))')
        faded=other.evaluate_all('(es)=>es.map(e=>+getComputedStyle(e).opacity)')
        assert active and all(e['opacity']==1 and e['width']>=3 and e['color']=='rgb(11, 92, 173)' for e in active),active
        assert faded and all(o<=.15 for o in faded),faded
        hover_results.append(dict(kind=kind,highlighted_edges=len(active),faded_edges=len(faded)))
    assert svg.locator('a .regshiny-map-source').count()>0
    assert svg.locator('a[title*="Release V10"]').count()>0
    page.mouse.move(0,0);page.wait_for_timeout(180)
    page.screenshot(path=str(root/'map.png'),full_page=True)
    page.locator('input[name="scenario"][value="before"]').check()
    page.wait_for_function('document.querySelectorAll(".regshiny-map-edge.pending").length === 12')
    assert not external, external
    (root/'metrics.json').write_text(json.dumps(dict(viewport=[1728,895],container_width=1260,scenarios=report,map_text_overflow=overflow,
        map_min_rendered_font_px=min(font_sizes),map_hover=hover_results,render_external_requests=external),indent=2))
    print('EVIDENCE_SUMMARY_BROWSER_OK: 1260px container; two-line readable previews; CSS Map tracing')
    browser.close()
