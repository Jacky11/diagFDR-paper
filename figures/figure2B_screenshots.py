from playwright.sync_api import sync_playwright
import sys
with sync_playwright() as p:
    b=p.chromium.launch()
    pg=b.new_page(viewport={'width':1250,'height':900}, device_scale_factor=2)
    pg.goto('file://' + __import__('os').path.abspath('out/report/diagFDR_report.html') + '')
    pg.wait_for_timeout(1500)
    pg.screenshot(path='full.png', full_page=True)
    # element screenshots
    for sel,name in [('#executive-summary','exec.png'),('#headline-summary-table','headline.png'),('#TOC','toc.png'),('h1.title','title.png')]:
        el=pg.query_selector(sel)
        print(sel, bool(el))
        if el: el.screenshot(path=name)
    b.close()
