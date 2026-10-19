from playwright.sync_api import sync_playwright
import sys
html, out = sys.argv[1], sys.argv[2]
w,h = int(sys.argv[3]), int(sys.argv[4])
with sync_playwright() as p:
    b=p.chromium.launch(); pg=b.new_page(viewport={'width':w,'height':h})
    pg.goto('file://'+html); pg.wait_for_timeout(800)
    pg.pdf(path=out, width=f'{w}px', height=f'{h}px', print_background=True, margin={'top':'0','bottom':'0','left':'0','right':'0'})
    pg.screenshot(path=out.replace('.pdf','.png'), full_page=False)
    b.close()
