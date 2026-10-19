import sys, io
## Relabelling helper for Figures 3-5: redact the original label text (PyMuPDF)
## and draw the replacement from a reportlab overlay.
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
FD='/usr/share/fonts/truetype/'
pdfmetrics.registerFont(TTFont('LS', FD+'liberation/LiberationSans-Regular.ttf'))
pdfmetrics.registerFont(TTFont('LSB', FD+'liberation/LiberationSans-Bold.ttf'))
pdfmetrics.registerFont(TTFont('DJI', FD+'dejavu/DejaVuSans-Oblique.ttf'))
pdfmetrics.registerFont(TTFont('DJBI', FD+'dejavu/DejaVuSans-BoldOblique.ttf'))

def width(segs, size):
    return sum(pdfmetrics.stringWidth(t, f, size*s) for t,f,s,r in segs)

def draw(c, segs, size, x, y, rot=0, anchor='left'):
    w=width(segs,size)
    c.saveState(); c.translate(x,y); c.rotate(rot)
    dx = -w/2 if anchor=='center' else 0
    for t,f,s,r in segs:
        c.setFont(f, size*s); c.drawString(dx, r*size, t); dx+=pdfmetrics.stringWidth(t,f,size*s)
    c.restoreState()

def D_mhalf(base='D', font='LS', ital='DJI'):
    # D_alpha^{-1/2}
    return [(base,font,1,0),('α',ital,0.75,-0.25),('−1/2',font,0.7,0.45)]

def apply(src, dst, ops):
    """White out and remove (redact) the original text under each 'white' box,
    then draw the replacement text from a reportlab overlay."""
    import pymupdf
    doc=pymupdf.open(src); pg=doc[0]
    W=pg.rect.width; H=pg.rect.height
    for op in ops:
        if op[0]=='white':
            _,x0,y0,x1,y1=op
            pg.add_redact_annot(pymupdf.Rect(x0,y0,x1,y1), fill=(1,1,1))
    pg.apply_redactions(images=0, graphics=0, text=0)
    buf=io.BytesIO(); c=canvas.Canvas(buf, pagesize=(W,H), initialFontName='LS')
    for op in ops:
        if op[0]=='text':
            _,segs,size,x,ytop,rot,anchor=op; c.setFillColorRGB(0,0,0)
            draw(c,segs,size,x,H-ytop,rot,anchor)
    c.save(); buf.seek(0)
    ov=pymupdf.open("pdf", buf.read())
    pg.show_pdf_page(pg.rect, ov, 0, overlay=True)
    doc.save(dst, garbage=3, deflate=True)
