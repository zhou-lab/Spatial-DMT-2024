#!/usr/bin/env python3
"""Draw the reviewed Figure 1 v5 as real, editable vector artwork.

All pictograms and matrices are schematic, not experimental observations.
The image-generation concepts and prompts are retained beside the final artwork.
"""
from pathlib import Path
from html import escape
import math
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/manuscript_figures/final_plots'
STEM = OUT / 'Figure1_framework_concept_v5'
W, H = 1800, 1030
C = dict(ink='#17283D', muted='#536477', line='#C9D3DC', blue='#1976C9',
         blue2='#83BCE7', bluebg='#F2F8FD', purple='#7850A5', purple2='#B39BCF',
         purplebg='#F8F5FC', teal='#13858A', tealbg='#F0F9F8', amber='#B57822',
         amberbg='#FCF7EF', gray='#8995A0')
parts = []
text_bounds = []

def add(s): parts.append(s)
def group(name): add(f'<g id="{name}">')
def end(): add('</g>')
def rect(x,y,w,h,fill='white',stroke=None,r=0,sw=1.6,dash=None):
    add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke or "none"}" stroke-width="{sw}"'+(f' stroke-dasharray="{dash}"' if dash else '')+'/>')
def path(d,stroke=C['ink'],sw=2,fill='none',dash=None):
    add(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" stroke-linecap="round" stroke-linejoin="round"'+(f' stroke-dasharray="{dash}"' if dash else '')+'/>')
def line(x1,y1,x2,y2,color=C['ink'],sw=2,dash=None):
    path(f'M{x1},{y1} L{x2},{y2}',color,sw,dash=dash)
def circle(x,y,r,fill,stroke=None,sw=1.5):
    add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke or "none"}" stroke-width="{sw}"/>')
def text(x,y,s,size=25,color=C['ink'],weight='normal',anchor='start'):
    text_bounds.append((x,y,s,size,anchor))
    add(f'<text x="{x}" y="{y}" font-size="{size}" fill="{color}" font-weight="{weight}" text-anchor="{anchor}">{escape(s)}</text>')
def center(x,y,s,**kw): text(x,y,s,anchor='middle',**kw)
def arrow(points,color=C['ink'],sw=2.6):
    d='M'+' L'.join(f'{x},{y}' for x,y in points)
    path(d,color,sw)
    x,y=points[-1]; px,py=points[-2]; a=math.atan2(y-py,x-px)
    k=9
    path(f'M{x-k*math.cos(a-.65)},{y-k*math.sin(a-.65)} L{x},{y} L{x-k*math.cos(a+.65)},{y-k*math.sin(a+.65)}',color,sw)
def reads(x,y,color,w=170):
    for i in range(3):
        yy=y+i*22
        rect(x,yy,49-i*4,8,color,r=3)
        line(x+53,yy+4,x+w-53,yy+4,C['gray'],1.8,'2 4')
        rect(x+w-49,yy,49-i*3,8,color,r=3)
def barcode(x,y,color,h=40):
    for i,ww in enumerate([3,5,2,3,6,2,4,2,5,3,2]):
        rect(x+i*6,y,ww,h,color)
def reference(x,y,w=225):
    line(x,y,x+w,y,C['gray'],6)
    for i in range(6): line(x+15+i*38,y-8,x+15+i*38,y+8,C['gray'],1.5)
def matrix(x,y,w,h,color):
    # A fixed illustrative pattern; white is not a measured zero.
    vals=[[.1,.65,.35,None,.85,.4,.1,.7],[.6,.1,.85,.45,None,.25,.7,.35],
          [.2,.85,.55,.15,.45,.9,.3,None],[None,.35,.1,.85,.55,.2,.7,.4]]
    rgb=tuple(int(color[i:i+2],16) for i in (1,3,5))
    for j,row in enumerate(vals):
        for i,v in enumerate(row):
            f='#EEF1F4' if v is None else '#'+''.join(f'{round(255+(c-255)*v):02x}' for c in rgb)
            rect(x+i*w/8,y+j*h/4,w/8,h/4,f,'white',sw=1)
    rect(x,y,w,h,'none',C['muted'],sw=1.3)
def tissue(x,y,scale=1,outline=False):
    add(f'<g transform="translate({x} {y}) scale({scale})">')
    d='M 7,85 C -2,49 24,13 57,11 C 82,-3 111,18 132,9 C 165,4 193,38 188,64 C 209,95 179,132 155,133 C 128,151 106,129 79,141 C 45,145 8,126 7,85 Z'
    path(d,C['muted'] if outline else '#B57599',1.8,fill='white' if outline else '#F6E5EE')
    if not outline:
        for row in range(5):
            for col in range(6):
                xx=28+col*27+(row%2)*8; yy=29+row*22
                if (col,row) in [(0,0),(5,0),(0,4),(5,4)]: continue
                add(f'<ellipse cx="{xx}" cy="{yy}" rx="11" ry="8.5" fill="#E6BED4" stroke="#BC86A8" stroke-width="1"/>')
                circle(xx+1,yy,3.4,'#966991')
    end()
def spatial_grid(x,y,color,small=False,variant=0):
    step=8 if small else 20
    radius=3 if small else 6
    n=8 if small else 9
    rgb=tuple(int(color[i:i+2],16) for i in (1,3,5))
    for j in range(n):
        for i in range(n):
            if small and ((i in (0,n-1) and j in (0,n-1)) or (i==0 and j==1)): continue
            v=.32+.64*((i*3+j*5+variant*2)%11)/10
            fill='#'+''.join(f'{round(255+(c-255)*v):02x}' for c in rgb)
            circle(x+i*step,y+j*step,radius,fill)
def report(x,y,w,h,title):
    rect(x,y,w,h,'white',C['muted'],r=8)
    rect(x+1,y+1,w-2,42,'#F0F3F6',r=7)
    line(x,y+43,x+w,y+43,C['line'],1.4)
    center(x+w/2,y+29,title,size=26,weight='bold')

add(f'<svg xmlns="http://www.w3.org/2000/svg" width="200mm" height="{200*H/W:.5f}mm" viewBox="0 0 {W} {H}">')
add('<title>Protocol-aware spatial DNA methylation and optional RNA processing framework</title>')
add('<desc>Editable conceptual schematic: protocol-aware barcode parsing produces matched spatial reads and an unmatched diagnostic track. Matched DNA is aligned, deduplicated with spatial barcodes and processed per pixel for methylation. Optional RNA is parsed and aligned independently to produce gene counts. MultiQC and spatial QC summarize processing. All pictograms and matrices are illustrative, not quantitative data.</desc>')
add('<g font-family="Liberation Sans, Arial, Helvetica, sans-serif">')
rect(0,0,W,H)
group('stage-headings')
for x,w,label in [(28,210,'Input'),(330,305,'Trimming'),(710,280,'Alignment'),(1065,305,'PileUp'),(1420,350,'Results')]:
    center(x+w/2,49,label,size=31,weight='bold')
    line(x,77,x+w,77,C['ink'],2)
for a,b in [(250,313),(650,692),(1005,1047),(1384,1407)]: arrow([(a,45),(b,45)],sw=2)
end()

group('input')
center(133,123,'Tissue section',weight='bold')
tissue(32,151,1.06)
arrow([(134,315),(134,345)],C['muted'])
center(133,377,'Spatial pixels',weight='bold')
rect(36,393,195,194,'white',C['line'],r=6)
spatial_grid(53,410,C['blue'])
center(133,624,'DNA reads',weight='bold',color=C['blue'])
reads(45,647,C['blue'],178)
center(133,793,'RNA reads',weight='bold',color=C['purple'])
center(133,820,'(optional)',size=23,color=C['purple'])
reads(45,855,C['purple'],176)
end()

group('protocol-parsing')
rect(330,108,305,205,'white',C['line'],r=8)
center(482.5,145,'Protocol-aware parsing',size=25,weight='bold')
for yy,label,color in [(155,'Spatial-DMT',C['blue']),(222,'SmC-seq',C['amber'])]:
    text(350,yy+20,label,size=26,weight='bold')
    rect(350,yy+40,52,10,color,r=2)
    rect(405,yy+40,52,10,C['blue2'] if color==C['blue'] else '#E4BD77',r=2)
    rect(460,yy+40,150,10,C['gray'],r=2)
arrow([(483,313),(483,336),(403,336),(403,357)],C['blue'])
arrow([(483,336),(563,336),(563,357)],C['amber'])
rect(330,363,146,108,C['bluebg'],C['blue'],r=7)
center(403,393,'Matched',color=C['blue'],weight='bold')
center(403,419,'spatial reads',size=23,color=C['blue'])
barcode(369,433,C['blue'],30)
rect(491,363,144,108,C['amberbg'],C['amber'],r=7)
center(563,393,'Unmatched',size=24,color=C['amber'],weight='bold')
center(563,419,'reads',size=23,color=C['amber'])
rect(515,438,89,8,C['amber'],r=3)
rect(515,454,71,8,C['amber'],r=3)
rect(330,503,305,66,C['tealbg'],C['teal'],r=7)
text(350,549,'λ',size=42,color=C['teal'])
text(395,543,'Lambda control',size=24,color=C['teal'],weight='bold')
center(483,795,'Barcode + UMI',color=C['purple'],weight='bold')
center(483,825,'extraction',color=C['purple'],weight='bold')
barcode(354,859,C['purple'],39)
arrow([(428,879),(451,879)],C['purple'],2)
rect(462,858,62,42,'white',C['purple'],r=4)
center(493,886,'UMI',size=23,color=C['purple'])
arrow([(536,879),(557,879)],C['purple'],2)
for j in range(3): rect(568,860+j*16,45,7,C['purple'],r=2)
end()

group('dna-alignment')
center(850,123,'DNA alignment',weight='bold',color=C['blue'])
center(850,170,'Reference genome',size=23,color=C['muted'])
reference(739,194,225)
for x,y,w in [(743,218,68),(780,239,88),(822,260,105)]: rect(x,y,w,9,C['blue'],r=3)
arrow([(850,288),(850,320)],C['blue'])
center(850,352,'Barcode-aware',size=25,weight='bold')
center(850,380,'deduplication',size=25,weight='bold')
for x,y,w in [(745,403,48),(763,417,70),(745,431,55)]: rect(x,y,w,7,'#AFBBC5',r=2)
arrow([(843,418),(864,418)],C['muted'],2)
rect(875,403,77,8,C['blue'],r=2);rect(888,425,57,8,C['blue'],r=2)
arrow([(850,453),(850,479)],C['blue'])
center(850,511,'Per-pixel alignments',size=25,weight='bold')
for j in range(3):
    yy=544+j*37
    circle(728,yy+3,5,C['blue'])
    line(746,yy+3,980,yy+3,'#D4DCE3',1.4)
    rect(755+j*19,yy-1,65,8,C['blue'],r=2)
    rect(853-j*7,yy-1,57+j*13,8,C['blue2'],r=2)
center(850,662,'Grouped by spatial barcode',size=22,color=C['muted'])
end()

group('rna-alignment')
center(850,794,'RNA alignment',weight='bold',color=C['purple'])
center(850,832,'Reference genome',size=23,color=C['muted'])
reference(738,857,225)
for x,y,w in [(751,880,69),(788,901,82),(827,922,109)]: rect(x,y,w,8,C['purple'],r=3)
end()

group('methylation-and-counting')
center(1217,123,'DNA methylation',weight='bold',color=C['blue'])
line(1098,207,1338,207,C['muted'],2)
for i in range(6):
    x=1118+i*40
    line(x,181,x,207,C['muted'],1.7)
    circle(x,172,8,C['blue'] if i in (0,2,3,5) else 'white',C['blue'],1.8)
arrow([(1217,224),(1217,257)],C['blue'])
center(1217,291,'Methylation matrix',weight='bold')
matrix(1100,316,248,180,C['blue'])
add('<g transform="translate(1081 406) rotate(-90)">');center(0,0,'Pixels',size=23);end()
center(1224,526,'CpG sites',size=23)
circle(1107,553,6,C['blue']); text(1119,560,'Methylated',size=20,color=C['muted'])
circle(1239,553,6,'white',C['blue']); text(1251,560,'Unmethylated',size=20,color=C['muted'])
center(1217,794,'Gene-count matrix',weight='bold',color=C['purple'])
matrix(1100,819,248,104,C['purple'])
add('<g transform="translate(1081 871) rotate(-90)">');center(0,0,'Pixels',size=23);end()
center(1224,952,'Genes',size=23)
end()

group('workflow-connectors')
arrow([(238,669),(279,669),(279,239),(322,239)],C['blue'],2.4)
arrow([(403,471),(403,492),(665,492),(665,244),(704,244)],C['blue'],2.4)
arrow([(990,581),(1030,581),(1030,407),(1063,407)],C['blue'],2.4)
arrow([(244,878),(321,878)],C['purple'],2.4)
arrow([(643,878),(706,878)],C['purple'],2.4)
arrow([(992,878),(1063,878)],C['purple'],2.4)
line(707,729,1370,729,C['line'],1.4,'5 5')
end()

group('multiqc-report')
report(1420,108,350,431,'MultiQC report')
rect(1421,152,118,386,'#F5F7FA')
line(1540,152,1540,539,C['line'],1.2)
for y,label in [(184,'General'),(206,'stats'),(255,'Curated'),(277,'metrics'),(328,'Unmatched'),(379,'BISCUIT')]:
    text(1430,y,label,size=21,color=C['blue'] if y in (184,206,379) else C['muted'])
for yy,label in [(184,'General Statistics'),(365,'Curated metrics')]:
    text(1549,yy,label,size=22,weight='bold')
    for i in range(5):
        rect(1546,yy+23+i*26,23,6,'#AAB3BC',r=1)
        for j,col in enumerate(['#9AC5E9','#5096CD','#91CBCB','#C8DCDE']):
            rect(1580+j*43,yy+19+i*26,36,14,col)
end()

group('spatial-qc-report')
report(1420,566,350,393,'Spatial QC report')
for j,(label,col) in enumerate([('Raw reads',C['blue']),('Mapping',C['teal']),('Duplication',C['teal']),('Methylation',C['blue']),('RNA UMIs',C['purple']),('RNA genes',C['purple'])]):
    xx=1478+(j%2)*175; yy=650+(j//2)*112
    center(xx+28,yy-12,label,size=24,weight='bold',color=col)
    spatial_grid(xx,yy,col,True,j)
    if j>=4: center(xx+28,yy+77,'(optional)',size=21,color=C['purple'])
end()

group('caption-note')
text(1770,1001,'Schematic; no quantitative data.',size=22,color=C['muted'],anchor='end')
end()
end();add('</svg>')
OUT.mkdir(parents=True,exist_ok=True)
svg='\n'.join(parts)+'\n'
ET.fromstring(svg)
assert '<image' not in svg and 'data:image' not in svg
STEM.with_suffix('.svg').write_text(svg)
for fmt in ['pdf','eps']:
    subprocess.run(['rsvg-convert','-d','72','-p','72','-f',fmt,'-o',str(STEM.with_suffix('.'+fmt)),str(STEM.with_suffix('.svg'))],check=True)
subprocess.run(['rsvg-convert','-w','4724','-h','2703','-o',str(STEM.with_suffix('.png')),str(STEM.with_suffix('.svg'))],check=True)
# Set print-resolution metadata without resampling the vector export.
from PIL import Image
im=Image.open(STEM.with_suffix('.png'))
im.save(STEM.with_suffix('.png'),dpi=(600,600))
subprocess.run(['rsvg-convert','-w','1800','-o',str(OUT/'figure1_v5_source/vector_preview.png'),str(STEM.with_suffix('.svg'))],check=True)
print(f'Created {STEM}.svg, .pdf, .eps, .png; SVG has no embedded raster images.')
