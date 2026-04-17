from pathlib import Path

outdir = Path('/Users/mariacuellar/Github/SBS-PediBIRN-analysis/analysis/figures-mediation-analysis')

def svg_header(title, subtitle, total_line, width=2200, height=1400):
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#7a7a7a"/>
    </marker>
    <style>
      .title {{ font-family: Helvetica, Arial, sans-serif; font-size: 40px; font-weight: 700; fill: #111; text-anchor: middle; }}
      .subtitle {{ font-family: Helvetica, Arial, sans-serif; font-size: 24px; fill: #333; text-anchor: middle; }}
      .box {{ fill: #f3f3f3; stroke: #b8b8b8; stroke-width: 3; rx: 10; ry: 10; }}
      .label {{ font-family: Helvetica, Arial, sans-serif; font-size: 28px; fill: #111; text-anchor: middle; dominant-baseline: middle; }}
      .edge {{ stroke: #7a7a7a; stroke-width: 3.5; fill: none; marker-end: url(#arrow); }}
      .edge-bold {{ stroke: #7a7a7a; stroke-width: 4.5; fill: none; marker-end: url(#arrow); }}
      .effect {{ font-family: Helvetica, Arial, sans-serif; font-size: 24px; fill: #222; }}
      .note {{ font-family: Helvetica, Arial, sans-serif; font-size: 22px; fill: #444; text-anchor: middle; }}
    </style>
  </defs>
  <text x="1100" y="70" class="title">{title}</text>
  <text x="1100" y="115" class="subtitle">Boos 2025 structure recreated with PediBIRN subset estimates</text>
  <text x="1100" y="155" class="subtitle">{total_line}</text>
'''

def box(x, y, w, h, lines):
    text = [f'<rect class="box" x="{x}" y="{y}" width="{w}" height="{h}"/>']
    center_x = x + w / 2
    if len(lines) == 1:
        ys = [y + h / 2]
    elif len(lines) == 2:
        ys = [y + h / 2 - 24, y + h / 2 + 24]
    else:
        start = y + h / 2 - 32 * (len(lines)-1) / 2
        ys = [start + i * 32 for i in range(len(lines))]
    for yy, line in zip(ys, lines):
        text.append(f'<text class="label" x="{center_x}" y="{yy}">{line}</text>')
    return '\n  '.join(text)

def edge(path, bold=False):
    cls = 'edge-bold' if bold else 'edge'
    return f'<path class="{cls}" d="{path}"/>'

def txt(x, y, text):
    return f'<text x="{x}" y="{y}" class="effect">{text}</text>'

def footer(note1, note2, height):
    return f'''  <text x="1100" y="{height-110}" class="note">{note1}</text>
  <text x="1100" y="{height-70}" class="note">{note2}</text>
</svg>
'''

# Model B
parts = [svg_header('Figure 2. Model B Mediation Replication', 'Inertial injury with or without contact injury', 'Total effect β = 27.420 | n = 312 complete cases')]
parts.append('  ' + box(90, 520, 420, 180, ['Inertial injury +/-', 'contact injury']))
parts.append('  ' + box(1690, 520, 420, 180, ['Severe retinal', 'injury']))
mediators = [
    ('Extra-axial blood', 760, 170, 420, 130, 'β = 1.201'),
    ('Clinical hypoxia or ischemia', 760, 390, 420, 130, 'β = 2.009'),
    ('Diffuse HIE', 760, 610, 420, 130, 'β = 5.046'),
    ('Seizure', 760, 830, 420, 130, 'β = 1.421'),
]
for label, x, y, w, h, beta in mediators:
    parts.append('  ' + box(x, y, w, h, [label]))
    parts.append('  ' + txt(690, y + 75, beta))
parts.extend([
    '  ' + edge('M 510 555 C 650 300, 690 250, 760 235'),
    '  ' + edge('M 510 595 C 640 470, 690 455, 760 455'),
    '  ' + edge('M 510 625 C 650 680, 690 675, 760 675'),
    '  ' + edge('M 510 665 C 640 910, 690 895, 760 895'),
    '  ' + edge('M 1180 235 C 1420 270, 1540 430, 1690 555'),
    '  ' + edge('M 1180 455 C 1440 470, 1540 500, 1690 585'),
    '  ' + edge('M 1180 675 C 1440 680, 1540 680, 1690 625'),
    '  ' + edge('M 1180 895 C 1430 880, 1540 760, 1690 665'),
    '  ' + edge('M 510 610 C 1030 610, 1200 610, 1690 610', True),
    '  ' + txt(880, 500, 'Direct effect β = 17.743'),
])
parts.append(footer('Effects shown are subset-based replication estimates on the log-odds scale.', 'Model B published structure: exposure = inertial injury +/- contact injury; mediators = extra-axial blood, clinical hypoxia or ischemia, diffuse HIE, seizure.', 1400))
(outdir / 'boos2025-figure2-modelB-replication.svg').write_text('\n'.join(parts))

# Model C
parts = [svg_header('Figure 3. Model C Mediation Replication', 'Diffuse HIE versus no HIE', 'Total effect β = 4.315 | n = 312 complete cases')]
parts.append('  ' + box(90, 520, 360, 180, ['Diffuse HIE', 'vs no HIE']))
parts.append('  ' + box(1750, 520, 360, 180, ['Severe retinal', 'injury']))
mediators = [
    ('Duration of altered consciousness', 700, 110, 520, 130, 'β = 0.663'),
    ('Extra-axial blood', 700, 330, 520, 130, 'β = 0.656'),
    ('Clinical hypoxia or ischemia', 700, 550, 520, 130, 'β = 0.898'),
    ('Seizure', 700, 770, 520, 130, 'β = 0.578'),
]
for label, x, y, w, h, beta in mediators:
    lines = [label] if label != 'Duration of altered consciousness' else ['Duration of altered', 'consciousness']
    parts.append('  ' + box(x, y, w, h, lines))
    parts.append('  ' + txt(600, y + 75, beta))
parts.extend([
    '  ' + edge('M 450 555 C 590 180, 620 175, 700 175'),
    '  ' + edge('M 450 585 C 600 360, 620 395, 700 395'),
    '  ' + edge('M 450 615 C 590 600, 620 615, 700 615'),
    '  ' + edge('M 450 645 C 600 860, 620 835, 700 835'),
    '  ' + edge('M 1220 175 C 1500 210, 1620 430, 1750 555'),
    '  ' + edge('M 1220 395 C 1490 405, 1620 470, 1750 585'),
    '  ' + edge('M 1220 615 C 1480 615, 1620 595, 1750 615'),
    '  ' + edge('M 1220 835 C 1490 810, 1620 740, 1750 645'),
    '  ' + edge('M 450 610 C 1080 610, 1250 610, 1750 610', True),
    '  ' + txt(1030, 500, 'Direct effect β = 1.520'),
])
parts.append(footer('Effects shown are subset-based replication estimates on the log-odds scale.', 'Model C published structure: exposure = diffuse HIE; mediators = duration of altered consciousness, extra-axial blood, clinical hypoxia or ischemia, seizure.', 1400))
(outdir / 'boos2025-figure3-modelC-replication.svg').write_text('\n'.join(parts))

# Model D
parts = [svg_header('Figure 4. Model D Mediation Replication', 'Diffuse HIE versus no HIE or focal HIE', 'Total effect β = 4.315 | n = 312 complete cases')]
parts.append('  ' + box(90, 520, 430, 180, ['Diffuse HIE', 'vs no HIE or focal HIE']))
parts.append('  ' + box(1760, 520, 360, 180, ['Severe retinal', 'injury']))
mediators = [
    ('Duration of altered consciousness', 760, 110, 520, 130, 'β = 0.663'),
    ('Clinical hypoxia or ischemia', 760, 330, 520, 130, 'β = 0.898'),
    ('Extra-axial blood', 760, 550, 520, 130, 'β = 0.656'),
    ('Seizure', 760, 770, 520, 130, 'β = 0.578'),
]
for label, x, y, w, h, beta in mediators:
    lines = [label] if label != 'Duration of altered consciousness' else ['Duration of altered', 'consciousness']
    parts.append('  ' + box(x, y, w, h, lines))
    parts.append('  ' + txt(660, y + 75, beta))
parts.extend([
    '  ' + edge('M 520 555 C 660 180, 690 175, 760 175'),
    '  ' + edge('M 520 585 C 670 360, 690 395, 760 395'),
    '  ' + edge('M 520 615 C 660 600, 690 615, 760 615'),
    '  ' + edge('M 520 645 C 670 860, 690 835, 760 835'),
    '  ' + edge('M 1280 175 C 1510 210, 1630 430, 1760 555'),
    '  ' + edge('M 1280 395 C 1510 405, 1630 470, 1760 585'),
    '  ' + edge('M 1280 615 C 1510 615, 1630 595, 1760 615'),
    '  ' + edge('M 1280 835 C 1510 810, 1630 740, 1760 645'),
    '  ' + edge('M 520 610 C 1120 610, 1300 610, 1760 610', True),
    '  ' + txt(1080, 500, 'Direct effect β = 1.520'),
])
parts.append(footer('Effects shown are subset-based replication estimates on the log-odds scale.', 'Model D published structure: exposure = diffuse HIE vs no HIE or focal HIE; mediators = duration of altered consciousness, clinical hypoxia or ischemia, extra-axial blood, seizure.', 1400))
(outdir / 'boos2025-figure4-modelD-replication.svg').write_text('\n'.join(parts))

# Model E
parts = [svg_header('Figure 5. Model E Mediation Replication', 'Clinical hypoxia or ischemia', 'Total effect β = 5.390 | n = 312 complete cases')]
parts.append('  ' + box(70, 520, 460, 180, ['Clinical hypoxia', 'or ischemia']))
parts.append('  ' + box(1760, 520, 360, 180, ['Severe retinal', 'injury']))
mediators = [
    ('Duration of altered consciousness', 760, 110, 520, 130, 'β = 0.806'),
    ('Diffuse HIE', 760, 330, 520, 130, 'β = 3.586'),
    ('Seizure', 760, 550, 520, 130, 'β = 1.077'),
    ('Extra-axial blood', 760, 770, 520, 130, 'β = -0.459'),
]
for label, x, y, w, h, beta in mediators:
    lines = [label] if label != 'Duration of altered consciousness' else ['Duration of altered', 'consciousness']
    parts.append('  ' + box(x, y, w, h, lines))
    parts.append('  ' + txt(670, y + 75, beta))
parts.extend([
    '  ' + edge('M 530 555 C 660 180, 690 175, 760 175'),
    '  ' + edge('M 530 585 C 670 360, 690 395, 760 395'),
    '  ' + edge('M 530 615 C 660 600, 690 615, 760 615'),
    '  ' + edge('M 530 645 C 670 860, 690 835, 760 835'),
    '  ' + edge('M 1280 175 C 1510 210, 1630 430, 1760 555'),
    '  ' + edge('M 1280 395 C 1510 405, 1630 470, 1760 585'),
    '  ' + edge('M 1280 615 C 1510 615, 1630 595, 1760 615'),
    '  ' + edge('M 1280 835 C 1510 810, 1630 740, 1760 645'),
    '  ' + edge('M 530 610 C 1120 610, 1300 610, 1760 610', True),
    '  ' + txt(1080, 500, 'Direct effect β = 0.381'),
])
parts.append(footer('Effects shown are subset-based replication estimates on the log-odds scale.', 'Model E published structure: exposure = clinical hypoxia or ischemia; mediators = duration of altered consciousness, diffuse HIE, seizure, extra-axial blood.', 1400))
(outdir / 'boos2025-figure5-modelE-replication.svg').write_text('\n'.join(parts))
