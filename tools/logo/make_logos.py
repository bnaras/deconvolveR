"""Assemble four deconvolveR hex-logo candidates from generated curve paths."""
import re, pathlib

d = pathlib.Path(__file__).parent
txt = (d / "curves.txt").read_text()
blocks = re.split(r"== .+? ==\n", txt)[1:]
G_LINE, G_AREA, F_LINE, F_AREA = [b.strip().replace("\n", " ") for b in blocks]

# pointy-top hexagon, viewBox 240x278, center (120,139)
HEX = "M120,3 L237.8,71 L237.8,207 L120,275 L2.2,207 L2.2,71 Z"
HEX_INNER = "M120,10 L231.7,74.5 L231.7,203.5 L120,268 L8.3,203.5 L8.3,74.5 Z"

FONT = "'Avenir Next','Avenir','Futura','Century Gothic','Trebuchet MS',sans-serif"

def svg(body, bg_defs=""):
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 240 278" role="img" aria-label="deconvolveR hex logo">
  <defs>
    {bg_defs}
    <clipPath id="hexclip"><path d="{HEX}"/></clipPath>
  </defs>
{body}
</svg>'''

def wordmark(y, fill, r_fill, size=25, weight=600):
    return (f'<text x="120" y="{y}" text-anchor="middle" font-family="{FONT}" '
            f'font-size="{size}" font-weight="{weight}" fill="{fill}" '
            f'letter-spacing="0.5">deconvolve<tspan fill="{r_fill}" font-weight="800">R</tspan></text>')

logos = {}

# 1. Hidden Signal - deep navy, ghost f, teal g
logos["hidden-signal"] = svg(f'''
  <path d="{HEX}" fill="url(#navy)"/>
  <g clip-path="url(#hexclip)">
    <path d="{F_AREA}" fill="#7FB3E8" fill-opacity="0.14"/>
    <path d="{F_LINE}" fill="none" stroke="#7FB3E8" stroke-opacity="0.55" stroke-width="2" stroke-dasharray="1 5" stroke-linecap="round"/>
    <path d="{G_AREA}" fill="#2EE6C8" fill-opacity="0.13"/>
    <path d="{G_LINE}" fill="none" stroke="#2EE6C8" stroke-width="3.5" stroke-linejoin="round" stroke-linecap="round"/>
    <line x1="30" y1="186" x2="210" y2="186" stroke="#7FB3E8" stroke-opacity="0.35" stroke-width="1"/>
  </g>
  {wordmark(224, "#EAF4FF", "#2EE6C8")}
  <text x="120" y="246" text-anchor="middle" font-family="{FONT}" font-size="11" font-weight="500" fill="#7FB3E8" letter-spacing="2">g-MODELING</text>
  <path d="{HEX}" fill="none" stroke="#2EE6C8" stroke-width="3"/>
''', '<linearGradient id="navy" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#123A63"/><stop offset="1" stop-color="#0A1F38"/></linearGradient>')

# 2. Cardinal - Stanford red, cream curves
logos["cardinal"] = svg(f'''
  <path d="{HEX}" fill="#8C1515"/>
  <g clip-path="url(#hexclip)">
    <path d="{F_LINE}" fill="none" stroke="#F7F0E3" stroke-opacity="0.5" stroke-width="2" stroke-dasharray="1 5" stroke-linecap="round"/>
    <path d="{G_AREA}" fill="#F7F0E3" fill-opacity="0.22"/>
    <path d="{G_LINE}" fill="none" stroke="#F7F0E3" stroke-width="3.5" stroke-linejoin="round" stroke-linecap="round"/>
    <line x1="30" y1="186" x2="210" y2="186" stroke="#F7F0E3" stroke-opacity="0.4" stroke-width="1"/>
  </g>
  {wordmark(226, "#F7F0E3", "#F7F0E3", weight=600)}
  <path d="{HEX_INNER}" fill="none" stroke="#F7F0E3" stroke-opacity="0.55" stroke-width="1.2"/>
  <path d="{HEX}" fill="none" stroke="#5E0E0E" stroke-width="3"/>
''')

# 3. Minimal Light - near-white, indigo fine lines
logos["minimal-light"] = svg(f'''
  <path d="{HEX}" fill="#FAFAF6"/>
  <g clip-path="url(#hexclip)">
    <g transform="translate(0,-16)">
      <path d="{F_AREA}" fill="#8B87D8" fill-opacity="0.16"/>
      <path d="{F_LINE}" fill="none" stroke="#8B87D8" stroke-width="1.8" stroke-dasharray="1 4.5" stroke-linecap="round"/>
      <path d="{G_LINE}" fill="none" stroke="#37339B" stroke-width="3" stroke-linejoin="round" stroke-linecap="round"/>
      <line x1="30" y1="186" x2="210" y2="186" stroke="#37339B" stroke-opacity="0.3" stroke-width="1"/>
      <text x="203" y="200" font-family="Georgia,serif" font-style="italic" font-weight="bold" font-size="19" fill="#37339B" fill-opacity="0.85">&#952;</text>
    </g>
  </g>
  {wordmark(214, "#2A2765", "#37339B", weight=500)}
  <path d="{HEX}" fill="none" stroke="#8C1515" stroke-width="5"/>
''')

# 4. Amber Glow - charcoal, glowing amber g
logos["amber-glow"] = svg(f'''
  <path d="{HEX}" fill="#1C1B22"/>
  <g clip-path="url(#hexclip)">
    <path d="{F_LINE}" fill="none" stroke="#9A97A6" stroke-opacity="0.6" stroke-width="2" stroke-dasharray="1 5" stroke-linecap="round"/>
    <path d="{G_LINE}" fill="none" stroke="#F5A623" stroke-width="6" stroke-opacity="0.28" filter="url(#soft)" stroke-linejoin="round"/>
    <path d="{G_AREA}" fill="#F5A623" fill-opacity="0.10"/>
    <path d="{G_LINE}" fill="none" stroke="#F5A623" stroke-width="3.2" stroke-linejoin="round" stroke-linecap="round"/>
    <line x1="30" y1="186" x2="210" y2="186" stroke="#9A97A6" stroke-opacity="0.4" stroke-width="1"/>
  </g>
  {wordmark(224, "#F2EFE8", "#F5A623")}
  <text x="120" y="246" text-anchor="middle" font-family="Georgia,serif" font-style="italic" font-size="13" fill="#9A97A6">f = g &#8727; noise</text>
  <path d="{HEX}" fill="none" stroke="#F5A623" stroke-width="3"/>
''', '<filter id="soft" x="-20%" y="-20%" width="140%" height="140%"><feGaussianBlur stdDeviation="3"/></filter>')

for name, content in logos.items():
    (d / f"logo-{name}.svg").write_text(content)
    print(f"logo-{name}.svg  {len(content)} bytes")
