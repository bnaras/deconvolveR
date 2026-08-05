"""Generate SVG path data for deconvolveR logo curves.

g(theta): spiky bimodal prior (mixture of two narrow normals).
f(x):     g convolved with a wide normal kernel (the observed, blurred density).
Both scaled into hex-sticker coordinates (viewBox 240x278, pointy-top hexagon).
"""
import math

N = 121


def norm(x, mu, sd):
    return math.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * math.sqrt(2 * math.pi))


def densities():
    xs = [i / (N - 1) for i in range(N)]
    g = [0.62 * norm(x, 0.33, 0.055) + 0.38 * norm(x, 0.72, 0.045) for x in xs]
    # convolution of the mixture with N(0, 0.13) kernel is analytic:
    f = [0.62 * norm(x, 0.33, math.hypot(0.055, 0.13))
         + 0.38 * norm(x, 0.72, math.hypot(0.045, 0.13)) for x in xs]
    return xs, g, f


def to_path(xs, ys, x0, x1, ybase, ymax_px, ymax_val, close=False):
    pts = []
    for x, y in zip(xs, ys):
        px = x0 + x * (x1 - x0)
        py = ybase - (y / ymax_val) * ymax_px
        pts.append((round(px, 1), round(py, 1)))
    d = "M{},{}".format(*pts[0]) + " L" + " ".join(f"{px},{py}" for px, py in pts[1:])
    if close:
        d += f" L{pts[-1][0]},{ybase} L{pts[0][0]},{ybase} Z"
    return d


xs, g, f = densities()
gmax = max(g)

X0, X1 = 38, 202          # horizontal extent inside hex
BASE = 186                # curve baseline y
print("== g curve (line) ==")
print(to_path(xs, g, X0, X1, BASE, 96, gmax))
print("== g curve (closed area) ==")
print(to_path(xs, g, X0, X1, BASE, 96, gmax, close=True))
print("== f curve (line) ==")
print(to_path(xs, f, X0, X1, BASE, 96, gmax))
print("== f curve (closed area) ==")
print(to_path(xs, f, X0, X1, BASE, 96, gmax, close=True))
