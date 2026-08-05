# Logo generation

The deconvolveR hex logo (`man/figures/logo.svg` / `logo.png`) is generated
by the two scripts here. The motif: the spiky prior g(θ) (solid indigo)
recovered from the smooth observed density f (dotted), which is g convolved
with normal noise.

To regenerate:

```sh
python3 gen_curves.py > curves.txt   # density curves as SVG path data
python3 make_logos.py                # writes logo-*.svg candidates here
```

`make_logos.py` produces four styled candidates; the package logo is the
`minimal-light` variant (near-white ground, indigo curves, cardinal border).
Install it with:

```sh
cp logo-minimal-light.svg ../../man/figures/logo.svg
rsvg-convert -w 480 logo-minimal-light.svg -o ../../man/figures/logo.png
Rscript -e 'pkgdown::build_favicons(overwrite = TRUE)'
```

This directory is listed in `.Rbuildignore` and is not shipped in the
package tarball.
