# Glider Polar Visualiser

A dependency-free static web app for plotting and comparing glider polar curves.
Polar data is fetched at runtime from the
[XCSoar `PolarStore.cpp`](https://github.com/XCSoar/XCSoar/blob/master/src/Polar/PolarStore.cpp)
source file, so there is no build step and no bundled dataset.

## Hosting

The site is served by GitHub Pages and is designed to work unchanged from either
a project-page subpath or a domain root:

- <https://polars.neale.dev/> — custom subdomain (canonical)
- <https://nealedj.github.io/polar-visualiser/> — GitHub Pages project path

Every asset reference in `index.html` is relative (`./style.css`, `./app.js`) and
`app.js` never constructs a same-origin URL, so no base path configuration is
needed for either origin.

### DNS

Add a `CNAME` record for the subdomain pointing at the GitHub Pages host:

```
polars   CNAME   nealedj.github.io.
```

### Repository settings

1. `Settings → Pages → Custom domain`: `polars.neale.dev`.
2. Tick **Enforce HTTPS** once GitHub has provisioned the certificate (this can
   take a few minutes after the DNS record propagates).

The `CNAME` file in the repository root records that custom domain, and must be
present on whichever branch Pages builds from (`main`). `.nojekyll` keeps Pages
from running the files through Jekyll.

Note that once a custom domain is configured, GitHub Pages issues a `301`
redirect from `nealedj.github.io/polar-visualiser/` to `polars.neale.dev` —
the old URL keeps working, it just resolves to the new one. Serving both
origins simultaneously without a redirect is not something GitHub Pages
supports.

## Local development

No tooling required, but the app uses ES modules so it must be served over HTTP
rather than opened as a `file://` URL:

```sh
python3 -m http.server 8000
```

Then open <http://localhost:8000/>.
