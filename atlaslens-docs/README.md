# AtlasLens documentation

This folder is a ready-to-fill documentation site for **AtlasLens**, built with
[MkDocs](https://www.mkdocs.org/) and the [Material theme](https://squidfunk.github.io/mkdocs-material/).
It produces a website that looks and behaves just like the SNEEP Read-the-Docs page
(sidebar, search, Next/Previous links, mobile-friendly).

You only ever edit plain **Markdown** files — no coding required.

---

## What's in here

```
atlaslens-docs/
├── mkdocs.yml            ← site config: title, colors, sidebar order
├── requirements.txt      ← what the site needs to build
├── .readthedocs.yaml     ← lets Read the Docs build it automatically
└── docs/
    ├── index.md          ← home page
    ├── getting-started.md
    ├── features.md
    ├── step-by-step.md   ← the main screenshot walkthrough
    ├── faq.md
    └── assets/           ← PUT YOUR SCREENSHOTS HERE
```

## How to fill it in

1. Open any `.md` file in `docs/` with a text editor.
2. Replace the `_[bracketed placeholder]_` text with your real words.
3. Drop your screenshots into `docs/assets/` using the filenames already referenced
   (see `docs/assets/README.txt`). They'll show up automatically.
4. To add a screenshot anywhere, write:  `![description](assets/your-image.png)`

## Preview it on your computer

You need [Python](https://www.python.org/downloads/) installed. Then, in a terminal,
from inside this folder:

```bash
pip install -r requirements.txt   # one-time setup
mkdocs serve                      # starts a live preview
```

Open the link it prints (usually http://127.0.0.1:8000). The preview updates
as you save your files.

## Publish it (free)

**Option A — Read the Docs (same host as SNEEP):**

1. Put this folder in a GitHub repository.
2. Sign in at [readthedocs.org](https://readthedocs.org) and "Import" the repo.
3. It builds automatically using `.readthedocs.yaml`. Every push updates the site.

**Option B — GitHub Pages (one command):**

```bash
mkdocs gh-deploy
```

This builds the site and pushes it to a `gh-pages` branch; GitHub serves it for free.

---

Need help writing the actual content or wiring up real screenshots? Just ask.
