#!/usr/bin/env python3
"""Append the DOI (or, without one, the URL) of every entry of Doxygen's citation list, read from the
bibliography file. Doxygen renders the list with its own BibTeX style, which prints neither field.

    python3 add_doi_links.py <citelist.html> <references.bib>
"""
import re
import sys


def fields(bib_text):
    """key -> (doi, url) for every entry of the bibliography."""
    out = {}
    for m in re.finditer(r"@\w+\{([^,]+),(.*?)\n\}", bib_text, re.S):
        key, body = m.group(1).strip(), m.group(2)
        doi = re.search(r"\bdoi\s*=\s*\{([^}]+)\}", body)
        url = re.search(r"\burl\s*=\s*\{([^}]+)\}", body)
        out[key] = (doi.group(1).strip() if doi else None, url.group(1).strip() if url else None)
    return out


def main(html_path, bib_path):
    html = open(html_path, encoding="utf-8").read()
    refs = fields(open(bib_path, encoding="utf-8").read())
    added = 0

    def link(key):
        doi, url = refs.get(key, (None, None))
        if doi:
            return f' <a href="https://doi.org/{doi}">doi: {doi}</a>.'
        if url:
            return f' <a href="{url}">{url}</a>.'
        return ""

    def repl(m):
        nonlocal added
        extra = link(m.group(1))
        if not extra or "doi.org" in m.group(0) or 'href="http' in m.group(0):
            return m.group(0)
        added += 1
        return m.group(0)[: -len("</dd>")].rstrip() + extra + "\n</dd>"

    # every entry: the anchor with its key, then the text up to the end of the definition
    new = re.sub(r'id="CITEREF_([^"]+)".*?</dd>', repl, html, flags=re.S)
    open(html_path, "w", encoding="utf-8").write(new)
    print(f"citation list: {added} DOI or URL links added")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
