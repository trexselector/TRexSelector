#!/usr/bin/env python3
"""Write PUBLICATIONS.md, the list of every paper of the T-Rex method family, from the "T-Rex
framework" section of references.bib, so the bibliography stays the single source.

    python3 tools/bib_to_publications.py            # writes PUBLICATIONS.md at the repository root
    python3 tools/bib_to_publications.py --check    # exit 1 when PUBLICATIONS.md is out of date

Entries are grouped into journal articles, conference papers, theses and preprints, newest first,
and printed in the reference style of the READMEs: initials and surnames, the title in quotes, the
venue in italics, volume, article or pages, year, and the DOI or URL as a link.
"""
import re
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIB = ROOT / "references.bib"
OUT = ROOT / "PUBLICATIONS.md"
# the parts of the page: a bibliography section (from its marker to the next) and the page heading
PARTS = [("% T-Rex framework", "% T-Rex for complex-valued data", "The real-valued framework"),
         ("% T-Rex for complex-valued data", "% Solver origins", "The complex-valued framework")]

GROUPS = [("article", "Journal articles"), ("inproceedings", "Conference papers"),
          ("phdthesis", "Theses"), ("online", "Preprints")]


def entries(text, section_start, section_end):
    """The entries of one bibliography section as (type, key, fields)."""
    start, end = text.index(section_start), text.index(section_end)
    out = []
    for m in re.finditer(r"@(\w+)\{([^,]+),\s*(.*?)\n\}", text[start:end], re.S):
        kind, key, body = m.group(1).lower(), m.group(2).strip(), m.group(3)
        fields = {}
        for f in re.finditer(r"(\w+)\s*=\s*\{((?:[^{}]|\{[^{}]*\})*)\}", body, re.S):
            fields[f.group(1).lower()] = re.sub(r"\s+", " ", f.group(2)).strip()
        out.append((kind, key, fields))
    return out


ACCENTS = {'\\"a': "ä", '\\"o': "ö", '\\"u': "ü", '\\"A': "Ä", '\\"O': "Ö", '\\"U': "Ü", "\\'e": "é",
           "\\`e": "è", "\\'a": "á", "\\^o": "ô", "\\ss": "ß"}


def plain(s):
    """TeX braces, accents and the few macros the entries use, removed."""
    s = re.sub(r"\\emph\{([^}]*)\}", r"*\1*", s)
    for tex, char in ACCENTS.items():
        s = s.replace("{" + tex + "}", char).replace(tex, char)
    s = s.replace("\\&", "&").replace("--", "–").replace("~", " ")
    return s.replace("{", "").replace("}", "")


def sentence_case(title):
    """The title in sentence case, as the READMEs print it: text in braces keeps its case, the first
    letter and the first letter after a colon stay capitals, everything else is lowered."""
    out, depth, start = [], 0, True
    for ch in title:
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        elif ch.isalpha():
            if depth == 0:
                out.append(ch if start else ch.lower())
                start = False
                continue
            start = False  # a protected acronym opens the title, the next word is not a start
        elif ch == ":":
            start = True
        out.append(ch)
    return plain("".join(out))


def authors(field):
    names = []
    for a in re.split(r"\s+and\s+", plain(field)):
        if "," in a:
            last, first = [x.strip() for x in a.split(",", 1)]
        else:
            parts = a.split()
            last, first = parts[-1], " ".join(parts[:-1])
        initials = " ".join(p[0] + "." for p in re.split(r"[\s-]+", first) if p)
        names.append(f"{initials} {last}".strip())
    if len(names) <= 2:
        return " and ".join(names)
    return ", ".join(names[:-1]) + ", and " + names[-1]


def line(kind, f):
    out = f"{authors(f['author'])}, \"{sentence_case(f['title'])},\""
    if kind == "article":
        out += f" *{plain(f['journal'])}*"
        if "volume" in f:
            out += f", vol. {f['volume']}"
        if "number" in f:
            out += f", no. {f['number']}"
        if "eid" in f:
            out += f", article {f['eid']}"
        elif "pages" in f:
            out += f", pp. {plain(f['pages'])}"
    elif kind == "inproceedings":
        out += f" in *{plain(f['booktitle'])}*"
        if "address" in f:
            out += f", {plain(f['address'])}"
        if "pages" in f:
            out += f", pp. {plain(f['pages'])}"
    elif kind == "phdthesis":
        out += f" Ph.D. dissertation, {plain(f['school'])}"
        if "address" in f:
            out += f", {plain(f['address'])}"
    elif kind == "online":
        if f.get("eprinttype", "").lower() == "arxiv":
            out += f" arXiv:{f['eprint']}"
    out += f", {f['year']}."
    if kind == "online" and "note" in f:
        note = plain(f["note"]).rstrip(".")
        out += f" {note[0].upper()}{note[1:]}."
    if "doi" in f:
        out += f" (<https://doi.org/{f['doi']}>)"
    elif f.get("eprinttype", "").lower() == "arxiv":
        out += f" (<https://arxiv.org/abs/{f['eprint']}>)"
    elif "url" in f:
        out += f" (<{f['url']}>)"
    return out


def render(text):
    parts = ["# Publications of the T-Rex method family", "",
             "<!-- Generated from the T-Rex sections of references.bib by tools/bib_to_publications.py. A run",
             "     of the script reproduces this page, so edits go into the bibliography or the script. -->", "",
             "Every paper of the Terminating-Random Experiments (T-Rex) framework and its extensions,",
             "newest first.", ""]
    for section_start, section_end, part_heading in PARTS:
        items = entries(text, section_start, section_end)
        parts += [f"## {part_heading}", ""]
        for kind, heading in GROUPS:
            group = sorted((f for k, _, f in items if k == kind), key=lambda f: (-int(f["year"]), f["title"]))
            if not group:
                continue
            parts += [f"### {heading}", ""]
            for f in group:
                # wrapped at the prose width of the repository, the link stays one token
                parts.extend(textwrap.wrap(line(kind, f), width=100, initial_indent="- ",
                                           subsequent_indent="  ", break_long_words=False,
                                           break_on_hyphens=False))
            parts.append("")
    return "\n".join(parts).rstrip("\n") + "\n"


if __name__ == "__main__":
    page = render(BIB.read_text(encoding="utf-8"))
    if "--check" in sys.argv:
        current = OUT.read_text(encoding="utf-8") if OUT.exists() else ""
        if current != page:
            print("PUBLICATIONS.md is out of date, run tools/bib_to_publications.py")
            sys.exit(1)
        print("PUBLICATIONS.md is up to date")
        sys.exit(0)
    OUT.write_text(page, encoding="utf-8")
    n = page.count("\n- ")
    print(f"wrote {OUT.relative_to(ROOT)} with {n} entries")
