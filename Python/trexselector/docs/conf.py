"""Sphinx configuration of the trexselector documentation: a user guide of MyST-NB pages executed at
build time, next to the API reference from the docstrings of the installed package."""

import importlib.metadata
import inspect
import re

import nbclient.client
import zmq

project = "trexselector"
author = "Team T-Rex"
copyright = "2026 Team T-Rex"
release = importlib.metadata.version("trexselector")

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
]

# The master bibliography of the repository, cited with {cite} roles, numbered in order of citation.
bibtex_bibfiles = ["../../../references.bib"]
bibtex_default_style = "unsrt"
bibtex_reference_style = "label"

# MyST-NB: every guide page is a notebook in Markdown form, executed at build time.
nb_execution_mode = "force"
nb_execution_timeout = 600
nb_execution_raise_on_error = True
myst_enable_extensions = ["dollarmath", "colon_fence", "deflist"]
myst_heading_anchors = 3

# ipykernel >= 7.4 warns that the kernel talks plain TCP; let the kernel manager provision CurveZMQ keys.
# MyST-NB has no option for it, and "auto" keeps plain TCP for a kernelspec without Curve support.
if zmq.has("curve"):
    _create_kernel_manager = nbclient.client.NotebookClient.create_kernel_manager

    def _create_encrypted_kernel_manager(self):
        km = _create_kernel_manager(self)
        km.transport_encryption = "auto"
        return km

    nbclient.client.NotebookClient.create_kernel_manager = _create_encrypted_kernel_manager

exclude_patterns = ["_build", "README.md", "**/.ipynb_checkpoints"]

html_theme = "furo"
html_title = "trexselector"
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#336C68",
        "color-brand-content": "#336C68",
    },
    "dark_css_variables": {
        "color-brand-primary": "#6BB2AD",
        "color-brand-content": "#6BB2AD",
    },
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_js_files = ["pages_nav.js"]

autodoc_default_options = {"members": True, "undoc-members": True, "show-inheritance": True}
autodoc_member_order = "groupwise"
autoclass_content = "both"
add_module_names = False
python_use_unqualified_type_names = True
maximum_signature_line_length = 70

# pybind11 writes the full runtime types into the signatures it puts in the docstrings; shorten them to
# what a reader needs: ndarray for the annotated numpy arrays, float / int for the typing protocols, bare
# class names instead of the internal module paths, and no "self".
_SUBS = [
    (re.compile(r"typing\.Annotated\[numpy\.typing\.NDArray\[numpy\.(float64|int32|int64|bool_?)\],"
                r"\s*\"[^\"]*\"(,\s*\"[^\"]*\")?\]"), "ndarray"),
    (re.compile(r"numpy\.typing\.NDArray\[numpy\.(float64|int32|int64|bool_?)\]"), "ndarray"),
    (re.compile(r"numpy\.ndarray\[[^\]]*\]"), "ndarray"),
    (re.compile(r"typing\.SupportsFloat \| typing\.SupportsIndex"), "float"),
    (re.compile(r"SupportsFloat \| SupportsIndex"), "float"),
    (re.compile(r"typing\.SupportsIndex"), "int"),
    (re.compile(r"\bSupportsIndex\b"), "int"),
    (re.compile(r"\bSupportsFloat\b"), "float"),
    (re.compile(r"trexselector\._core\.[A-Za-z0-9_.]*\.([A-Za-z0-9_]+)"), r"\1"),
    (re.compile(r"trexselector\.[A-Za-z0-9_.]*\.([A-Za-z0-9_]+)"), r"\1"),
    (re.compile(r"\btyping\."), ""),
    (re.compile(r"\bcollections\.abc\."), ""),
    (re.compile(r"\bnumpy\."), ""),
    (re.compile(r"Annotated\[ndarray(?:,\s*\"[^\"]*\")*\]"), "ndarray"),
    (re.compile(r"Annotated\[([A-Za-z_][A-Za-z0-9_]*)(?:,\s*\"[^\"]*\")*\]"), r"\1"),
    (re.compile(r"= ([A-Z][A-Za-z0-9_]*)\((?:[^()]|\([^()]*\))*\)"), r"= \1()"),
    (re.compile(r"Annotated\[ArrayLike,\s*(\w+),\s*'[^']*'\]"), r"ArrayLike[\1]"),
    (re.compile(r"SupportsInt \| int"), "int"),
    (re.compile(r"\bSupportsInt\b"), "int"),
    (re.compile(r"= <([A-Za-z0-9_]+) object at 0x[0-9a-fA-F]+>"), r"= \1()"),
]


def _shorten(text):
    if not text:
        return text
    for pattern, repl in _SUBS:
        text = pattern.sub(repl, text)
    return text


def _process_signature(app, what, name, obj, options, signature, return_annotation):
    if signature:
        signature = _shorten(signature)
        signature = re.sub(r"^\(self(: [^,)]*)?(, )?", "(", signature)
    return signature, _shorten(return_annotation)


_RUBRIC = re.compile(r"^\.\. rubric:: References\s*$")
_INITIAL = re.compile(r"^([A-Z])\. ")


def _reflow_references(lines):
    # the binding docstrings list their references as "J. Machkour, ..." with the continuation lines
    # indented by two spaces; docutils reads the initial as a list enumerator and the continuation as a
    # block quote. Join each reference into one paragraph, escape the initial, separate the entries.
    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        out.append(line)
        i += 1
        if not _RUBRIC.match(line):
            continue
        out.append("")
        entries = []
        while i < len(lines):
            cur = lines[i]
            if cur.startswith(".. ") or (not cur.strip() and entries and i + 1 < len(lines)
                                          and lines[i + 1].strip() and not lines[i + 1].startswith(" ")
                                          and not _INITIAL.match(lines[i + 1])):
                break
            if not cur.strip():
                i += 1
                continue
            if cur.startswith(" ") and entries:
                entries[-1] += " " + cur.strip()
            else:
                entries.append(cur.strip())
            i += 1
        for entry in entries:
            out.append(_INITIAL.sub(r"\1\\. ", entry))
            out.append("")
    return out


_ABS = re.compile(r"\|([A-Za-z_][A-Za-z0-9_]*)\|")
_OVERLOAD = re.compile(r"^(\d+)\. (\w+\(.*\) -> .*)$")


def _process_docstring(app, what, name, obj, options, lines):
    # napoleon reads a leading "Label: text" of a property as "type: description"; no property of the
    # package declares its type that way, so the property keeps its docstring as written
    if what == "property" and any(line.startswith(":type: ") for line in lines):
        lines[:] = (inspect.getdoc(obj) or "").splitlines()
    # pybind11 repeats the signature as the first docstring line of overloaded or wrapped callables and
    # lists the overloads as "1. f(self: T, ...) -> R", shown here as code without self;
    # |A| in a docstring is an absolute value, not a reStructuredText substitution
    for i, line in enumerate(lines):
        if "typing.Annotated" in line or "trexselector._core" in line or "Supports" in line:
            lines[i] = _shorten(line)
        overload = _OVERLOAD.match(lines[i])
        if overload:
            signature = re.sub(r"\(self(: [^,)]*)?(, )?", "(", overload.group(2))
            lines[i] = f"{overload.group(1)}. ``{signature}``"
        elif "|" in lines[i]:
            lines[i] = _ABS.sub(r"\\|\1\\|", lines[i])
    lines[:] = _reflow_references(lines)


def setup(app):
    app.connect("autodoc-process-signature", _process_signature)
    app.connect("autodoc-process-docstring", _process_docstring)
