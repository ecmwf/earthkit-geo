#!/usr/bin/env python3

# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Print RST tables of the predefined grids, areas and projections.

The entries are read from the catalogue files installed with eckit and grouped
by what they define, never by the file they come from: a catalogue without a
title of its own joins the merged "Additional ..." table of its kind. Used
through the ``module-output`` directive::

    .. module-output:: generate_catalogue_rst grids
    .. module-output:: generate_catalogue_rst areas
    .. module-output:: generate_catalogue_rst projections

Without arguments all three kinds are printed, in that order. The files are
read, not validated: the tables only show what the catalogues contain.
"""

import sys
from pathlib import Path

import yaml

# Directory holding the catalogue files. When it is None the directory shipped
# with eckit is used. Assigning to this variable (e.g. with monkeypatch) points
# the generator at another set of files, which is how the tests do it.
CATALOGUE_DIR = None

# The top level keys of a catalogue file, and the kind of entry each lists.
# A file contributes to the kind of the keys it has.
ENTRY_KEYS = {
    "grid_names": "grids",
    "area_names": "areas",
    "projection_names": "projections",
}

# The kinds, in the order they are printed when none is asked for.
KINDS = ("grids", "areas", "projections")

# Catalogue files which are not listed: FESOM grids are not yet fully
# supported. The names are file stems, normalised with _normal().
EXCLUDED_FILES = ("fesom",)

# The tag of the empty value an entry name gets when it carries an anchor.
NULL_TAG = "tag:yaml.org,2002:null"

# The suffix of the top level keys which list entries by uid only. Their items
# carry no data, so they are not listed.
UID_SUFFIX = "_uids"

# Entry fields left out of the tables: internal or too long to be readable.
HIDDEN = ("uid", "url", "url_connectivity", "file")

# Characters RST would try to interpret inside a plain text table cell.
MARKUP = ("`", "*", "_", "|", "\\")

# Section title per catalogue, describing what its entries define. Titles never
# show a file name: a catalogue which is not listed here contributes its rows
# to the merged table at the end of its kind (see _catalogues and _title).
TITLES = {
    "grid": "Regional projected grids",
    "hydro": "Hydrological grids",
    "orca": "ORCA ocean grids",
    "icon-ch": "ICON grids of MeteoSwiss",
    "icon-dwd": "ICON grids of DWD",
    "icon-mpim": "ICON grids of MPI-M",
    "area": "Named areas",
    "projection": "Named projections",
}

# One sentence per catalogue, shown above its table.
INTROS = {
    "grid": "Regional grids, given as a map projection and regular ``x`` and ``y`` coordinate ranges.",
    "hydro": "Hydrological grids, given as a grid ``type``, a number of cells ``n`` and the area they cover.",
    "orca": "The ORCA ocean configurations, with one entry per grid arrangement.",
    "area": "Named areas which can be used in place of an explicit area specification.",
    "projection": "Named projections which can be used in place of an explicit projection specification.",
}

ICON_INTRO = "The ICON grids, given as a shape and, where available, an ICON type or a bounding box."

# Above the merged table of the catalogues without a section of their own.
ADDITIONAL_INTRO = "Further predefined {kind}."


def catalogue_dir(path=None):
    """Return the directory of the catalogue files, or None when it is unknown."""
    if path is not None:
        return Path(path)

    if CATALOGUE_DIR is not None:
        return Path(CATALOGUE_DIR)

    try:
        import eckitlib
    except Exception:
        return None

    return Path(eckitlib.__file__).parent / "share" / "eckit" / "geo"


def _key(value):
    """Convert a YAML key into a printable name."""
    # Booleans need care: the "false" key of grid.yaml is loaded as False.
    if value is None:
        return ""

    if isinstance(value, bool):
        return "true" if value else "false"

    return str(value)


def _node_name(node):
    """The printable name of a scalar node, or None when it is not a name."""
    if not isinstance(node, yaml.ScalarNode):
        return None

    return _key(node.value)


def _text(value):
    """Compact, printable text for a YAML value."""
    if value is None:
        return ""

    if isinstance(value, bool):
        return "true" if value else "false"

    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_text(v) for v in value) + "]"

    if isinstance(value, dict):
        # A short description of a nested mapping: what it is, not every detail.
        for name in ("type", "name"):
            if name in value:
                return _text(value[name])
        return ", ".join(f"{_key(k)}={_text(v)}" for k, v in value.items())

    return str(value)


def _cell(text, literal=True):
    """Format one table cell, using an inline literal unless that is unsafe."""
    text = str(text).strip()

    if not text:
        return "-"

    if text.startswith("`") or text.endswith("`") or "``" in text:
        literal = False
    elif not literal and any(m in text for m in MARKUP):
        # Plain text with markup characters could break the table: quote it.
        literal = True

    if literal:
        return f"``{text}``"

    return text


def _entry(item):
    """Return the (name, fields) of one catalogue entry.

    The catalogue files use two layouts. In grid.yaml, hydro.yaml, area.yaml and
    projection.yaml an item is ``{NAME: {...fields...}}``. In the ORCA and ICON
    files the YAML anchor sits on the entry name, so an item is one flat
    dict where the name is the key with a null value::

        - ORCA2_T: &ORCA2_T
          type: ORCA
          name: ORCA2

    An item may also hold only a name, as the aliases of ICON-CH and the
    identifiers of the grid_uids lists do. :func:`_rows` gives an alias the fields
    of the entry it refers to.
    """
    if not isinstance(item, dict):
        return _key(item), {}

    if len(item) == 1:
        name, fields = next(iter(item.items()))
        if isinstance(fields, dict):
            return _key(name), {_key(k): v for k, v in fields.items()}
        return _key(name), {}

    names = [k for k, v in item.items() if v is None]
    if names:
        first = names[0]
        fields = {k: v for k, v in item.items() if k != first}
        return _key(first), {_key(k): v for k, v in fields.items()}

    fields = {_key(k): v for k, v in item.items()}
    return _text(fields.pop("name", "")), fields


def _kind_entries(data, kind):
    """Return the (name, fields) entries one file lists for one kind.

    None when the file has no list of that kind, [] when its list is empty.
    """
    if not isinstance(data, dict):
        return None

    for key, k in ENTRY_KEYS.items():
        if k != kind:
            continue
        items = data.get(key)
        if isinstance(items, list):
            return [_entry(item) for item in items]

    return None


def _anchor_of(item):
    """Return the (name, node) of one entry of the composed node graph.

    The node is the value of the entry name key, which is where the anchor of an
    entry sits: in ``- ORCA2_T: &ORCA2_T`` the anchor belongs to the empty value
    of the ORCA2_T key, not to a mapping of fields.
    """
    if not isinstance(item, yaml.MappingNode) or not item.value:
        return None, None

    pairs = item.value

    if len(pairs) == 1:
        key, node = pairs[0]
        return _node_name(key), node

    for key, node in pairs:
        if node.tag == NULL_TAG:
            return _node_name(key), node

    return None, None


def _aliases(text):
    """Return {alias name: target name} for the entries which are other names.

    A catalogue file can list the same entry under several names::

        - ICON-CH1-V1: &ICON_CH1_V1
          shape: [1147980]

        - ICON-CH1: *ICON_CH1_V1

    Because the anchor is on the null value of the target name key, safe_load
    gives the alias item nothing but its name. In the composed node graph the
    alias value IS the anchored node, so the two items share that node object,
    and the first one to use it is the target (YAML needs the anchor before the
    alias).
    """
    try:
        root = yaml.compose(text, Loader=yaml.SafeLoader)
    except Exception:
        return {}

    if not isinstance(root, yaml.MappingNode):
        return {}

    aliases = {}

    for key, items in root.value:
        if not isinstance(items, yaml.SequenceNode) or (_node_name(key) or "").endswith(UID_SUFFIX):
            continue

        seen = {}
        for item in items.value:
            name, node = _anchor_of(item)
            if name is None or node is None:
                continue

            target = seen.get(id(node))
            if target is not None and target != name:
                aliases[name] = target
            else:
                seen.setdefault(id(node), name)

    return aliases


def _alias_target(name, aliases, fields_of):
    """Return the (name, fields) of the entry the alias called name refers to.

    Unresolvable aliases give (None, None), so that they are shown as they are.
    """
    target = aliases.get(name)
    seen = {name}

    while target is not None and target not in seen:
        seen.add(target)

        if target not in fields_of:
            # The anchor refers to something which is not a listed entry.
            return None, None

        fields = fields_of[target]
        if fields:
            return target, fields

        # The target is an alias itself: follow the chain.
        target = aliases.get(target)

    return None, None


def _rows(entries, aliases):
    """Return the (name, fields, alias) table rows of one catalogue file."""
    fields_of = {}
    for name, fields in entries:
        fields_of.setdefault(name, fields)

    rows = []
    for name, fields in entries:
        target, target_fields = _alias_target(name, aliases, fields_of) if name in aliases else (None, None)

        if target is None:
            rows.append((name, fields, None))
        else:
            # An alias keeps its own fields when it has some.
            rows.append((name, fields or target_fields, target))

    return rows


def _name_of(name, fields):
    return name


def _field(key):
    """Column showing one entry field."""
    return lambda name, fields: _text(fields.get(key))


def _area(name, fields):
    if "area" in fields:
        return _text(fields["area"])

    if "west" in fields:
        return f"west = {_text(fields['west'])}"

    return ""


def _details(name, fields):
    return "; ".join(f"{k} = {_text(v)}" for k, v in fields.items() if k != "type" and k not in HIDDEN)


def _generic_columns():
    """The columns of a table which is not of one specific catalogue kind.

    The name, the type, and all the other fields as one details column.
    """
    return [
        ("Name", 3, _name_of, "name"),
        ("Type", 2, _field("type"), "literal"),
        ("Details", 5, _details, "literal"),
    ]


def _columns(stem, rows):
    """The (header, weight, getter, style) tuples of the table of a catalogue.

    The stem is None for the merged table of the catalogues without a section
    of their own, which uses the generic columns. The style tells how a cell
    is printed: "name" is the entry name, "literal" puts a value in an inline
    literal and "text" keeps prose as it is.
    """
    if stem is None:
        return _generic_columns()

    keys = set()
    for _, fields, _ in rows:
        keys.update(fields)

    name = ("Name", 3, _name_of, "name")

    if stem == "grid":
        return [
            name,
            ("Projection", 3, _field("projection"), "literal"),
            ("x", 2, _field("x"), "literal"),
            ("y", 2, _field("y"), "literal"),
        ]

    if stem == "hydro":
        return [
            name,
            ("Type", 2, _field("type"), "literal"),
            ("n", 1, _field("n"), "literal"),
            ("Area", 4, _area, "literal"),
        ]

    if stem == "orca":
        return [
            name,
            ("Model", 2, _field("name"), "literal"),
            ("Arrangement", 1, _field("arrangement"), "literal"),
            ("Shape", 2, _field("shape"), "literal"),
        ]

    if stem.startswith("icon"):
        columns = [name]

        if "description" in keys:
            columns.append(("Description", 6, _field("description"), "text"))

        columns.append(("Shape", 2, _field("shape"), "literal"))

        if "icon_type" in keys:
            columns.append(("Type", 3, _field("icon_type"), "literal"))
        elif "bounding_box" in keys:
            columns.append(("Bounding box", 3, _field("bounding_box"), "literal"))
        else:
            columns.append(("Type", 3, _field("type"), "literal"))

        return columns

    # area.yaml, projection.yaml and any other listed catalogue
    return _generic_columns()


def _widths(weights):
    """Turn column weights into the integers of the list-table widths option."""
    total = float(sum(weights))
    values = [max(1, int(round(100 * w / total))) for w in weights]

    drift = 100 - sum(values)
    if drift:
        largest = values.index(max(values))
        values[largest] += drift

    return values


def _row(cells):
    return "   * - " + "\n     - ".join(cells)


def _table(columns, rows):
    headers = [c[0] for c in columns]
    widths = " ".join(str(w) for w in _widths([c[1] for c in columns]))

    lines = [
        ".. list-table::",
        "   :header-rows: 1",
        f"   :widths: {widths}",
        "",
        _row(headers),
        "",
    ]

    for row in rows:
        lines.append(_row([_text_of(column, row) for column in columns]))
        lines.append("")

    return lines


def _text_of(column, row):
    """The printed cell of one column for one row."""
    name, fields, alias = row
    _, _, getter, style = column

    text = _cell(getter(name, fields), literal=style != "text")

    if alias and style == "name":
        text += f" (alias of {_cell(alias)})"

    return text


def _title(stem, kind):
    """The section title of one catalogue, which is never its file name.

    The stem is None for the merged table of the catalogues without a title of
    their own; every other catalogue is listed in TITLES.
    """
    if stem is None:
        return f"Additional {kind}"

    return TITLES[stem]


def _intro(stem, kind):
    """The sentence shown above the table of one catalogue."""
    if stem is None:
        return ADDITIONAL_INTRO.format(kind=kind)

    if stem in INTROS:
        return INTROS[stem]

    return ICON_INTRO


def _summary(rows, kind):
    """The line under a table, counting its entries."""
    if rows:
        number = len(rows)
        word = kind if number != 1 else kind[:-1]  # grids -> grid
        text = f"The table lists {number} {word}."
    else:
        text = f"There are no {kind} to list."

    return [text, ""]


def _note(text):
    lines = [".. note::", ""]
    lines += [f"    {line}" for line in str(text).split("\n")]
    lines += [""]

    return lines


def _normal(name):
    """Normalised name: lower case, without a .yaml suffix."""
    name = str(name).strip().lower()
    return name[: -len(".yaml")] if name.endswith(".yaml") else name


def _files(directory):
    """The catalogue files to list, sorted, without the excluded ones."""
    files = sorted(directory.glob("*.yaml"), key=lambda p: (p.name.lower(), p.name))

    return [f for f in files if _normal(f.stem) not in EXCLUDED_FILES]


def _read(path):
    """Return the (text, data, error) of one catalogue file."""
    try:
        text = path.read_text(encoding="utf-8")
    except Exception as e:
        return None, None, str(e)

    try:
        return text, yaml.safe_load(text), None
    except Exception as e:
        return text, None, str(e)


def _catalogues(files):
    """Group the table rows of the catalogue files by kind.

    Return ({kind: [(stem, rows)]}, [(file name, reason)]), with the sections
    in file order. A catalogue which has no title of its own (its stem is not
    in TITLES) does not get a section: its rows join the merged table at the
    end of the kind, whose stem is None and whose title shows no file name.
    The second list holds the files which could not be read and the ones
    listing none of the kinds; their names are the only file names which may
    be printed, in a note of the combined listing.
    """
    groups = {kind: [] for kind in KINDS}
    merged = {kind: None for kind in KINDS}
    unlisted = []

    for path in files:
        text, data, error = _read(path)

        if error is not None:
            reason = error.strip().split("\n")[0] or "no reason given"
            unlisted.append((path.name, f"could not be read: {_cell(reason)}"))
            continue

        listed = False
        for kind in KINDS:
            entries = _kind_entries(data, kind)
            if entries is None:
                continue
            listed = True
            rows = _rows(entries, _aliases(text))
            stem = _normal(path.stem)
            if stem in TITLES:
                groups[kind].append((stem, rows))
            else:
                merged[kind] = (merged[kind] or []) + rows

        if not listed:
            unlisted.append((path.name, "lists no grids, areas or projections"))

    for kind in KINDS:
        if merged[kind] is not None:
            groups[kind].append((None, merged[kind]))

    return groups, unlisted


def _render_section(stem, kind, rows, titled):
    """The lines of one table, with a title when the kind has several.

    The stem is None for the merged table of the unknown catalogues.
    """
    lines = []

    if titled:
        title = _title(stem, kind)
        lines += [title, "-" * len(title), ""]

    lines.append(_intro(stem, kind))
    lines.append("")

    if rows:
        lines += _table(_columns(stem, rows), rows)

    return lines + _summary(rows, kind)


def _kinds_requested(args):
    """Return the (wanted, unknown) kinds named by the arguments."""
    wanted = []
    unknown = []

    for arg in args:
        kind = _normal(arg)

        if not kind:
            continue

        if kind in KINDS:
            if kind not in wanted:
                wanted.append(kind)
        elif kind not in unknown:
            unknown.append(kind)

    return wanted, unknown


def _join(lines):
    """Assemble the printed lines into the RST of the page section."""
    # The leading blank line separates a section title from the content
    # which precedes the directive in the page.
    return "\n".join([""] + list(lines))


def render(*args, path=None):
    """Return the RST listing the predefined grids, areas and projections.

    The optional arguments select the kinds to print (``grids``, ``areas`` and
    ``projections``) and the order; without them all three are printed. An
    unknown kind gives a note, not an error. ``path`` overrides the directory
    of the catalogue files; without it, and unless :data:`CATALOGUE_DIR` is
    set, the files installed with eckit are used. A missing directory gives a
    note too, never an error.
    """
    directory = catalogue_dir(path)

    if directory is None:
        return _join(_note("The catalogue files are installed with eckit, which is not available here."))

    if not directory.is_dir():
        return _join(_note(f"No catalogue directory was found at ``{directory}``."))

    wanted, unknown = _kinds_requested(args)

    if not wanted and not unknown:
        wanted = list(KINDS)

    lines = []

    if unknown:
        names = ", ".join(f"``{u}``" for u in unknown)
        kinds = ", ".join(f"``{k}``" for k in KINDS)
        lines += _note(f"No kind of entry matches {names}. The kinds are {kinds}.")

    files = _files(directory)

    if not files:
        return _join(lines + _note(f"The directory ``{directory}`` contains no catalogue files."))

    groups, unlisted = _catalogues(files)

    for kind in wanted:
        if len(wanted) > 1:
            # A plain RST comment separating the kinds in a combined listing.
            lines += [f".. {kind}", ""]

        sections = groups[kind]

        for stem, rows in sections:
            lines += _render_section(stem, kind, rows, titled=len(sections) > 1)

        if not sections:
            lines += _note(f"No catalogue file lists any {kind}.")

    if set(wanted) == set(KINDS) and unlisted:
        lines += _note("\n".join(f"``{name}`` {reason}" for name, reason in unlisted))

    return _join(lines)


def execute(*args, path=None):
    print(render(*args, path=path))


if __name__ == "__main__":
    execute(*sys.argv[1:])
