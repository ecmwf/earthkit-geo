#!/usr/bin/env python3

# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Tests for the catalogue contents generator used by the documentation."""

import io

import pytest

from earthkit.geo.sphinxext import generate_catalogue_rst as mod

# The sections of the grids page built from the files installed with eckit,
# in order. FESOM is excluded (not yet fully supported), so it is not here.
GRID_SECTIONS = (
    "Regional projected grids",
    "Hydrological grids",
    "ICON grids of MeteoSwiss",
    "ICON grids of DWD",
    "ICON grids of MPI-M",
    "ORCA ocean grids",
)

# Nested layout: one item is {NAME: {...fields...}}. grid.yaml also has a
# "false" key, which YAML loads as the boolean False.
GRID_YAML = """
grid_names:
    - LAEA-EFAS-5km:
          x: [0, 4995000, 5000]
          y: [0, -4745000, -5000]
          false: [1818500.021, -2287499.955]
          projection: {type: laea, lat_0: 52, lon_0: 10, figure: grs80, reverse: true}

    - swisslv95:
          x: [2439000, 2867000, 1000]
          y: [1040500, 1333500, 1000]
          projection: swisslv95
"""

HYDRO_YAML = """
grid_names:
    - EFAS-5:
          type: arakawa_c
          n: 5400
          area: [72.25, -25.25, 22.75, 50.25]

    - CAMA-15min:
          type: arakawa_c
          n: 360
          west: -180
"""

# Flat layout: the anchor on the entry name makes the item one dict where the
# name is the key with a null value. An entry can also be an alias, holding
# nothing but its name. The uids are names without any data.
ORCA_YAML = """
grid_names:
  - ORCA2_T: &ORCA2_T
    type: ORCA
    name: ORCA2
    intgrid: O42
    arrangement: T
    shape: [149, 182]
    uid: "d5bde4f52ff3a9bea5629cd9ac514410"
    url: https://example.org/ORCA2_T.atlas

  - ORCA1_F: &ORCA1_F
    type: ORCA
    name: ORCA1
    arrangement: F
    shape: [292, 362]
    uid: "174487fbace54b00d959d971e88b71e7"

  - ORCA2_LIM: *ORCA2_T

grid_uids:
  - d5bde4f52ff3a9bea5629cd9ac514410:
  - 174487fbace54b00d959d971e88b71e7:
"""

FESOM_YAML = """
grid_names:
  - pi_C: &pi_C
    type: FESOM
    name: pi
    arrangement: C
    shape: [5839]
    uid: "e548b74fa53eef5ab412c6061330f043"
    url_connectivity: grid/fesom/pi.fesom-conn.ek
"""

ICON_YAML = """
grid_names:
  - icon-grid-0001-r02b05-r: &icon
    description: "Reduced radiation grid (global, R02B05). 80 km resolution."
    type: ICON
    name: icon-grid-0001-r02b05-r
    file: "icon_grid_0001_R02B05_R.nc"
    shape: [81920]
    icon_type: hrz_radiation
    arrangement: C
"""

ICON_CH_YAML = """
grid_names:
  - ICON-CH1-V1: &ch
    type: ICON
    name: ICON-CH1
    shape: [1147980]
    bounding_box: [50.6, -0.9, 41.9, 17.8]
    arrangement: C

  - ICON-CH1: *ch
"""

# An entry can be another name for an entry listed before it: the anchor is on the
# empty value of the target name, so the alias item holds nothing but its name.
ALIAS_YAML = """
grid_names:
  - BASE: &base
    type: ICON
    name: base
    shape: [1147980]
    bounding_box: [50.6, -0.9, 41.9, 17.8]

  - SECOND: *base

  - THIRD: &third
    type: ICON
    description: "Third grid. Fine resolution."
    shape: [20]

  - FOURTH: *third
"""

# Nested layout, with an anchor on the mapping of the fields.
NESTED_ALIAS_YAML = """
grid_names:
  - LAEA-A: &a
      x: [0, 10, 1]
      projection: {type: laea}

  - LAEA-B: *a
"""

# Aliases which cannot be resolved: the anchor of an identifier list, and an entry
# whose target has no fields of its own.
BROKEN_ALIAS_YAML = """
grid_uids:
  - "aaaa1111": &anchor

grid_names:
  - NO-DATA: *anchor

  - EMPTY: &empty

  - CHAIN: *empty
"""

AREA_YAML = """
area_names:
  - globe-prime-meridian:
      type: bounding_box
      bounding_box: [90., 0., -90., 360.]
"""

PROJECTION_YAML = """
projection_names:
    - swisslv95:
          type: proj
          source: "EPSG:2056"
"""

WEIRD_YAML = """
grid_names:
  - odd-one:
      type: whatever
      extra: [1, 2]

  - odd-two:
      type: other
"""

EXTRA_YAML = """
grid_names:
  - spare-five:
      type: yet-other

  - spare-six:
      type: last
"""

BOOLEAN_YAML = """
area_names:
  - false:
      type: bounding_box
      bounding_box: [1., 2., 3., 4.]

  - true:
      type: circle
      centre: [1., 2.]
"""


def _headings(text):
    """The title lines of the generated RST: each is followed by dashes."""
    lines = text.splitlines()
    return [lines[i - 1] for i, line in enumerate(lines) if line and set(line) == {"-"} and i > 0]


def _catalogue(tmp_path, files, name=None):
    """Write a directory of catalogue files and return its path."""
    directory = tmp_path if name is None else tmp_path / name
    directory.mkdir(parents=True, exist_ok=True)

    for filename, text in files.items():
        (directory / filename).write_text(text.lstrip("\n"), encoding="utf-8")

    return directory


@pytest.fixture
def check_rst():
    """Return a function parsing RST with docutils.

    It gives back the number of tables and the error level messages, with the
    docutils ``list-table`` directive registered so that the tables are really
    parsed.
    """
    nodes = pytest.importorskip("docutils.nodes")
    core = pytest.importorskip("docutils.core")
    directives = pytest.importorskip("docutils.parsers.rst.directives")

    try:
        from docutils.parsers.rst.directives.tables import ListTable
        from docutils.utils import Reporter
    except ImportError:  # pragma: no cover
        pytest.skip("docutils has no list-table directive")

    registry = getattr(directives, "_directives", {})
    saved = dict(registry)
    directives.register_directive("list-table", ListTable)

    def check(text):
        doctree = core.publish_doctree(
            text,
            settings_overrides={"report_level": 1, "halt_level": 6, "warning_stream": io.StringIO()},
        )
        errors = []
        for message in doctree.findall(nodes.system_message):
            if (message.get("level") or 0) >= Reporter.ERROR_LEVEL:
                errors.append(message.astext().strip())

        tables = len(list(doctree.findall(nodes.table)))
        return tables, errors

    try:
        yield check
    finally:
        registry.clear()
        registry.update(saved)


def test_entry_layouts(tmp_path):
    d = _catalogue(tmp_path, {"grid.yaml": GRID_YAML, "ORCA.yaml": ORCA_YAML})
    out = mod.render("grids", path=d)

    # section titles describe the entries, never the files
    assert "\nRegional projected grids\n------------------------\n" in out
    assert "\nORCA ocean grids\n----------------\n" in out
    assert "grid.yaml" not in out
    assert "ORCA.yaml" not in out

    # nested layout: name, and the projection reduced to its type
    assert "``LAEA-EFAS-5km``" in out
    assert "``swisslv95``" in out
    assert "``laea``" in out
    assert "``[0, 4995000, 5000]``" in out

    # the false key of grid.yaml does not crash it and is not shown as False
    assert "False" not in out
    assert "1818500" not in out

    # flat layout: the name is the key holding nothing, the other keys are fields
    assert "``ORCA2_T``" in out
    assert "``ORCA1_F``" in out
    assert "``ORCA2``" in out  # the model, from the name field
    assert "``T``" in out  # the arrangement
    assert "``[149, 182]``" in out  # the shape

    # an entry which is an alias of another one is listed with the data it refers to
    assert "``ORCA2_LIM`` (alias of ``ORCA2_T``)" in out

    # uids, urls and identifier lists are not mentioned at all
    assert "d5bde4f52ff3a9bea5629cd9ac514410" not in out
    assert "example.org" not in out
    assert "grid_uids" not in out
    assert "identifiers" not in out

    assert out.count(".. list-table::") == 2
    assert "The table lists 2 grids." in out
    assert "The table lists 3 grids." in out


def test_single_table_has_no_title(tmp_path):
    d = _catalogue(tmp_path, {"area.yaml": AREA_YAML})
    out = mod.render("areas", path=d)

    # one table only: the page title is enough, no section title is printed
    assert ".. list-table::" in out
    assert "Named areas which can be used in place of an explicit area specification." in out
    assert "The table lists 1 area." in out
    assert _headings(out) == []


def test_hydro_columns(tmp_path):
    d = _catalogue(tmp_path, {"hydro.yaml": HYDRO_YAML})
    out = mod.render("grids", path=d)

    assert "   * - Name\n     - Type\n     - n\n     - Area\n" in out
    assert "``EFAS-5``" in out
    assert "``arakawa_c``" in out
    assert "``5400``" in out
    assert "``[72.25, -25.25, 22.75, 50.25]``" in out
    assert "``west = -180``" in out


def test_orca_columns(tmp_path):
    d = _catalogue(tmp_path, {"ORCA.yaml": ORCA_YAML})
    out = mod.render("grids", path=d)

    assert "   * - Name\n     - Model\n     - Arrangement\n     - Shape\n" in out
    assert "``ORCA2_T``" in out
    assert "``[149, 182]``" in out


def test_fesom_is_excluded(tmp_path):
    d = _catalogue(tmp_path, {"grid.yaml": GRID_YAML, "FESOM.yaml": FESOM_YAML})

    # neither when its kind is asked for nor in a combined listing
    for out in (mod.render("grids", path=d), mod.render(path=d)):
        assert "``LAEA-EFAS-5km``" in out
        assert "pi_C" not in out
        assert "fesom" not in out.lower()


def test_icon_columns(tmp_path):
    d = _catalogue(tmp_path, {"ICON-DWD.yaml": ICON_YAML})
    out = mod.render("grids", path=d)

    # the description is shown when present, the type when there is one
    assert "   * - Name\n     - Description\n     - Shape\n     - Type\n" in out
    assert "``icon-grid-0001-r02b05-r``" in out
    assert "``hrz_radiation``" in out
    assert "``[81920]``" in out
    assert "Reduced radiation grid (global, R02B05). 80 km resolution." in out
    assert "icon_grid_0001_R02B05_R.nc" not in out

    # a file without descriptions has no description column, and the icon type
    # column falls back to the bounding box
    d = _catalogue(tmp_path, {"ICON-CH.yaml": ICON_CH_YAML}, name="no-description")
    out = mod.render("grids", path=d)
    assert "   * - Name\n     - Shape\n     - Bounding box\n" in out
    assert "Description" not in out
    assert "``ICON-CH1`` (alias of ``ICON-CH1-V1``)" in out
    assert out.count("``[1147980]``") == 2
    assert "``[50.6, -0.9, 41.9, 17.8]``" in out


def test_titles_per_kind(tmp_path):
    files = {"grid.yaml": GRID_YAML, "ORCA.yaml": ORCA_YAML, "ICON-CH.yaml": ICON_CH_YAML, "area.yaml": AREA_YAML}
    d = _catalogue(tmp_path, files)

    out = mod.render("grids", path=d)

    # known catalogues have a fixed title, and sections are ordered by file name
    headings = _headings(out)
    assert headings == ["Regional projected grids", "ICON grids of MeteoSwiss", "ORCA ocean grids"], headings

    # the entries of the other kinds are not printed here
    assert "globe-prime-meridian" not in out

    # without arguments all the kinds are printed, grids first
    out = mod.render(path=d)
    assert "``globe-prime-meridian``" in out
    assert out.index("Regional projected grids") < out.index("globe-prime-meridian")


def test_alias_entries(tmp_path):
    d = _catalogue(tmp_path, {"ICON-CH.yaml": ALIAS_YAML})
    out = mod.render("grids", path=d)

    # the alias rows show the data of the entry they refer to
    assert "``SECOND`` (alias of ``BASE``)" in out
    assert "``FOURTH`` (alias of ``THIRD``)" in out
    assert out.count("``[1147980]``") == 2
    assert out.count("Third grid. Fine resolution.") == 2

    # the plain entries are not marked
    assert "``BASE``\n" in out
    assert "``THIRD``\n" in out
    assert out.count("alias of") == 2

    # a row which is an alias still counts as one entry
    assert "The table lists 4 grids." in out


def test_nested_layout_alias(tmp_path):
    d = _catalogue(tmp_path, {"grid.yaml": NESTED_ALIAS_YAML})
    out = mod.render("grids", path=d)

    assert "``LAEA-B`` (alias of ``LAEA-A``)" in out
    assert out.count("``[0, 10, 1]``") == 2
    assert out.count("``laea``") == 2


def test_unresolvable_alias(tmp_path, check_rst):
    d = _catalogue(tmp_path, {"ICON-BROKEN.yaml": BROKEN_ALIAS_YAML})
    out = mod.render("grids", path=d)

    # an alias whose data cannot be found is shown as it is, without a mark
    assert "``NO-DATA``" in out
    assert "``CHAIN``" in out
    assert "alias of" not in out
    assert "The table lists 3 grids." in out

    # the unknown file joins the merged table, which is alone and titled by
    # content only: the stem never shows
    assert _headings(out) == []
    assert "broken" not in out.lower()

    # the identifier list is not mentioned
    assert "grid_uids" not in out
    assert "aaaa1111" not in out

    tables, errors = check_rst(out)
    assert errors == []
    assert tables == 1


def test_no_alias_is_unchanged(tmp_path):
    d = _catalogue(
        tmp_path,
        {"grid.yaml": GRID_YAML, "ORCA.yaml": ORCA_YAML.replace("  - ORCA2_LIM: *ORCA2_T\n", "")},
    )
    out = mod.render("grids", path=d)

    assert "alias" not in out


def test_area_projection_columns(tmp_path):
    d = _catalogue(tmp_path, {"area.yaml": AREA_YAML, "projection.yaml": PROJECTION_YAML})
    out = mod.render(path=d)

    assert "   * - Name\n     - Type\n     - Details\n" in out
    assert "``globe-prime-meridian``" in out
    assert "``bounding_box = [90.0, 0.0, -90.0, 360.0]``" in out
    assert "``source = EPSG:2056``" in out
    assert "The table lists 1 area." in out
    assert "The table lists 1 projection." in out


def test_unknown_files_are_merged(tmp_path, check_rst):
    d = _catalogue(tmp_path, {"ORCA.yaml": ORCA_YAML, "weird.yaml": WEIRD_YAML, "extra.yaml": EXTRA_YAML})
    out = mod.render("grids", path=d)

    # the two unknown catalogues share one table at the end, titled by content
    assert _headings(out) == ["ORCA ocean grids", "Additional grids"], _headings(out)

    # no file name anywhere in the listing (the field named "extra" is content)
    assert "weird" not in out.lower()
    assert "yaml" not in out.lower()

    # the merged table is generic: name, type and the other fields as details
    assert "   * - Name\n     - Type\n     - Details\n" in out
    assert "Further predefined grids." in out
    assert "``odd-one``" in out
    assert "``whatever``" in out
    assert "``extra = [1, 2]``" in out
    assert "``spare-five``" in out
    assert "``spare-six``" in out

    # one table for both files, with one summary counting all its rows
    assert out.count(".. list-table::") == 2
    assert "The table lists 4 grids." in out
    assert "The table lists 3 grids." in out  # the known catalogue next to it

    tables, errors = check_rst(out)
    assert errors == []
    assert tables == 2


def test_merged_table_alone_has_no_title(tmp_path):
    d = _catalogue(tmp_path, {"weird.yaml": WEIRD_YAML})
    out = mod.render("grids", path=d)

    # the merged table is the only one of its kind: the single-table rule holds
    assert _headings(out) == []
    assert "Further predefined grids." in out
    assert "   * - Name\n     - Type\n     - Details\n" in out
    assert "The table lists 2 grids." in out


def test_boolean_keys(tmp_path):
    d = _catalogue(tmp_path, {"area.yaml": BOOLEAN_YAML})
    out = mod.render("areas", path=d)

    # a "false" key is loaded as the boolean False, and must print as false
    assert "``false``" in out
    assert "``true``" in out
    assert "False" not in out
    assert "The table lists 2 areas." in out


def test_unknown_kind(tmp_path):
    d = _catalogue(tmp_path, {"grid.yaml": GRID_YAML, "area.yaml": AREA_YAML})

    out = mod.render("shapefile", path=d)
    assert ".. note::" in out
    assert "``shapefile``" in out
    assert "``grids``" in out  # the available kinds are listed
    assert ".. list-table::" not in out

    # the kinds which do exist are still shown
    out = mod.render("grids", "shapefile", path=d)
    assert "``LAEA-EFAS-5km``" in out
    assert "globe-prime-meridian" not in out
    assert ".. note::" in out
    assert "``shapefile``" in out


def test_missing_directory(tmp_path, capsys):
    out = mod.render(path=tmp_path / "nowhere")
    assert ".. note::" in out
    assert "No catalogue directory" in out

    mod.execute(path=tmp_path / "nowhere")
    assert ".. note::" in capsys.readouterr().out


def test_empty_directory(tmp_path):
    out = mod.render(path=tmp_path)
    assert ".. note::" in out
    assert "no catalogue files" in out

    out = mod.render("grids", path=tmp_path)
    assert ".. note::" in out
    assert "no catalogue files" in out


def test_unreadable_file(tmp_path, check_rst):
    broken = "grid_names:\n  - name: [unclosed\n invalid: :\n"
    files = {"grid.yaml": GRID_YAML, "broken.yaml": broken, "empty.yaml": "", "list.yaml": "- a\n- b\n"}
    d = _catalogue(tmp_path, files)

    # a combined listing notes the files it could not group
    out = mod.render(path=d)

    assert "``broken.yaml`` could not be read" in out
    assert "``empty.yaml`` lists no grids, areas or projections" in out
    assert "``list.yaml`` lists no grids, areas or projections" in out
    assert "``LAEA-EFAS-5km``" in out

    # a kind listing leaves them out, and shows no file names
    out = mod.render("grids", path=d)
    assert "broken" not in out
    assert "yaml" not in out.lower()
    assert "``LAEA-EFAS-5km``" in out

    tables, errors = check_rst(mod.render(path=d))
    assert errors == []
    assert tables == 1


def test_execute_prints(tmp_path, monkeypatch, capsys):
    _catalogue(tmp_path, {"grid.yaml": GRID_YAML, "area.yaml": AREA_YAML})
    monkeypatch.setattr(mod, "CATALOGUE_DIR", tmp_path)

    mod.execute()
    out = capsys.readouterr().out

    assert ".. list-table::" in out
    assert "   :header-rows: 1" in out
    assert "``LAEA-EFAS-5km``" in out
    assert "``globe-prime-meridian``" in out

    # the arguments select the kind
    mod.execute("areas")
    out = capsys.readouterr().out
    assert "``globe-prime-meridian``" in out
    assert "LAEA-EFAS-5km" not in out


def test_rst_valid(tmp_path, check_rst):
    files = {
        "grid.yaml": GRID_YAML,
        "ORCA.yaml": ORCA_YAML,
        "ICON-DWD.yaml": ICON_YAML,
        "ICON-CH.yaml": ICON_CH_YAML,
        "ICON-X.yaml": ALIAS_YAML,
    }
    d = _catalogue(tmp_path, files)
    out = mod.render("grids", path=d)

    tables, errors = check_rst(out)
    assert errors == []
    # ICON-X is unknown: four titled sections plus one merged table
    assert tables == len(files)

    # notes are valid RST too
    tables, errors = check_rst(mod.render(path=tmp_path / "nowhere"))
    assert errors == []
    assert tables == 0

    tables, errors = check_rst(mod.render("nope", path=d))
    assert errors == []


def test_installed_catalogues():
    pytest.importorskip("eckitlib")

    grids = mod.render("grids")
    areas = mod.render("areas")
    projections = mod.render("projections")
    out = mod.render()

    # names from both entry layouts, on the page of their kind
    for name in ("LAEA-EFAS-5km", "GLOFAS-4", "ORCA2_T", "icon-grid-0001-r02b05-r", "ICON-CH1", "swisslv95"):
        assert f"``{name}``" in grids, name
    assert "``globe-prime-meridian``" in areas
    assert "``globe-antimeridian``" in areas
    assert "``swisslv95``" in projections
    assert "GLOFAS" not in areas
    assert "GLOFAS" not in projections

    # FESOM is not shown anywhere
    assert "pi_C" not in out
    assert "fesom" not in out.lower()

    # the six grid sections, titled by what they define, and no file names
    assert _headings(grids) == list(GRID_SECTIONS)
    assert _headings(areas) == []
    assert _headings(projections) == []
    for page in (grids, areas, projections):
        for heading in _headings(page):
            assert "yaml" not in heading.lower()

    # the identifier lists are not mentioned
    assert "grid_uids" not in grids

    # the false key of grid.yaml is a name, not a Python boolean
    assert "False" not in grids

    # ICON-CH lists two entries as aliases of others
    assert "``ICON-CH1`` (alias of ``ICON-CH1-V1``)" in grids
    assert "``ICON-CH2`` (alias of ``ICON-CH2-V1``)" in grids

    # a row which is not an alias stays as it is
    assert "``ORCA2_T``\n" in grids
    assert "``LAEA-EFAS-5km``\n" in grids

    assert out.count(".. list-table::") == grids.count(".. list-table::") + 2


def test_installed_rst_valid(check_rst):
    pytest.importorskip("eckitlib")

    for out in (mod.render("grids"), mod.render("areas"), mod.render("projections"), mod.render()):
        tables, errors = check_rst(out)
        assert errors == []

    assert mod.render("grids").count(".. list-table::") == len(GRID_SECTIONS)
