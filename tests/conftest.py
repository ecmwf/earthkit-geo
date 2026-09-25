# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import pytest

SKIP = {
    "short": ["long_test"],
    "long": [],
}


def pytest_addoption(parser):
    help_str = "NAME: short, long. Runs a subset of tests.\n"
    for k, v in SKIP.items():
        if v:
            help_str += f"'{k}': skip tests marked as {','.join(v)}.\n"
        else:
            help_str += f"'{k}': do not skip tests.\n"

    parser.addoption(
        "-E",
        action="store",
        metavar="NAME",
        default="short",
        help=help_str,
    )


_ekd_temporary_config = None


def pytest_sessionstart(session):
    global _ekd_temporary_config

    # earthkit-data is an optional dependency: skip if it is not installed.
    try:
        from earthkit.data.core.config import CONFIG as EKD_CONFIG
    except ImportError:
        return

    # Run the whole test session under a temporary earthkit-data config using
    # a temporary cache, so tests never touch (or are affected by) the user's
    # real earthkit-data cache. Also disable the out-of-date check: test
    # fixtures are static, so there is no need to make a network request to
    # check for updates every time a cached file is used.
    _ekd_temporary_config = EKD_CONFIG.temporary({
        "cache-policy": "temporary",
        "check-out-of-date-urls": False,
    })
    _ekd_temporary_config.__enter__()


def pytest_sessionfinish(session, exitstatus):
    global _ekd_temporary_config

    if _ekd_temporary_config is not None:
        _ekd_temporary_config.__exit__(None, None, None)
        _ekd_temporary_config = None


def pytest_runtest_setup(item):
    # print(f"config {item.config.option}")
    flag = item.config.getoption("-E")
    marks_to_skip = SKIP[flag]

    marks_in_items = list([m.name for m in item.iter_markers()])

    if marks_to_skip is None:
        if flag not in marks_in_items:
            pytest.skip(f"test is skipped because custom pytest option : -E {flag}")
        return

    for m in marks_in_items:
        if m in marks_to_skip:
            pytest.skip(f"test is skipped because custom pytest option: -E {flag}")

    if "matrix_db" in marks_in_items:
        from earthkit.geo.grids._regrid.backends.precomputed.db import SYS_DB

        SYS_DB._clear_index()

    tmp_cache = "tmp_cache" in marks_in_items

    # settings
    from earthkit.geo import config

    # ensure settings are not saved automatically
    config.autosave = False

    # ensure all the tests use the default settings
    if tmp_cache:
        # ensure these tests use a temporary cache
        config.reset()
        config.set("cache-policy", "temporary")
    else:
        config.reset()
        config.set("cache-policy", "user")
