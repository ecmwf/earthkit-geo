#!/usr/bin/env python3

# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


import pytest

from earthkit.geo.utils.url import join_url_path


@pytest.mark.parametrize(
    "root,args,expected",
    [
        ("http://example.com", [], "http://example.com"),
        ("http://example.com/", [], "http://example.com/"),
        ("http://example.com", ["path"], "http://example.com/path"),
        ("http://example.com/", ["path"], "http://example.com/path"),
        ("http://example.com", ["path", "to", "file"], "http://example.com/path/to/file"),
        ("http://example.com/", ["path", "to", "file"], "http://example.com/path/to/file"),
    ],
)
def test_join_url_path(root, args, expected):
    assert join_url_path(root, *args) == expected
