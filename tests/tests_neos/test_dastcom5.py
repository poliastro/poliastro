import os
from unittest import mock

import numpy as np
import pytest

from poliastro.examples import iss
from poliastro.neos import dastcom5


@mock.patch("poliastro.neos.dastcom5.np.fromfile")
@mock.patch("poliastro.neos.dastcom5.open")
def test_asteroid_db_is_called_with_right_path(mock_open, mock_np_fromfile):
    dastcom5.asteroid_db()
    mock_open.assert_called_with(dastcom5.AST_DB_PATH, "rb")


@mock.patch("poliastro.neos.dastcom5.np.fromfile")
@mock.patch("poliastro.neos.dastcom5.open")
def test_comet_db_is_called_with_right_path(mock_open, mock_np_fromfile):
    dastcom5.comet_db()
    mock_open.assert_called_with(dastcom5.COM_DB_PATH, "rb")


@mock.patch("poliastro.neos.dastcom5.orbit_from_record")
@mock.patch("poliastro.neos.dastcom5.record_from_name")
def test_orbit_from_name(mock_record_from_name, mock_orbit_from_record):
    name = "example_name"
    mock_orbit_from_record.return_value = iss
    mock_record_from_name.return_value = [1]

    assert dastcom5.orbit_from_name(name) == [iss]
    mock_record_from_name.assert_called_with(name)


@mock.patch("poliastro.neos.dastcom5.np.fromfile")
@mock.patch("poliastro.neos.dastcom5.open")
def test_read_headers(mock_open, mock_np_fromfile):
    dastcom5.read_headers()
    mock_open.assert_any_call(
        os.path.join(dastcom5.DBS_LOCAL_PATH, "dast5_le.dat"), "rb"
    )
    mock_open.assert_any_call(
        os.path.join(dastcom5.DBS_LOCAL_PATH, "dcom5_le.dat"), "rb"
    )


@mock.patch("poliastro.neos.dastcom5.read_headers")
@mock.patch("poliastro.neos.dastcom5.np.fromfile")
@mock.patch("poliastro.neos.dastcom5.open")
def test_read_record(mock_open, mock_np_fromfile, mock_read_headers):
    mocked_ast_headers = np.array(
        [(3184, -1, b"00740473", b"00496815")],
        dtype=[
            ("IBIAS1", np.int32),
            ("IBIAS0", np.int32),
            ("ENDPT2", "|S8"),
            ("ENDPT1", "|S8"),
        ],
    )
    mocked_com_headers = np.array([(99999,)], dtype=[("IBIAS2", "<i4")])

    mock_read_headers.return_value = mocked_ast_headers, mocked_com_headers
    dastcom5.read_record(740473)
    mock_open.assert_called_with(
        os.path.join(dastcom5.DBS_LOCAL_PATH, "dast5_le.dat"), "rb"
    )
    dastcom5.read_record(740473 + 1)
    mock_open.assert_called_with(
        os.path.join(dastcom5.DBS_LOCAL_PATH, "dcom5_le.dat"), "rb"
    )


@mock.patch("poliastro.neos.dastcom5.os.makedirs")
@mock.patch("poliastro.neos.dastcom5.zipfile")
@mock.patch("poliastro.neos.dastcom5.os.path.isdir")
@mock.patch("poliastro.neos.dastcom5.urllib.request")
def test_download_dastcom5_raises_error_when_folder_exists(
    mock_request, mock_isdir, mock_zipfile, mock_makedirs
):
    mock_isdir.side_effect = lambda x: x == os.path.join(
        dastcom5.POLIASTRO_LOCAL_PATH, "dastcom5"
    )
    with pytest.raises(FileExistsError):
        dastcom5.download_dastcom5()
    mock_isdir.assert_called_once_with(
        os.path.join(dastcom5.POLIASTRO_LOCAL_PATH, "dastcom5")
    )


@mock.patch("poliastro.neos.dastcom5.urllib.request")
@mock.patch("poliastro.neos.dastcom5.os.makedirs")
@mock.patch("poliastro.neos.dastcom5.zipfile")
@mock.patch("poliastro.neos.dastcom5.os.path.isdir")
def test_download_dastcom5_creates_folder(
    mock_isdir, mock_zipfile, mock_makedirs, mock_request
):
    mock_isdir.return_value = False
    mock_zipfile.is_zipfile.return_value = False
    dastcom5.download_dastcom5()
    mock_makedirs.assert_called_once_with(dastcom5.POLIASTRO_LOCAL_PATH)


@mock.patch("poliastro.neos.dastcom5.zipfile")
@mock.patch("poliastro.neos.dastcom5.os.path.isdir")
@mock.patch("poliastro.neos.dastcom5.urllib.request.urlretrieve")
def test_download_dastcom5_downloads_file(mock_request, mock_isdir, mock_zipfile):
    mock_isdir.side_effect = lambda x: x == dastcom5.POLIASTRO_LOCAL_PATH
    mock_zipfile.is_zipfile.return_value = False
    dastcom5.download_dastcom5()
    mock_request.assert_called_once_with(
        dastcom5.FTP_DB_URL + "dastcom5.zip",
        os.path.join(dastcom5.POLIASTRO_LOCAL_PATH, "dastcom5.zip"),
        dastcom5._show_download_progress,
    )


@mock.patch("poliastro.neos.dastcom5.string_record_from_name")
def test_issue_902(mock_string_record_from_name):
    name = "1P"
    mock_string_record_from_name.return_value = [
        "90000001 Halley",
        "90000002 Halley",
    ]
    assert dastcom5.record_from_name(name) == [
        90000001,
        90000002,
    ]
