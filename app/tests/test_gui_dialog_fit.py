"""Dialogs must open tall enough for their own text.

The HPC export and Build BLAST database dialogs set a fixed starting size.
Word-wrapped help text is allowed to squash below its natural height, so a
starting size that is too small makes the HPC dialog's rows overlap and hides
the Build dialog's last help paragraph.
"""
import importlib.util
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    yield app


def _opens_at_least_its_size_hint(dlg, qapp):
    dlg.show()
    qapp.processEvents()
    try:
        assert dlg.height() >= dlg.sizeHint().height(), (
            dlg.height(), dlg.sizeHint().height())
    finally:
        dlg.close()


def test_build_blast_db_dialog_fits_its_text(qapp, tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path))
    from flexprobe_gui.tools_dialog import BuildDBDialog
    _opens_at_least_its_size_hint(BuildDBDialog(), qapp)


def test_hpc_export_dialog_fits_its_text(qapp, tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path))
    from flexprobe.models import RunConfig
    from flexprobe_gui.hpc_dialog import HpcExportDialog
    from flexprobe_gui.run_history import RunRecord
    cfg = RunConfig(fasta_path=str(tmp_path / "in.fa"),
                    output_dir=str(tmp_path / "out"))
    _opens_at_least_its_size_hint(
        HpcExportDialog([(RunRecord(name="r"), cfg, None)]), qapp)
