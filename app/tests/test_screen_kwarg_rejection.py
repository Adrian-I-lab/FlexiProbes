"""The screen constructors must refuse keywords neither engine accepts.

WHY THIS TEST EXISTS. ``make_screen`` passes ONE kwarg set to whichever engine
is selected, so each engine has to tolerate the other's parameters. A bare
``**_ignored`` would buy that tolerance by swallowing everything -- including
typos and settings the caller believed were taking effect.

For example, ``KmerScreen(off_target_min_mismatch=7)`` must not be absorbed
with no error. That one is the dangerous case rather than a cosmetic one, because ``off_target_min_mismatch`` IS a real setting elsewhere
in the codebase, so a caller can reasonably believe the search budget moved
when it did not.

The tests below are deliberately paired. Rejection alone would pass if the
guard refused everything, which would break ``make_screen``; acceptance alone
would pass if the guard did nothing at all. Both directions have to hold.
"""
import pytest

from flexprobe.screen import BlastScreen, KmerScreen, make_screen


class TestRejectsUnknownKeywords:
    """The guard must be able to say no."""

    def test_kmer_rejects_a_real_setting_that_is_not_its_own(self):
        # The measured case. off_target_min_mismatch belongs to the
        # specificity layer, not to the screen, and passing it here never did
        # anything.
        with pytest.raises(TypeError, match="off_target_min_mismatch"):
            KmerScreen(off_target_min_mismatch=7)

    def test_kmer_rejects_outright_nonsense(self):
        with pytest.raises(TypeError, match="not_a_real_parameter"):
            KmerScreen(not_a_real_parameter=1)

    def test_blast_rejects_outright_nonsense(self):
        with pytest.raises(TypeError, match="not_a_real_parameter"):
            BlastScreen(not_a_real_parameter=1)

    def test_the_error_names_what_is_accepted(self):
        # An error that says only "no" sends the reader to the source. This one
        # has to say what would have worked.
        with pytest.raises(TypeError) as exc:
            KmerScreen(not_a_real_parameter=1)
        assert "max_mismatch" in str(exc.value)


class TestStillAcceptsTheOtherEngineParameters:
    """The guard must NOT break the cross-engine tolerance make_screen needs.

    These are the positive controls. If any of them starts raising, every
    ``make_screen`` call breaks, so a guard that passed only the tests above
    would be worse than no guard.
    """

    def test_kmer_tolerates_a_blast_parameter(self):
        KmerScreen(min_pident=90.0)

    def test_blast_tolerates_a_kmer_parameter(self):
        BlastScreen(cache_dir="/tmp/does-not-need-to-exist")

    def test_each_engine_accepts_its_own(self):
        KmerScreen(max_mismatch=3)
        BlastScreen(min_pident=95.0)

    def test_make_screen_passes_a_mixed_set_to_either_engine(self):
        # The actual call shape this tolerance exists for.
        make_screen("kmer", min_pident=90.0, k=12)
        make_screen("blast", cache_dir="/tmp", min_pident=90.0)

    def test_no_keywords_at_all_is_fine(self):
        KmerScreen()
        BlastScreen()
