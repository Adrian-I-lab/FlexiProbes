"""Multi-strain consensus + per-strain binding scoring."""
from flexprobe.consensus import (
    StrainAlignment, build_consensus, conserved_start_positions,
    resolve_to_reference, score_probe_against_all_strains,
    strain_unique_start_positions,
)
import flexprobe.consensus as C


def test_external_msa_returns_none_without_tools(monkeypatch):
    # Audit #25: with no MAFFT/MUSCLE anywhere, _external_msa yields None so the
    # caller falls back to the built-in progressive aligner. Patched at
    # _find_msa_tool rather than shutil.which: since the resolver also searches
    # the common install dirs, a disabled which() no longer means "no tool".
    monkeypatch.setattr(C, "_find_msa_tool", lambda _n: None)
    assert C._external_msa([("a", "ACGT"), ("b", "ACGTA")]) is None


def test_find_msa_tool_looks_beyond_path(monkeypatch, tmp_path):
    # The bug this guards: a macOS .app launched from Finder inherits only a
    # minimal PATH, so shutil.which() misses a Homebrew MAFFT and the run drops
    # silently to the built-in aligner. The resolver must search the common
    # install dirs too. Both branches are asserted so the check is shown taking
    # two different values on inputs we control.
    import shutil
    monkeypatch.setattr(shutil, "which", lambda _n: None)   # stripped PATH
    stub = tmp_path / "mafft"
    stub.write_text("#!/bin/sh\nexit 0\n")
    stub.chmod(0o755)

    monkeypatch.setattr(C, "_COMMON_MSA_DIRS", (str(tmp_path),))
    assert C._find_msa_tool("mafft") == str(stub)            # found off-PATH

    monkeypatch.setattr(C, "_COMMON_MSA_DIRS", (str(tmp_path / "nope"),))
    monkeypatch.delenv("FLEXPROBE_MAFFT_DIR", raising=False)
    assert C._find_msa_tool("mafft") is None                 # and absent when gone


def test_find_msa_tool_honours_env_override(monkeypatch, tmp_path):
    # $FLEXPROBE_MAFFT_DIR is the explicit escape hatch, mirroring
    # $FLEXPROBE_BLAST_DIR. It must win even when which() finds nothing.
    import shutil
    monkeypatch.setattr(shutil, "which", lambda _n: None)
    monkeypatch.setattr(C, "_COMMON_MSA_DIRS", ())
    stub = tmp_path / "mafft"
    stub.write_text("#!/bin/sh\nexit 0\n")
    stub.chmod(0o755)
    monkeypatch.setenv("FLEXPROBE_MAFFT_DIR", str(tmp_path))
    assert C._find_msa_tool("mafft") == str(stub)


def test_external_msa_reads_tool_output_in_input_order(monkeypatch):
    # Simulate MAFFT being present and writing an aligned FASTA (reordered);
    # _external_msa must restore input order and equal length, and report the
    # tool it used so the run's provenance survives into the outputs.
    import subprocess
    import types
    monkeypatch.setattr(C, "_find_msa_tool",
                        lambda n: "/usr/bin/mafft" if n == "mafft" else None)

    def fake_run(cmd, stdout=None, **kw):
        if "--version" in cmd:
            return types.SimpleNamespace(stdout="v7.526 (2024/Apr/26)",
                                         stderr="")
        # tool reorders records (idx 1 then 0); we must restore input order
        stdout.write(">1\nAC-GT\n>0\nACGGT\n")
        return None
    monkeypatch.setattr(subprocess, "run", fake_run)
    aligned, tool = C._external_msa([("a", "ACGGT"), ("b", "ACGT")])
    assert aligned == ["ACGGT", "AC-GT"]       # input order, equal length
    assert tool == "MAFFT v7.526 (2024/Apr/26)"


def test_build_consensus_records_which_aligner_ran(monkeypatch):
    # A finished result must say which aligner produced its alignment, MAFFT
    # or the built-in fallback, which is NOT a true multiple alignment. A
    # transient log line is not enough.
    monkeypatch.setattr(C, "_find_msa_tool", lambda _n: None)
    # Lengths must differ by MORE than length_tolerance (5) to reach an
    # aligner at all; 16 vs 6 does, 10 vs 6 would take the trim path.
    aln = build_consensus([("a", "ACGTACGTACGTACGT"), ("b", "ACGTAC")])
    assert aln.aligner == "built-in progressive"
    assert aln.used_pairwise is True
    # Near-equal lengths take the trim path, where no alignment is performed.
    same = build_consensus([("a", "ACGTACGTAC"), ("b", "ACGTACGTAA")])
    assert same.aligner == "trim-to-shortest"
    assert same.used_pairwise is False


def test_resolve_to_reference_fills_ambiguous_bases():
    # Any N / IUPAC code defaults to the first (reference) sequence's base, so
    # the resulting probe window is concrete ACGT with no N's (#3).
    assert resolve_to_reference("ACNTACGTSC", "ACGTACGTACGT") == "ACGTACGTAC"
    # Already-clean windows pass through unchanged (just uppercased).
    assert resolve_to_reference("acgt", "TTTT") == "ACGT"
    # If the reference is itself ambiguous at a position, keep the original.
    assert resolve_to_reference("ANGT", "ANGT") == "ANGT"


def test_consensus_iupac_at_variant():
    a = "ACGTACGTAC"
    b = "ACGTACGTAG"   # differs at last position (C vs G -> S)
    aln = build_consensus([("A", a), ("B", b)])
    assert aln.consensus[:9] == "ACGTACGTA"
    assert aln.consensus[9] == "S"          # C/G ambiguity
    assert aln.variant_positions == [9]


def test_conserved_windows_exclude_variants():
    a = "ACGT" * 20            # 80 bp
    b = a[:40] + "T" + a[41:]  # one SNP at position 40
    aln = build_consensus([("A", a), ("B", b)])
    starts = conserved_start_positions(aln, tile_length=50, tile_step=1)
    # No conserved 50-mer can contain position 40.
    assert all(not (s <= 40 < s + 50) for s in starts)


def test_strain_unique_requires_min_mismatches():
    a = "ACGT" * 20
    # strain B differs by 10 bases, positioned so that in the start=0 window
    # they sit at the LHS/nick boundary (offsets 15-24) rather than out at the
    # window's far edge — otherwise the junction-duplex veto
    # (see test_strain_unique_per_half_distal_mismatches_still_seal_the_duplex
    # below) correctly excludes every window, since 20+ contiguous bases would
    # still span the nick untouched. This fixture is meant to test the plain
    # count threshold, not that veto, so the mismatches are placed where they
    # also break the duplex.
    blist = list(a)
    for i in range(15, 25):
        blist[i] = "T" if blist[i] != "T" else "A"
    b = "".join(blist)
    aln = build_consensus([("A", a), ("B", b)])
    starts = strain_unique_start_positions(
        aln, target_strain_idx=0, tile_length=50, tile_step=1,
        min_mismatches_vs_others=5,
    )
    assert starts          # at least one window distinguishes strain A
    assert all(mm >= 5 for _s, mm in starts)


# --- per-half strain discrimination ---------------------------------------
# The >=5 threshold is borrowed from the 10x off-target rule, which scopes it
# PER PROBE HALF ("at least five mismatches in at least one of the LHS or RHS
# probes", CG000621 Rev D / CG000839 Rev A). Counting it over the whole 50-mer
# accepts windows the engine's own off-target classifier says would still
# hybridise and ligate on the non-target strain.

_REF_50 = ("ACGT" * 13)[:50]


def _mutate(seq, positions):
    """Return ``seq`` with every index in ``positions`` changed to another base."""
    out = list(seq)
    for p in positions:
        out[p] = "T" if out[p] != "T" else "A"
    return "".join(out)


def _one_window(ref, other):
    """Alignment of two equal-length 50-mers — exactly one 50 bp window at 0."""
    aln = build_consensus([("A", ref), ("B", other)])
    # Guard the fixture: if the aligner introduced gaps or extra columns the
    # mismatch positions below would no longer mean what the test claims.
    assert aln.length == 50
    assert all("-" not in s for s in aln.aligned_seqs)
    return aln


def test_strain_unique_per_half_rejects_junction_split_mismatches():
    # The failure mode measured on the OROV panel: 3 LHS + 2 RHS mismatches sum
    # to 5, but neither half reaches 5, so both halves still bind the other
    # strain and the pair ligates there. Not strain-discriminating.
    other = _mutate(_REF_50, [2, 7, 13, 30, 41])       # 3 in LHS, 2 in RHS
    aln = _one_window(_REF_50, other)
    common = dict(target_strain_idx=0, tile_length=50, tile_step=1,
                  min_mismatches_vs_others=5)

    # Legacy whole-window count accepts it (3 + 2 == 5)...
    legacy = strain_unique_start_positions(aln, per_half=False, **common)
    assert legacy == [(0, 5)]

    # ...the per-half count does not (max(3, 2) == 3).
    assert strain_unique_start_positions(aln, per_half=True, **common) == []


def test_strain_unique_per_half_is_the_default():
    # Not a passthrough: omitting per_half must give the per-half answer, not
    # the legacy one. Same fixture as above.
    other = _mutate(_REF_50, [2, 7, 13, 30, 41])
    aln = _one_window(_REF_50, other)
    assert strain_unique_start_positions(
        aln, target_strain_idx=0, tile_length=50, tile_step=1,
        min_mismatches_vs_others=5,
    ) == []


def test_strain_unique_per_half_keeps_mismatches_within_one_half():
    # 5 mismatches all inside the LHS: the LHS cannot bind the other strain, so
    # no ligation product forms. Qualifies under both readings.
    other = _mutate(_REF_50, [1, 5, 9, 14, 20])        # 5 LHS, 0 RHS
    aln = _one_window(_REF_50, other)
    common = dict(target_strain_idx=0, tile_length=50, tile_step=1,
                  min_mismatches_vs_others=5)
    assert strain_unique_start_positions(aln, per_half=True, **common) == [(0, 5)]
    assert strain_unique_start_positions(aln, per_half=False, **common) == [(0, 5)]


def test_strain_unique_reported_mm_is_in_the_gate_unit():
    # The second tuple element ranks windows, so it must be the quantity the
    # gate compared: best-half count in per-half mode, whole-window otherwise.
    other = _mutate(_REF_50, [1, 5, 9, 14, 20, 22, 30, 41, 47])   # 6 LHS, 3 RHS
    aln = _one_window(_REF_50, other)
    common = dict(target_strain_idx=0, tile_length=50, tile_step=1,
                  min_mismatches_vs_others=5)
    assert strain_unique_start_positions(aln, per_half=True, **common) == [(0, 6)]
    assert strain_unique_start_positions(aln, per_half=False, **common) == [(0, 9)]


def test_strain_unique_per_half_follows_tile_size():
    # The half boundary must track the tile length, not a hard-coded 25 — the
    # same invariant score_probe_against_all_strains holds (audit H1). On a
    # 40-mer the halves are 20 bp, so positions 20..24 are RHS, not LHS.
    ref = ("ACGT" * 10)[:40]
    other = _mutate(ref, [16, 17, 18, 21, 22])     # 40-mer: 3 LHS + 2 RHS
    aln = build_consensus([("A", ref), ("B", other)])
    assert aln.length == 40 and all("-" not in s for s in aln.aligned_seqs)
    common = dict(target_strain_idx=0, tile_length=40, tile_step=1,
                  min_mismatches_vs_others=5)
    assert strain_unique_start_positions(aln, per_half=True, **common) == []
    assert strain_unique_start_positions(aln, per_half=False, **common) == [(0, 5)]


def test_strain_unique_per_half_agrees_with_reported_half_columns():
    # The invariant the fix exists to restore: every window the gate calls
    # strain-unique must reach the threshold in a half, as measured by the same
    # scorer that fills the lhs_mm__* / rhs_mm__* columns of results.csv.
    ref = "ACGTTGCA" * 20                                  # 160 bp
    other = _mutate(ref, [
        60, 62, 64, 66, 68,              # clustered: 5 fit inside one 25 bp half
        88, 95, 102, 109, 116, 123,      # spread 7 bp apart: >=5 total, <=4 per half
    ])
    aln = build_consensus([("A", ref), ("B", other)])
    assert aln.length == 160 and all("-" not in s for s in aln.aligned_seqs)
    starts = strain_unique_start_positions(
        aln, target_strain_idx=0, tile_length=50, tile_step=1,
        min_mismatches_vs_others=5, per_half=True,
    )
    assert starts                                          # fixture is useful
    target = aln.aligned_seqs[0]
    for start, _mm in starts:
        score = {s.strain_id: s for s in score_probe_against_all_strains(
            aln, start, tile_length=50,
            probe_seq=target[start:start + 50])}["B"]
        assert max(score.lhs_mismatches, score.rhs_mismatches) >= 5

    # And the converse: the legacy gate admits windows that fail exactly that
    # check — i.e. this fixture actually exercises the difference.
    legacy = strain_unique_start_positions(
        aln, target_strain_idx=0, tile_length=50, tile_step=1,
        min_mismatches_vs_others=5, per_half=False,
    )
    assert len(legacy) > len(starts)


def test_strain_unique_per_half_counts_gap_columns_not_shortcuts_them():
    # A gap column is one position the probe cannot base-pair, so it counts as
    # one mismatch — the same thing score_probe_against_all_strains does. The
    # legacy mode instead treats ANY gap as automatically discriminating, which
    # scores a 2 nt indel as if it were a whole mismatched window.
    #
    # Real case from the OROV Large panel (R08, alignment column 56): the other
    # strain carries a 2 bp deletion at the start of the window plus 1 LHS and
    # 1 RHS SNP — 3 LHS / 1 RHS, nowhere near 5 in a half, yet the shortcut
    # admitted it as strain-unique.
    ref = "TGCTCAACCAATATCGGAATAGGATATTGCACTGCCGTGAACCTGAGATA"
    other = "--CTCAACCAATATCGAAATAGGATATTGCACTGCCGTGAACCTGAAATA"
    aln = StrainAlignment(
        fragment_name="t",
        strain_ids=["A", "B"],
        aligned_seqs=[ref, other],
        consensus=ref,
        variant_positions=[],
    )
    common = dict(target_strain_idx=0, tile_length=50, tile_step=1,
                  min_mismatches_vs_others=5)
    # 2 gap columns + 1 SNP in the LHS, 1 SNP in the RHS -> max half is 3.
    assert strain_unique_start_positions(aln, per_half=True, **common) == []
    # The legacy shortcut passes it on the strength of the gap alone.
    assert strain_unique_start_positions(aln, per_half=False, **common) == [(0, 50)]


def test_strain_unique_per_half_still_passes_a_large_deletion():
    # The shortcut is not needed to handle real absence: a deletion long enough
    # to matter produces >=5 uncovered positions in a half on its own. Placed
    # at offsets 17-24, adjacent to the nick, so it also breaks the junction
    # duplex. An 8 bp deletion clear at the FAR end of the
    # half would leave 17+ contiguous matched bases spanning the nick, which is
    # the over-certifying case that veto exists to catch (see
    # test_strain_unique_per_half_distal_mismatches_still_seal_the_duplex).
    ref = "TGCTCAACCAATATCGGAATAGGATATTGCACTGCCGTGAACCTGAGATA"
    other = ref[:17] + "-" * 8 + ref[25:]
    aln = StrainAlignment(
        fragment_name="t",
        strain_ids=["A", "B"],
        aligned_seqs=[ref, other],
        consensus=ref,
        variant_positions=[],
    )
    assert strain_unique_start_positions(
        aln, target_strain_idx=0, tile_length=50, tile_step=1,
        min_mismatches_vs_others=5, per_half=True,
    ) == [(0, 8)]


# --- within-half position: the junction-duplex veto ------------------------
# Per-half fixes WHICH half the mismatches must fall in, but says nothing
# about WHERE in that half. Ligation only cares about the nick, so a raw
# per-half count cannot tell a mismatch that breaks the nick-spanning duplex
# from one that leaves 20 clean bases across it untouched. These two tests are
# a pair a count-only gate cannot distinguish: IDENTICAL mismatch count (5),
# opposite position, opposite verdict.
#
# Direction, stated plainly and matching this codebase's own physics
# elsewhere (constants.mismatch_position_weight, the graded ΔTm ramp;
# specificity.junction_clean_duplex / OFF_TARGET_JUNCTION_DUPLEX_BP, the
# off-target twin this reuses): a mismatch AT the nick is what actually stops
# the ligase, so nick-adjacent mismatches are what should let the gate certify
# discrimination, not distal ones. A window with all its mismatches at the
# probe's terminal end leaves the nick region a perfectly matched duplex, and
# that is the condition under which the OTHER strain still hybridises and
# ligates there regardless of the raw count — see
# test_strain_unique_per_half_distal_mismatches_still_seal_the_duplex.

def test_strain_unique_per_half_distal_mismatches_still_seal_the_duplex():
    # 5 mismatches — meets min_mismatches_vs_others — but all bunched at the
    # FAR/outer end of the LHS (offsets 0-4), 20 bp from the nick. The other
    # strain's window therefore presents an unbroken 20 bp matched run on the
    # LHS side of the nick and a full 25 bp on the RHS side — comfortably
    # >= OFF_TARGET_JUNCTION_DUPLEX_BP (9) on BOTH sides, so the duplex seals
    # and the pair would still hybridise + ligate on strain B. A count-only
    # gate would certify this window strain-unique; the veto does not.
    other = _mutate(_REF_50, [0, 1, 2, 3, 4])
    aln = _one_window(_REF_50, other)
    common = dict(target_strain_idx=0, tile_length=50, tile_step=1,
                  min_mismatches_vs_others=5)
    assert strain_unique_start_positions(aln, per_half=True, **common) == []

    # The escape hatch restores the pre-fix (count-only) admission exactly —
    # same "0 disables" convention as off_target_junction_duplex_bp.
    assert strain_unique_start_positions(
        aln, per_half=True, junction_duplex_bp=0, **common
    ) == [(0, 5)]


def test_strain_unique_per_half_nick_adjacent_mismatches_break_the_duplex():
    # SAME mismatch count (5) as the fixture directly above, but moved to
    # offsets 20-24 — immediately against the nick (offset 24 is the LHS
    # probe's last base). The very first base the scan walks backward from the
    # nick is already mismatched, so clean_left == 0 and no duplex can seal no
    # matter how clean the rest of the half is. Genuinely discriminating, and
    # the gate must still say so — this is the "different position, different
    # behaviour" half of the pair the old gate could not tell apart.
    other = _mutate(_REF_50, [20, 21, 22, 23, 24])
    aln = _one_window(_REF_50, other)
    common = dict(target_strain_idx=0, tile_length=50, tile_step=1,
                  min_mismatches_vs_others=5)
    assert strain_unique_start_positions(aln, per_half=True, **common) == [(0, 5)]


def test_strain_unique_junction_duplex_veto_only_applies_in_per_half_mode():
    # The legacy whole-window reading stays bit-for-bit reproducible (its own
    # docstring promise): a duplex check keyed to a per-half nick has no
    # meaning without the per-half split, so it must not touch that path even
    # though this fixture's mismatches would seal a duplex under per-half
    # scoring too.
    other = _mutate(_REF_50, [0, 1, 2, 3, 4])
    aln = _one_window(_REF_50, other)
    assert strain_unique_start_positions(
        aln, target_strain_idx=0, tile_length=50, tile_step=1,
        min_mismatches_vs_others=5, per_half=False,
    ) == [(0, 5)]


def test_binding_score_counts_mismatches():
    a = "ACGT" * 20
    blist = list(a); blist[30] = "T" if a[30] != "T" else "A"
    b = "".join(blist)
    aln = build_consensus([("A", a), ("B", b)])
    scores = score_probe_against_all_strains(aln, 0, tile_length=50, probe_seq=a[:50])
    by = {s.strain_id: s for s in scores}
    assert by["A"].total_mismatches == 0
    assert by["B"].total_mismatches == 1
    assert by["B"].estimated_dTm < 0   # mismatch penalizes Tm


def test_binding_score_half_length_follows_tile_size():
    # Audit H1: the LHS/RHS split (and junction window) must follow the tile
    # size, not a hard-coded 25. A mismatch at position 22 is in the LHS half of
    # a 50-mer (half=25) but the RHS half of a 40-mer (half=20).
    a = "ACGT" * 20
    blist = list(a); blist[22] = "T" if a[22] != "T" else "A"
    b = "".join(blist)
    aln = build_consensus([("A", a), ("B", b)])
    s50 = {s.strain_id: s for s in score_probe_against_all_strains(
        aln, 0, tile_length=50, probe_seq=a[:50])}["B"]
    assert (s50.lhs_mismatches, s50.rhs_mismatches) == (1, 0)
    s40 = {s.strain_id: s for s in score_probe_against_all_strains(
        aln, 0, tile_length=40, probe_seq=a[:40])}["B"]
    assert (s40.lhs_mismatches, s40.rhs_mismatches) == (0, 1)


def test_progressive_aligner_preserves_gapped_strains():
    # Audit M1: re-lifting earlier strains as a later strain inserts gaps must
    # not drop or shift their bases. The built-in progressive aligner is used
    # for >2 differing-length strains when no MAFFT/MUSCLE is on PATH.
    from flexprobe.consensus import _pairwise_align
    recs = [("A", "ACGTACGT"), ("B", "ACGTTACGT"), ("C", "ACGACGT")]
    aligned = _pairwise_align(recs)
    assert len({len(s) for s in aligned}) == 1          # all columns aligned
    for (_sid, orig), a in zip(recs, aligned):
        assert a.replace("-", "") == orig               # no base dropped/shifted


def test_iupac_consensus_not_a_mismatch():
    # A consensus 'S' (C/G) should not be a mismatch against either C or G.
    a = "ACGTACGTAC"
    b = "ACGTACGTAG"
    aln = build_consensus([("A", a), ("B", b)])
    scores = score_probe_against_all_strains(aln, 0, tile_length=10)  # use consensus
    for s in scores:
        assert s.total_mismatches == 0


# ---- bundled MAFFT on Windows (simulated) ----
#
# Every Windows run recorded "built-in progressive" while every macOS run used
# the bundled MAFFT. The driver was launched as ["bash", mafft.bat]; Windows has
# no bash on PATH, the OSError was caught, and the run fell back silently. The
# vendored Windows package's own launchers (mafft.bat, mafft-signed.ps1) both
# run usr\bin\bash.exe on /usr/bin/mafft after setting MAFFT_BINARIES, TMPDIR
# and a /usr/bin-first PATH, so that is what must be invoked.

def _fake_mafft(seen, fail=None):
    import types

    def fake_run(cmd, stdout=None, env=None, **kw):
        seen.append((list(cmd), dict(env or {})))
        if "--version" in cmd:
            return types.SimpleNamespace(stdout="v7.526 (2024/Apr/26)",
                                         stderr="")
        if fail is not None:
            raise fail
        stdout.write(">0\nACGGT\n>1\nAC-GT\n")
        return None
    return fake_run


def _bundle(monkeypatch, tmp_path):
    root = tmp_path / "mafft"
    driver = root / "mafft.bat"
    monkeypatch.setattr(C, "_bundled_mafft", lambda: (str(driver), str(root)))
    return root, driver


def test_bundled_mafft_runs_through_vendored_bash_on_windows(monkeypatch,
                                                             tmp_path):
    import subprocess
    root, _driver = _bundle(monkeypatch, tmp_path)
    monkeypatch.setattr(C, "_is_windows", lambda: True)
    monkeypatch.setenv("TMP", str(tmp_path))
    monkeypatch.setenv("MAFFT_BINARIES", "/some/other/mafft")
    seen = []
    monkeypatch.setattr(subprocess, "run", _fake_mafft(seen))

    aligned, tool = C._external_msa([("a", "ACGGT"), ("b", "ACGT")])

    bash = str(root / "usr" / "bin" / "bash.exe")
    run_cmd, run_env = seen[0]
    assert run_cmd[:2] == [bash, "/usr/bin/mafft"]
    assert run_cmd[2:4] == ["--auto", "--quiet"]
    assert run_env["MAFFT_BINARIES"] == "/usr/lib/mafft"
    assert run_env["PATH"].startswith("/usr/bin/:")
    assert run_env["TMPDIR"] == run_env["MAFFT_TMPDIR"] == str(tmp_path)
    assert seen[1][0][:2] == [bash, "/usr/bin/mafft"]     # version probe too
    assert aligned == ["ACGGT", "AC-GT"]
    assert tool == "MAFFT v7.526 (2024/Apr/26) (bundled)"


def test_bundled_mafft_command_unchanged_off_windows(monkeypatch, tmp_path):
    # The other value: macOS and Linux keep ["bash", mafft.bat], with any
    # inherited MAFFT_BINARIES removed so it cannot point at another install.
    import subprocess
    _root, driver = _bundle(monkeypatch, tmp_path)
    monkeypatch.setattr(C, "_is_windows", lambda: False)
    monkeypatch.setenv("MAFFT_BINARIES", "/some/other/mafft")
    seen = []
    monkeypatch.setattr(subprocess, "run", _fake_mafft(seen))

    C._external_msa([("a", "ACGGT"), ("b", "ACGT")])

    run_cmd, run_env = seen[0]
    assert run_cmd[:2] == ["bash", str(driver)]
    assert "MAFFT_BINARIES" not in run_env


def test_failed_bundled_mafft_is_recorded_in_consensus_provenance(monkeypatch,
                                                                  tmp_path):
    # A bundled aligner that exists but fails must not reach only a
    # logger.warning. The fallback must survive into the manifest.
    import subprocess
    from flexprobe.pipeline import _consensus_provenance
    _bundle(monkeypatch, tmp_path)
    monkeypatch.setattr(C, "_is_windows", lambda: True)
    monkeypatch.setattr(subprocess, "run", _fake_mafft(
        [], fail=OSError("[WinError 2] The system cannot find the file")))

    aln = build_consensus([("a", "ACGTACGTACGTACGT"), ("b", "ACGTAC")])

    assert aln.aligner == "built-in progressive"
    reason = getattr(aln, "aligner_fallback", "")
    assert "bundled MAFFT" in reason and "WinError 2" in reason
    assert _consensus_provenance(aln).get("aligner_fallback") == reason


def test_no_fallback_key_when_no_aligner_was_found(monkeypatch):
    # Absent, not empty, when nothing failed: no tool at all is not a failure,
    # and an always-present key would change every multi-strain manifest.
    from flexprobe.pipeline import _consensus_provenance
    monkeypatch.setattr(C, "_bundled_mafft", lambda: None)
    monkeypatch.setattr(C, "_find_msa_tool", lambda _n: None)
    aln = build_consensus([("a", "ACGTACGTACGTACGT"), ("b", "ACGTAC")])
    assert aln.aligner == "built-in progressive"
    assert "aligner_fallback" not in _consensus_provenance(aln)
