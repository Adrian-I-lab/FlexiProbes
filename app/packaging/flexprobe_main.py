"""PyInstaller entry point — launches the flexprobe GUI."""
import multiprocessing


def _launch() -> None:
    # Imported lazily so multiprocessing worker children never load Qt:
    # freeze_support() below short-circuits them before this runs.
    from flexprobe_gui.app import main
    main()


if __name__ == "__main__":
    # CRITICAL for the frozen .app. The parallel tiling stage uses a
    # ProcessPoolExecutor, and on macOS multiprocessing's "spawn" start method
    # re-launches THIS executable once per worker. Without freeze_support(),
    # each worker process would fall straight through to main() and open a whole
    # new GUI window — running with --workers 8 spawned 8 copies of the app.
    # freeze_support() makes those worker processes run the pickled task and
    # exit instead of launching the UI. It is a no-op in the normal parent
    # process, so a regular launch is unaffected.
    multiprocessing.freeze_support()
    _launch()
