"""Reading the index, for the guards whose subject is the tree itself.

Two guards in this suite ask git what is tracked rather than walking the
working tree: the identifier guard in ``test_sanitized_identifiers``, which
scans every tracked text file for collaborator names, and the dead-name guards
in ``test_qa_seahub_rename``, which take both their subjects and their
references from tracked modules. Each has its own reason to prefer the index,
recorded at its call site.

They had a copy of this plumbing each, and only one of the two scrubbed git's
own variables out of the environment -- so the other was correct by coincidence
of the directory it happened to run in, and a fix to either reached one guard.
The point of one copy is that the failure mode below is silent, so a guard
whose plumbing is wrong does not report anything; it passes.
"""

from __future__ import annotations

from functools import cache
import os
from pathlib import Path
import subprocess

import pytest


def git_free_env() -> dict[str, str]:
    """The environment without git's own variables.

    Git hands its own state down to its hooks, and what it hands down depends
    on the checkout. Measured, from git 2.54: a plain one exports
    ``GIT_INDEX_FILE=.git/index`` and no ``GIT_DIR``; a linked worktree -- which
    is how this repo is usually worked on -- exports both, absolute, naming
    ``.git/worktrees/<name>``. A subprocess that inherits either then answers
    about the repository the variables came from rather than the one it was
    asked about, and does it quietly:

    * An absolute ``GIT_DIR`` makes git skip repository discovery, so ``-C``
      some directory stops selecting anything. The listing describes the
      inherited repository, with names relative to *its* top level, at status 0
      -- and pointed at a throwaway repo built by a test, it answers about the
      real repository entirely.
    * A relative ``GIT_INDEX_FILE`` resolves against the discovered top level
      rather than against git's own directory, so it is harmless below the root
      and fatal in a linked worktree, where ``.git`` is a file and
      ``.git/index`` is therefore not a path. What it does silently is name an
      index that is not there -- a partial commit points it at a temporary one
      -- because git reads a missing index as an empty one: no files, status 0.

    Either way the failure arrives as a guard that checked nothing, or checked
    somewhere else, rather than as an error. The ``pre-commit`` hook that runs
    this suite before every commit is that environment, and it caught both
    halves of this suite's git usage in turn: the listing below, and the
    miniature repo ``test_qa_seahub_rename`` builds to pin it.

    Every ``GIT_*``, not the eight or so that name a location. The wider cut
    also drops ``GIT_CONFIG_COUNT`` and friends, which some container images
    use to inject ``safe.directory`` -- but losing those makes git exit 128 and
    fail the guard out loud, while missing one location variable makes it read
    the wrong index in silence. Between the two, take the one that shouts.
    """
    return {
        name: value for name, value in os.environ.items() if not name.startswith("GIT_")
    }


@cache
def tracked_files(directory: Path) -> tuple[Path, ...]:
    """Every tracked file under ``directory`` that still exists on disk.

    ``git ls-files`` rather than ``rglob``, because a working tree holds plenty
    that is not published: real S3 listings, QA outputs and scratch CSVs, an
    unfinished module, a virtualenv full of third-party code. None of that may
    fail somebody's unrelated commit, and none of it vouches for anything.
    ``git add`` is where a guard built on this takes effect. What each caller
    gains by that is its own, and is recorded where it calls in.

    ``directory / name`` for every name, sorted by name, and no filtering by
    what a file is: the callers want different halves of the listing -- one
    takes it by suffix, the other by which modules it ships -- and a helper that
    guessed which would have to be told anyway. Absolute in, absolute out;
    passing a relative directory works and gets its own cache entry for the same
    tree, which no caller wants.

    Run with ``-C directory`` and no pathspec. ``ls-files`` lists what it finds
    under the directory it runs in and prints it relative to the same place, so
    the join above holds from anywhere inside a checkout, at the top level or
    below it. This reverses an earlier decision -- the repo root with a
    ``bcp`` pathspec, chosen because running at the root was what made a hook's
    inherited variables survivable -- and it is :func:`git_free_env` that makes
    the reversal safe: with those variables gone, the directory is the only
    thing deciding which repository and which subtree are listed. Measured, both
    ways: with an absolute ``GIT_DIR`` inherited, ``-C`` selects nothing, the
    names come back relative to the inherited repository's top level, and every
    join then points at a file that is not there.

    The directory picking the repository is a real difference from the pathspec
    form, not only a safer one: point this at a submodule or a nested checkout
    and it answers out of that repository's index, where the pathspec form
    reported the gitlink and nothing under it. No caller does, and both would
    have to be argued afresh if one did.

    Which is why an empty listing fails here instead of reading as a clean
    tree. No files and status 0 is the shape most of those failures take, and it
    is indistinguishable from success to any caller. A git failure fails the
    test too, rather than skipping it -- a guard that passes without having
    checked anything is worse than one that errors. What this cannot check is
    that anything survived a *caller's* filter, which it knows nothing about; a
    caller narrowing the listing to nothing has to say so itself.

    A listed path that is no longer on disk is dropped, since no caller can read
    one, but a listing where *every* path is missing fails: that is a sparse
    checkout or a wholesale staged deletion, and reporting it as nothing tracked
    would send the reader looking for a listing scoped to the wrong place. Which
    direction a dropped path errs in depends on what the caller does with it, so
    that too is argued at the call sites.
    """
    try:
        result = subprocess.run(
            ["git", "-C", str(directory), "ls-files", "-z"],
            capture_output=True,
            text=True,
            check=False,
            env=git_free_env(),
        )
    except FileNotFoundError:  # pragma: no cover - git absent
        pytest.fail("no git binary, so tracked files cannot be listed")
    if result.returncode != 0:
        pytest.fail(
            f"git ls-files exited {result.returncode} in {directory}, so this "
            f"guard cannot see what is tracked: {result.stderr.strip()}"
        )

    names = sorted(name for name in result.stdout.split("\0") if name)
    if not names:
        pytest.fail(
            f"git listed nothing tracked under {directory}, so this guard would "
            "pass without having read anything"
        )
    paths = tuple(
        path for path in (directory / name for name in names) if path.is_file()
    )
    if not paths:
        pytest.fail(
            f"git listed {len(names)} tracked files under {directory} and none of "
            "them are on disk -- a sparse checkout, or every one staged for deletion"
        )
    return paths
