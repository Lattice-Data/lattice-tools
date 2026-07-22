---
name: github-pr-review
description: "Review a GitHub pull request and write a short, readable PR description for it. Produces two things: (1) a prioritized code review covering correctness and logic, security, cleanliness, style, tests, and docs, with blocking issues first and nits collapsed; and (2) a tight 2-4 sentence PR summary plus a bullet list of changes, short enough that people actually read it. Pulls the PR via the gh CLI from a PR number or URL, or works from a diff the user pastes. Use whenever the user wants to review a PR or pull request, review or critique a diff or code changes, check a PR before merging, get feedback on their code, or write, summarize, or improve a PR description. Trigger on mentions of PR, pull request, \"review this diff\", \"gh pr\", a GitHub PR link, or \"summarize these changes\". Never posts to GitHub automatically; it presents the outputs and gives the commands to post if asked."
---


# GitHub PR Review

Review a pull request and produce two artifacts:

1. **A code review** - prioritized findings on correctness/logic, security, cleanliness, style, tests, and docs.
2. **A PR description** - 2-4 sentences of what-and-why plus a short bullet list of changes, kept short enough to actually get read.

Produce both by default. If the user clearly wants only one (e.g. "just write me a description"), produce only that one.

Guiding principle: be useful to a reviewer who has limited time. Say what matters, in priority order, and stop. No praise filler, no restating the obvious.

## Step 1 - Get the PR

**Primary path: the `gh` CLI.** Assume it is installed and authenticated. Given a PR number or URL, run:

```bash
gh pr view <pr> --json number,title,body,author,baseRefName,headRefName,url,additions,deletions,changedFiles,files,commits
gh pr diff <pr>
```

`gh pr view` gives metadata (title, existing body, branch names, file list, commit messages). `gh pr diff` gives the unified diff to review. `<pr>` can be a number (`482`), a URL, or `owner/repo#number`.

If `gh` is missing or not authenticated, say so in one line, then fall back.

**Fallback: a pasted diff.** If the user pastes a diff (or the output of `git diff <base>...<head>`), review that directly. If the user gives only a URL and `gh` is unavailable, ask them to paste the diff or run `gh pr diff <pr>` themselves - do not try to fetch a PR over the web, since private PRs need auth and web fetches are unreliable here.

## Step 2 - Read for review

Read the whole diff before writing anything. Then judge, in this order of importance:

1. **Correctness and logic** - bugs, wrong conditionals, off-by-one, unhandled null/empty/error cases, incorrect assumptions, mismatch between the code and what the PR says it does, concurrency/ordering hazards, resource leaks. This is the top priority.
2. **Security and data safety** - injected/unsanitized input, secrets or credentials committed, unsafe deserialization, missing authz checks, destructive operations (deletes, overwrites, migrations) without guards. Weigh this heavily for anything touching storage, credentials, or user data.
3. **Cleanliness and maintainability** - dead code, duplication, unclear names, functions doing too much, needless complexity, leftover debug output or commented-out code.
4. **Tests** - is new logic covered? Do the tests actually assert the behavior, or just run it? Are edge cases tested?
5. **Docs and comments** - missing or misleading docstrings/comments where the code is non-obvious. Do not ask for comments on self-explanatory code.
6. **Style and format** - only deviations from the repo's own conventions. Anything a formatter or linter auto-fixes goes in Nits, collapsed - do not spend real estate on it.

Respect existing conventions over personal preference. If the repo has a linter config, style guide, or CONTRIBUTING file and it is visible, defer to it. Do not impose a style the project does not use.

When the diff alone is not enough to judge correctness (e.g. a changed function calls something not shown), read the full file from the local checkout if the repo is available, or `git show <headRef>:<path>`. If you still cannot tell, say so and mark it as a Question rather than guessing.

When investigating files beyond the diff, limit deep reading to the files most relevant to correctness and security. Skip full investigation of the rest, and note in a Question that they were not reviewed in full rather than reading everything.

For large diffs: prioritize files with real logic changes. Skim generated code, lockfiles, and vendored dependencies, and say you skimmed them rather than reviewing line by line.

## Step 3 - Write the code review

Use severity tags so the reader can triage. Order sections by severity. **Omit any section that is empty** - do not print empty headers.

- `[BLOCKING]` - must fix before merge: a bug, security issue, data-loss risk, or logic that does not do what it should.
- `[SHOULD-FIX]` - a real issue but not a merge blocker: missing test for new logic, a fragile edge case, a maintainability problem.
- `[NIT]` - style or preference. Collapse these; one short line each, grouped.
- `[Q]` - intent is unclear and you are inferring or guessing. State the assumption and ask.

For each finding: point to `path:line`, say what is wrong, why it matters in one clause, and the fix if it is not obvious. Be specific; a finding the author cannot act on is noise.

Distinguish confidence in the wording: a definite bug ("this returns None when the list is empty, which the caller dereferences") reads differently from a possible one ("if `items` can be empty here, this dereferences None - can it?"). Do not state guesses as facts.

Open with a one-line verdict so the reader knows the outcome before the details.

**Review template:**

```markdown
## Review: <title> (#<number>)

**Verdict:** <one line, e.g. "One blocking bug in the retry path; rest is sound." or "No blocking issues; two small suggestions.">

### Blocking
- **`src/foo.py:42`** - <what is wrong>. <why it matters>. <fix if non-obvious>.

### Should fix
- **`src/bar.py:88`** - <issue>.

### Nits
- `src/baz.py:12` <one-liner>. `src/baz.py:30` <one-liner>.

### Questions
- **`src/qux.py:15`** - <assumption>; is that intended?
```

If nothing is wrong, say so plainly in the verdict and keep the body short. Do not invent issues to fill space.

## Step 4 - Write the PR description

Short and readable is the whole point. Lead with what changed and why (2-4 sentences), then a tight bullet list of the actual changes. Cover the "why" in plain terms - the problem this solves or the reason for the change - because the diff already shows the "what" in detail.

Do not pad. No "This PR..." throat-clearing if it can be cut. Skip a bullet list entirely for a one-line change; just write the sentence.

**Description template:**

```markdown
## PR description

<2-4 sentences: what changed and why it was needed.>

**Changes:**
- <concrete change>
- <concrete change>
- <concrete change>
```

## Step 5 - Post the code review as a comment

Post findings; never decide approval. This skill reports blocking issues, should-fix items, nits, and questions - the merge decision (approve, request changes, or nothing) is a human call, always.

**Post the code review as a PR comment.** Write the Step 3 review to `review.md`, then:

```bash
gh pr comment <pr> --body-file review.md
```

Post the code review only. Do not post the PR description here, do not edit the PR body, and never run `gh pr review --approve` or `gh pr review --request-changes` from this skill - approving or blocking a PR is not this skill's decision to make, regardless of how many blocking findings there are.

**Interactive runs (a person is driving, not CI).** Do not post automatically. Present both artifacts inline in fenced code blocks and offer the commands below; run them only if the user explicitly says to.

```bash
# Post the review as a comment:
gh pr comment <pr> --body-file review.md

# Set the PR description (this REPLACES the existing body):
gh pr edit <pr> --body-file pr-description.md
```

For multi-line bodies use `--body-file` (write the text to a file first) rather than `--body`. Note that `gh pr edit --body` overwrites the description; to keep the existing body, append to it instead of replacing.

## Output style

Match the tone of a direct, senior reviewer: concise, specific, no flattery. Facts, inferences, and guesses stay visibly distinct (the severity tags and hedged wording do this work). Prefer fewer high-signal findings over an exhaustive list. When in doubt about whether something is worth mentioning, it probably is not.
