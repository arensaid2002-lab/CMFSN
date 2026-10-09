# Rules for Claude in CMFSN

## Branch
- Before starting any task, ask me which branch to work on. Don't assume it's the
  branch currently checked out.
- Check that branch is checked out (`git branch --show-current`) before editing.
  Don't switch branches yourself without my OK.
- Never work on `main`. If `main` is checked out, stop and tell me.

## Safety
- Never delete, move or rename files without asking me first.
- Never overwrite an existing file I wrote without showing me the change first.
- Only work inside this repository folder.
- Before starting, check `git status`. If there are uncommitted changes I didn't
  mention, ask before touching them.

## Git
- Never commit, push, merge, rebase or reset unless I ask. I commit myself in
  GitHub Desktop.
- Never use destructive commands (`git reset --hard`, `git clean`,
  `git push --force`, `git checkout -- <file>`, `git restore`).

## Code
- C++20 with CMake. The build must have zero warnings and all tests must pass
  before you say a task is done.
- Don't install software or add dependencies without asking.
- Keep changes small: one feature per task, so I can review and commit each step.

## Learning
- Don't write the code for me: explain, give hints and review. I type the code and
  run the builds myself. Only write code if I explicitly ask for it.
- I'm learning C++ and git. Explain what you changed and why, in simple terms.
- Comment the code, especially the math (NACA equations, units, coordinate
  conventions).
