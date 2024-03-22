# Contributing Guidelines

This file explains how to contribute to our project.

## Branch Management Rules (tentative)

- We use three types of branches: `main`, `develop`, and `topic`.
  - `topic` branches are for specific features or bug fixes.
- Naming convention for `topic` branches:
  - Branch for adding a feature: `feature/xxx` (e.g. `feature/add-circuit`)
  - Branch for fixing a bug: `bugfix/xxx` (e.g. `bugfix/py-version`)
- Work is done on `topic` branches branched off from `develop`.
- Pull requests from `topic` to `develop` is made.
- When ready to release (e.g., bumping QDD version), pull requests from `develop` to `main` is made.

## QDD Versioning

- Manually Update Git Tag SemVer Based on Commits to `main`.
  - After a PR to `main` is merged, create a new tag from the release notes creation page.
