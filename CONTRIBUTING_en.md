# Contributing Guidelines

This file explains about how to contribute to our project!

## Branch Management Rules

- We use three types of branches: `main`, `develop`, and `topic`.
  - `topic` branches are for specific features or bug fixes.
- Naming convention for `topic` branches (tentative):
  - `feature/add-circuit`
  - `bugfix/py-version`
- Work should be done on `topic` branches branched off from `develop`.
- Pull requests from `topic` to `develop` should be made.
- When ready to release (e.g., bumping QDD version), pull requests from `develop` to `main` should be made.
  - Upon merging to `main`, the following actions will be triggered by GitHub Actions:
    - Building the wheel
    - Publishing to PyPI

## Q&A

### Can GitHub's Free Plan Handle This?

- GitHub Actions are free for public repositories.
  - [About billing for GitHub Actions](https://docs.github.com/en/billing/managing-billing-for-github-actions/about-billing-for-github-actions)

> **GitHub Actions is free for public repositories**
>
> We take pride in our Open Source legacy, and are happy to provide free CI/CD for public repositories. Check out the doc to see which runners are included.
>
- Since QDD is a public repository, it's covered under this.

### QDD Versioning

- Option 1) Manually Update SemVer Tags Based on Commits to `main`
  - [About SemVer](https://qiita.com/usamik26/items/c8911219b610101e69a9)
  - After a PR to `main` is merged, create a new tag from the release notes creation page.
    - [git tag と GitHub の Release 機能でプロっぽさを出してみよう](https://qiita.com/tommy_aka_jps/items/5b39e4b27364c759aa53#%E3%81%A1%E3%81%AA%E3%81%BF%E3%81%AB)

- Option 2) Automatically Update Git Tags with GitHub Actions When Making Commits to `main`
  - [Reference Article](https://dev.classmethod.jp/articles/howto-use-github-tag/)
  - The prefix rule for commit messages must be followed.
  - Verified by the following repository.
    - [Validation Repository](https://github.com/NT-marlowe/branch-management-test/tree/main)
  - GitHub's release creation Action is unnecessary as QDD is published to PyPI.
