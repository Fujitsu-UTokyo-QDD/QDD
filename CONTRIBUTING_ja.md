# ブランチ管理規則

このファイルは、このプロジェクトへのコミットの仕方を説明するものです。

## ブランチ管理ルール

- 3つの種類のブランチを使用します: `main`, `develop`, `topic`.
  - `topic` ブランチは、特定の機能やバグ修正用です。
- `topic` ブランチの命名規則（暫定的）:
  - `feature/add-circuit`
  - `bugfix/py-version`
- 作業は、`develop` ブランチから分岐した `topic` ブランチで行います。
- `topic` から `develop` へのプルリクエストを作成してください。
- リリースの準備が整った場合（例: QDDのバージョンを上げる場合）、`develop` から `main` へのプルリクエストを作成してください。
  - `main` にマージされると、GitHub Actionsによって以下のアクションがトリガーされます:
    - wheelのビルド
    - PyPIへの公開

## Q&A

### GitHubの無料プランで対応可能か？

- GitHub Actionsはpublic repositoryでは無料です。
  - [GitHub Actionsの公式説明](https://docs.github.com/en/billing/managing-billing-for-github-actions/about-billing-for-github-actions)

> **GitHub Actions is free for public repositories**
>
> We take pride in our Open Source legacy, and are happy to provide free CI/CD for public repositories. Check out the doc to see which runners are included.
>
- QDDはpublic repositoryなので、これに該当します。

### QDDのバージョニング方法

- オプション1) `main` へのコミットの内容に応じてSemVerタグを手動で更新します
  - [SemVerについて](https://qiita.com/usamik26/items/c8911219b610101e69a9)
  - `main` へのPRがマージされた後、release noteの作成画面から新しいタグを作成してください。
    - [git tag と GitHub の Release 機能でプロっぽさを出してみよう](https://qiita.com/tommy_aka_jps/items/5b39e4b27364c759aa53#%E3%81%A1%E3%81%AA%E3%81%BF%E3%81%AB)

- オプション2) `main` へのコミット時に、GitHub Actionsにより自動的にGit Tagを更新します
  - [参考記事](https://dev.classmethod.jp/articles/howto-use-github-tag/)
  - コミットメッセージのprefixのルールに従う必要があります。
  - 以下のリポジトリで検証済み。
    - [検証リポジトリ](https://github.com/NT-marlowe/branch-management-test/tree/main)
  - QDDはPyPIに公開されるため、GitHubのreleaseを作成するActionは不要です。
