# ブランチ管理規則

このファイルは、このプロジェクトへのコミットの仕方を説明するものです。

## ブランチ管理ルール　（暫定）

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

## QDDのバージョニング方法

- `main` へのコミットの内容に応じてGit TagのSemVerを手動で更新します
  - `main` へのPRがマージされた後、release noteの作成画面から新しいタグを作成してください。
