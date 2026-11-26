# 旧アップロードの照合

確認した旧HEAD: `a67b399c95fd8078cab715e54d2970e4227e7d3b`。

現在取得できる履歴は、READMEの初期コミットと、2026-06-30の`Add files via upload`の2件です。一括追加された時点のファイルに、名前と中身の取り違えがあります。その後に内容を改変したコミットは確認されていません。履歴だけからアップロード操作の具体的な原因や実行者の操作手順を断定することはできません。

原本とのSHA-256または改行だけを正規化したハッシュの照合で、例えば次が確認できました。

| 旧アップロード名 | 中に入っていた原本の内容 |
| --- | --- |
| `calcollQ_TAKAKI_newton_q_xy.f90` | `dof_layout.f90` |
| `calweight_TAKAKI_2D_q.f90` | `calcollQ_TAKAKI_newton_q_xy.f90` |
| `COMV1T_TAKAKI_q.f90` | `calcollQ_TAKAKI_newton_q_xyz.f90` |
| `COMV3d_TAKAKI_q_Mxyz.f90` | `calweight_TAKAKI_2D_q.f90` |
| `Straight_Duct_gene_TAKAKI_3D_LMN.f90` | `COMV1T_TAKAKI_q.f90` |
| `test_dof_layout.f90` | 定常NS主プログラム |
| `test_sparse_triplet_builder.f90` | 定常Stokes主プログラム |

全ファイルのハッシュと照合結果は`previous_upload_audit.json`に保存しています。これはGitHubによる自動改変を立証する資料ではなく、現在確認できる登録内容とローカル原本の対応を示す資料です。

今回の登録対象は元ファイル名のままコピーし、登録対象のバイト列と原本のSHA-256の完全一致を確認します。変更前のGit履歴は保持し、履歴を書き換えるforce pushは使用しません。
