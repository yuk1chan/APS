# 環境
- Juliaが使える環境。
- Juliaのversionは1.0.1。
- Plotsパッケージを用いている。

# 使い方
コマンドライン上で`julia tdse.jl`を実行すれば、計算結果がgifとして出力される。

また、jupyter notebookでも使える。

# 結果
## aについて
- 波束はポテンシャルで反射されるものもあれば、ポテンシャルを透過するものもある。
- 全ての時刻において波束はガウス型という訳ではない。
  ポテンシャルで反射や透過をする部分があるため。
- 入射波がポテンシャルに衝突する時刻を$`t_{inc}`$
  反射波が初期位置に戻るまでの時刻を$`t_{ref}`$とする。
  $`t_{inc}`$ = $`t_{ref}`$は成立しない。
  それは、透過波が存在しているためだと考えられる。
  透過波の