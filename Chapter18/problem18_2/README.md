# 使い方

`./eigen.sh`

を実行することで、input_data.txtにある値を初期値として波動関数を計算する。
計算してgnuplotで波動関数をphi.pngに出力している。

input_data.txtの内容は、

| 要素            |  値  |
|:---------------|-----:|
| $`V_a`$        | 150  |
| a              | 1.0  |
| $`V_b`$        | 10   |
| b              | 0.1  |
| $`\Delta x`$   | 0.01 |
| xmax           | 1.5  |
| parity         |  1   |
| $`E`$          | 1.0  |

となっている。


# 結果


$`V_a = 150, a = 1.0`$

| $`V_b`$| $`b`$ | 推定固有値 $`E`$ |
|:-------|------|---------:|
| 10     |  0.1 |  2.603   |
| 20     |  0.1 |  3.377   |

くらいで見てみると良いと思う。