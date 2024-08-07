<!--
 * @Author: cy 1348987670@qq.com
 * @Date: 2024-07-26 11:32:09
 * @LastEditors: cy 1348987670@qq.com
 * @LastEditTime: 2024-07-26 11:51:01
 * @FilePath: \ZJU\README.md
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
-->

## Overview

This code can be used to efficiently &eta;-threshold maintainance in dynamic uncertain graphs.

## Function

1. **Initial_threshold_compute_map()**: the &eta;-threshold decomposition algorithm for uncertain graphs.
2. **insert_threshold_compare()**: integrates all the maintenance algorithms for edge insertion.

- insert_threshold_compute_candidate(): our proposed basic &eta;-threshold maintenance algorithm.
- insert_threshold_compute_restriction_point(): the basic &eta;-threshold maintenance algorithm with contraction optimization.
- insert_threshold_compute_update_range(): the basic &eta;-threshold maintenance algorithm with candidate pruning optimization.
- insert_threshold_compute_batchUP(): the basic &eta;-threshold maintenance algorithm with batch update optimization.
- insert_threshold_compute_opt(): the basic &eta;-threshold maintenance algorithm with all three optimizations.

3. **delete_threshold_compare()**: integrates all the maintenance algorithms for edge deletion.

- The related algorithms are the same as above.

4. **increase_threshold_compare()** and **decrease_threshold_compare()**: integrates all the algorithms for edge probability change.
