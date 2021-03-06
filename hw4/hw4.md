# 第三次作业报告

董莫 027

## 任务

模仿 PowerPoint 写一个曲线设计与编辑工具

- 输入有序点列（型值点），实时生成分段的三次样条曲线
- 可修改拖动型值点的位置（保持整条曲线 $C^2$）
- 可编辑型值点处的切线信息，成为 $G^1$ 或 $G^0$ 

## 目的

- 学习三次样条函数的求解
- 了解曲线设计和编辑工具的原理

## 要求

- 使用 Utopia 框架或其他语言（matlab，python 等）和框架
- deadline：2020 年 11 月 14 日晚

## 思路

本次作业实现的是分段三次样条插值。重点是确保过程中不同级别的光滑（$C0,C1,C2,G0,G1$）。

**分段三次样条插值**：曲线的三次样条插值和函数的分段三次样条插值类似，相邻的每两个插值点之间用一个三次多项式拟合，共有4个参数，n段共4n个参数，每段首尾的函数值可列两个约束方程，曲线中间的插值点处要确保两边的一阶导数和二阶导数连续，每个点处两个约束方程，再加上两个边界条件（固定边界或者自由边界），一共4n个方程，正好可以求解。因为约束具有局部性，所以可以通过适当的变换将矩阵的非零元集中到对角线附近。可采用追赶法求解。

**参数连续性**：参数连续性取决于参数化的方法，$C0$连续只需要插值点两侧没有断开即可。$C1$连续需要插值对点两侧的参数求导相等，在我的作业中，我采用的是最常见的均匀参数化方法，在求解方程组的时候，约束条件即为点的左右两侧一阶导数相等，所以能够保证$C1$连续，$C2$连续也是同理。在程序里可以通过切线长来显示一阶导数的大小，$C1$连续需要保证两边切线长度大小都相等。

**几何连续性**：$G0$连续和$C0$连续要求相同，$G1$连续需要在每一点切线方向相同，所以只需要保证插值点两侧的切线方向相同即可，而不必要保证大小。

## 程序实现

- 本次作业使用了c++语言，配合基于opengl的imgui实现图形操作界面，使用vs 2019 IDE编译运行。

- 使用方法：
  - 鼠标点击绘制顶点，每次绘制一个顶点会自动连接到已有曲线的最末端
  - 勾选“edit curve”进入曲线编辑模式，此时界面上会显示各个插值点，以及插值点处的切线
  - 编辑模式下可以点击移动插值点，或者切线的端点，再次点击停止移动
  - 移动插值点的时候曲线会随之产生变化，保证$C2$连续性
  - 当一个顶点的切线被编辑后，切线方向会被保留，此时整体不具有$C2$连续性
  - 左边（t减小方向）的切线移动会强制使右边的切线与其长度相同方向相反，保证$G1$连续
  - 点击“clear”清除所有顶点

## 运行结果

视频里演示了如下情况：
- 点击添加新的顶点
- 移动调整顶点的位置
- 调整顶点的切线，使其仅具有$G1$连续
- 调整顶点的切线，使其具有尖锐的角
- 清除