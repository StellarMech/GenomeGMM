# GenomeGMM

## 1. 简介
GenomeGMM是一个用于预估多倍体基因组大小的R包。它支持读取两列格式的基因组文件，结合多倍体倍性，自动预估基因组大小并生成可视化PDF。

## 2. 安装
- 详情请见[installation.md](https://github.com/StellarMech/GenomeGMM/blob/main/installation.md)

## 3. 使用方法

### 3.1 函数说明

主函数为 `GenomeGMM`，参数如下：

- `data_file`：输入的k-mer频率直方图文件（txt格式），每行两列，第一列为k-mer频率，第二列为该频率的k-mer数。
- `n`：高斯分量数（整数），建议与理论主峰数一致。

### 3.2 示例

```r
library(GenomeGMM)
GenomeGMM("4merged.histo.txt", n = 4)
```


## 4. 结果说明
- 控制台输出包括：估算的基因组大小
- 生成的图像文件展示了k-mer分布直方图、GMM拟合曲线及各分量。
- 图中右下角显示基因组大小预估值（bp、Mb、Gb）。

---

如需更多帮助，请参考包内文档或联系作者。