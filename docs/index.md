# SpatialPerturb 文档

SpatialPerturb 是一个面向空间转录组与 Perturb-seq 联合分析的 Python 工具包。

当前版本提供：

- 统一的小写包入口 `spatialperturb`
- 基础 CLI，用于检查已安装版本
- 将基因集字典转换为二值 signature matrix 的工具函数

安装：

```bash
pip install SpatialPerturb
```

快速开始：

```python
import spatialperturb as sp
from spatialperturb import build_signature_matrix

print(sp.__version__)
print(build_signature_matrix({"program": ["STAT1", "IRF1"]}))
```

更多接口见 [API 参考](api.md)。
