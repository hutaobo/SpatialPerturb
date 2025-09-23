import os, sys
# 关键：让 Sphinx 能 import 到 src 布局下的包
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../src"))

project = "SpatialPerturb"
author = "Taobo Hu"
extensions = [
    "myst_parser",               # Markdown (MyST)
    "sphinx.ext.autodoc",        # 自动 API
    "sphinx.ext.napoleon",       # 支持 NumPy/Google 风格 docstring
    "sphinx_autodoc_typehints",  # 渲染类型注解
    "sphinx.ext.viewcode",       # 显示源码链接
]
html_theme = "furo"
# 若文档源文件既有 .md 也有 .rst，可显式声明（可选）
# source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
