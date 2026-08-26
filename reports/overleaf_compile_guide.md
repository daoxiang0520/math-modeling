# Overleaf 导入与编译指南（C 题论文）

> 适用对象：`outputs/overleaf_cumcm_*.zip`（由框架自动打包）或按相同结构手动整理的目录。
> 目标：在 Overleaf 上编译 2025 高教社杯 C 题《NIPT 的时点选择与胎儿的异常判定》国赛 LaTeX 论文。

## 1. 快速开始（三步）

1. **新建项目并导入 zip**：Overleaf 右上角 `New Project → Upload Project`，选择
   `outputs/overleaf_cumcm_*.zip` 上传。项目根目录结构如下：

   ```text
   项目根/
   ├── figures/          # 论文引用的全部图（PNG/PDF/SVG）
   └── paper/
       ├── main.tex      # 主文档（封面/摘要/目录 + 章节 input）
       ├── references.tex
       └── sections/
           ├── 1_restatement.tex
           ├── 2_analysis.tex
           ├── 3_assumptions.tex
           ├── 4_symbols.tex
           ├── 5_problem1.tex
           ├── 6_problem2.tex
           ├── 7_problem3.tex
           ├── 8_problem4.tex
           ├── 8_sensitivity.tex
           ├── 9_evaluation.tex
           └── A_code.tex
   ```

2. **设置主文档**：`Menu → Main document` 选择 `paper/main.tex`。
   图片路径 `../figures/...` 按 `paper/` 的相对位置解析到项目根 `figures/`，
   因此必须保持该两级目录结构。

3. **选择 XeLaTeX 并编译**：`Menu → Compiler` 选择 **XeLaTeX**（`main.tex` 首行已有
   `% !TEX program = xelatex`，Overleaf 通常自动识别）。点击 `Recompile`，
   **至少编译两遍**（第一遍生成目录/交叉引用，第二遍解析 `\ref`、`\cite`、目录页码）；
   或开启 `Menu → Auto Compile` 后手动多 Recompile 一次。

## 2. 字体设置（关键）

模板默认 `\documentclass[fontset=mac, ...]{ctexart}`，依赖 macOS 字体
（Songti SC / Heiti SC / Kaiti SC）。**Overleaf 运行在 Linux（TeX Live）上，没有这些字体**，
直接编译会出现字体缺失或回退混乱。

修改 `paper/main.tex` 第 24 行：

```latex
% 原：\documentclass[fontset=mac, 12pt, a4paper]{ctexart}
\documentclass[fontset=fandol, 12pt, a4paper]{ctexart}
```

`fontset=fandol` 使用 TeX Live 自带的 Fandol 中文字体（宋/黑/楷/仿宋），Overleaf 开箱即用。
若希望用 Windows 字体可改用 `fontset=windows`（Overleaf 需自行上传字体，一般不需要）。

## 3. 宏包依赖

论文只依赖 TeX Live 标准宏包，Overleaf 默认已全部安装：
`ctex`（ctexart）、`geometry`、`amsmath`、`graphicx`、`fontspec`、`titlesec`、
`tocloft`、`booktabs`、`array`、`xcolor`、`listings`。无需上传任何 `.cls` / `.sty` / `.bib` 文件。

## 4. 常见问题排查

- **编译报错找不到字体**（`fontspec` font not found）：未改 `fontset=fandol`，见第 2 节。
- **图片不显示或 `File ... not found`**：`figures/` 未上传，或主文档未设为 `paper/main.tex`，
  或图片路径被改动。
- **目录页码是 `??`**：只编译了一遍；多 Recompile 一次。
- **引用显示 `图??` / `表??`**：`\ref` 无对应 `\label`；检查 `sections/*.tex` 的 `\label{fig:...}`。
- **中文乱码或字体变形**：确保编译器为 XeLaTeX（pdfLaTeX 无法处理 ctex 字体方案）。
- **编译很慢**：论文含 40+ 张图；Overleaf 免费版编译队列较长，属正常。
- **页数/排版与本地不同**：Overleaf 用 Fandol 字体，字宽与本地（mac/Windows）略有差异，属正常跨平台偏差。

## 5. 本地编译对照（可选）

本机安装 TeX Live / MacTeX 后：

```bash
cd paper
xelatex -interaction=nonstopmode main.tex
xelatex -interaction=nonstopmode main.tex   # 第二遍生成目录
```

结果 PDF 输出到 `paper/main.pdf`。Overleaf 与本地结果应一致（仅字体细节略有差异）。

## 6. 论文结构速览（对应国赛要求）

- `main.tex`：封面（标题/摘要/关键词）、目录、章节装载
- `1_restatement`：问题重述；`2_analysis`：问题分析；`3_assumptions`：模型假设
- `4_symbols`：符号说明；`5_problem1`–`8_problem4`：四题模型建立与求解（各以「问题小结」收尾）
- `8_sensitivity`：灵敏度分析；`9_evaluation`：模型评价与推广；`A_code`：附录代码
