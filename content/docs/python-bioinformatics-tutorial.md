---
title: "Python生物信息学入门教程"
date: 2024-12-01
description: "从基础语法到生物数据处理的完整Python教程"
weight: 1
---

# Python生物信息学入门教程

## 第一章：Python基础

### 1.1 环境配置

```bash
# 安装Anaconda
wget https://repo.anaconda.com/archive/Anaconda3-latest-Linux-x86_64.sh
bash Anaconda3-latest-Linux-x86_64.sh

# 创建生物信息学环境
conda create -n bioinfo python=3.9
conda activate bioinfo
```

### 1.2 必要库安装

```bash
# 核心科学计算库
conda install numpy pandas matplotlib seaborn scipy

# 生物信息学专用库
conda install -c bioconda biopython pysam bcftools samtools

# 机器学习库
conda install scikit-learn
```

## 第二章：生物序列处理

### 2.1 使用Biopython处理FASTA文件

```python
from Bio import SeqIO
import pandas as pd

# 读取FASTA文件
sequences = []
for record in SeqIO.parse("genome.fasta", "fasta"):
    sequences.append({
        'id': record.id,
        'length': len(record.seq),
        'gc_content': (record.seq.count('G') + record.seq.count('C')) / len(record.seq)
    })

# 创建DataFrame
df = pd.DataFrame(sequences)
print(df.head())
```

### 2.2 序列质量分析

```python
def analyze_sequence_quality(fastq_file):
    """分析FASTQ文件的序列质量"""
    from Bio import SeqIO
    import numpy as np
    
    qualities = []
    lengths = []
    
    for record in SeqIO.parse(fastq_file, "fastq"):
        qualities.extend(record.letter_annotations["phred_quality"])
        lengths.append(len(record.seq))
    
    return {
        'mean_quality': np.mean(qualities),
        'mean_length': np.mean(lengths),
        'total_reads': len(lengths)
    }
```

## 第三章：基因组变异分析

### 3.1 使用pysam处理BAM文件

```python
import pysam

# 打开BAM文件
bamfile = pysam.AlignmentFile("sample.bam", "rb")

# 统计比对情况
mapped = 0
unmapped = 0

for read in bamfile.fetch():
    if read.is_unmapped:
        unmapped += 1
    else:
        mapped += 1

print(f"Mapped reads: {mapped}")
print(f"Unmapped reads: {unmapped}")
print(f"Mapping rate: {mapped/(mapped+unmapped)*100:.2f}%")
```

### 3.2 变异频率统计

```python
def count_variants(vcf_file):
    """统计VCF文件中的变异类型"""
    import pysam
    
    variants = {'SNP': 0, 'InDel': 0}
    
    vcf = pysam.VariantFile(vcf_file)
    for record in vcf.fetch():
        if len(record.ref) == 1 and len(record.alts[0]) == 1:
            variants['SNP'] += 1
        else:
            variants['InDel'] += 1
    
    return variants
```

## 第四章：转录组数据分析

### 4.1 表达矩阵处理

```python
import pandas as pd
import numpy as np
from scipy import stats

# 读取表达矩阵
expression = pd.read_csv("expression_matrix.csv", index_col=0)

# 数据标准化
expression_norm = np.log2(expression + 1)

# 差异表达分析
def differential_expression(control_samples, treatment_samples):
    """简单的差异表达分析"""
    results = []
    
    for gene in expression_norm.index:
        control_expr = expression_norm.loc[gene, control_samples]
        treatment_expr = expression_norm.loc[gene, treatment_samples]
        
        # t检验
        t_stat, p_value = stats.ttest_ind(control_expr, treatment_expr)
        
        # 计算fold change
        fc = np.mean(treatment_expr) - np.mean(control_expr)
        
        results.append({
            'gene': gene,
            'log2_fc': fc,
            'p_value': p_value,
            't_statistic': t_stat
        })
    
    return pd.DataFrame(results)
```

## 第五章：数据可视化

### 5.1 基因组数据可视化

```python
import matplotlib.pyplot as plt
import seaborn as sns

# 设置图形样式
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

# 绘制质量分布图
def plot_quality_distribution(qualities):
    """绘制测序质量分布图"""
    plt.figure(figsize=(10, 6))
    plt.hist(qualities, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('Quality Score')
    plt.ylabel('Frequency')
    plt.title('Sequencing Quality Distribution')
    plt.grid(True, alpha=0.3)
    plt.show()

# 绘制火山图
def plot_volcano(de_results, fc_threshold=1, p_threshold=0.05):
    """绘制差异表达火山图"""
    plt.figure(figsize=(10, 8))
    
    # 筛选显著差异基因
    significant = (np.abs(de_results['log2_fc']) > fc_threshold) & \
                 (de_results['p_value'] < p_threshold)
    
    # 绘制散点图
    plt.scatter(de_results['log2_fc'], -np.log10(de_results['p_value']), 
               c=significant, cmap='RdYlBu', alpha=0.6)
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 P-value')
    plt.title('Volcano Plot')
    plt.axhline(y=-np.log10(p_threshold), color='red', linestyle='--')
    plt.axvline(x=fc_threshold, color='red', linestyle='--')
    plt.axvline(x=-fc_threshold, color='red', linestyle='--')
    plt.show()
```

## 练习题

1. 编写函数计算DNA序列的GC含量
2. 从FASTQ文件中提取质量得分最高的100条序列
3. 比较两个样本的基因表达相关性
4. 绘制基因在染色体上的分布图

## 参考资料

- [Biopython官方文档](https://biopython.org/)
- [Pandas用户指南](https://pandas.pydata.org/)
- [生物信息学Python实战](https://github.com/bioinformatics-python)
