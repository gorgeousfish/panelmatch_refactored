# PanelMatch: 面板数据匹配和因果推断分析工具 (修复版)

## 简介

PanelMatch是一个用于面板数据匹配和因果推断分析的Stata命令。该命令可以根据处理历史找到相似的控制单位，并计算处理效应。此版本（v1.0.1）修复了原版中的多项稳定性和兼容性问题。

## 主要修复内容

1. **文件加载机制修复**：
   - 去除了`panelmatch.ado`中的`include`语句，改用Stata的`adopath`机制
   - 添加了对函数定义的清理，避免重复定义错误

2. **硬编码路径问题修复**：
   - 在`panelestimate.ado`中使用`findfile`命令替代硬编码路径
   - 在测试脚本中也使用了相同的改进

3. **处理历史变量生成改进**：
   - 将动态长度的字符串`str\`=\`lag'+1'`替换为固定的足够大的`str32`
   - 避免了使用过大的滞后期可能导致的长度限制问题

4. **Mata函数命名空间修复**：
   - 为所有Mata函数添加了`panelmatch_`前缀，避免与其他包冲突
   - 降低了Mata版本需求，从16.0降到14.0，提高了兼容性

5. **检查与诊断功能增强**：
   - 添加了更详细的参数显示和诊断信息
   - 增加了数据集状态检查

## 安装方法

### 从GitHub安装

```stata
net install panelmatch_refactored, from(https://raw.githubusercontent.com/gorgeousfish/panelmatch_refactored/main) replace
```

### 手动安装

1. 下载此仓库中的所有文件
2. 将所有文件放置于您的个人ado文件夹中（可通过`sysdir`命令查看）
3. 重启Stata

## 使用方法

### 基本用法

```stata
panelmatch treatvar, unitvar(id_var) timevar(time_var) lag(3)
```

### 使用精炼方法

```stata
panelmatch treatvar, unitvar(id_var) timevar(time_var) lag(3) ///
                     refinement(maha) covsformula(x1 x2 x3) sizematch(5)
```

### 计算处理效应

```stata
panelmatch treatvar, unitvar(id_var) timevar(time_var) lag(3) ///
                     refinement(maha) covsformula(x1 x2 x3) ///
                     outcome(y_var) lead(0/3) savematched(matched_data) replace

panelestimate using matched_data.dta
```

## 注意事项

- 此修复版兼容Stata 14及以上版本
- 对于大型数据集，匹配过程可能需要较长时间

## 引用信息

如使用此包进行研究，请引用原始文献和此修复版。

## 问题反馈

如有问题或建议，请通过GitHub issues提交。

## 项目结构
- `src/` - 源代码目录
  - `panelmatch.ado` - 主命令实现
  - `panelmatch_aux.ado` - 辅助函数（修复版1.0.1）
  - `panelmatch_mata.do` - Mata矩阵计算函数（修复版）
  - `panelmatch_ps_model.ado` - 倾向得分模型相关函数
  - `panelmatch_error.ado` - 错误处理模块
  - `panelestimate.ado` - 处理效应估计函数

## 特性
1. 支持基于处理历史的匹配
2. 多种精炼方法：无精炼、马氏距离、倾向得分
3. 丰富的效应估计选项
4. 详细的结果报告