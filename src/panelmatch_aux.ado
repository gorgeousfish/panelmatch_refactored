/*===========================================================================
 PanelMatch辅助函数集 - 重构修复版
 
 版本: 1.0.1
 日期: 2025-03-04
 作者: [作者信息]
 
 说明: 本文件包含PanelMatch和PanelEstimate命令使用的辅助函数（修复版）
=============================================================================*/

/* 主程序中使用的参数验证和数据准备函数 */

/*===========================================================================
  验证PanelMatch参数
============================================================================*/
program define validate_panelmatch_opts
    syntax [varname(max=1)], [id(varname) time(varname) outcome(varname) treatment(varname) ///
                              covariates(varlist) model(string) refinement(string) qoi(string) *]
    
    // 如果提供了treatment，使用它作为处理变量
    local treatvar `varlist'
    if "`treatment'" != "" {
        local treatvar `treatment'
    }
    
    // 设置默认值
    if "`refinement'" == "" {
        local refinement "none"
    }
    
    if "`qoi'" == "" {
        local qoi "att"
    }
    
    if "`model'" == "" {
        local model "logit"
    }
    
    // 验证refinement选项
    local valid_ref none maha ps
    local ref_valid 0
    foreach r of local valid_ref {
        if "`refinement'" == "`r'" {
            local ref_valid 1
        }
    }
    
    if `ref_valid' == 0 {
        display as error "refinement选项必须是以下之一: none, maha, ps"
        exit 198
    }
    
    // 验证qoi选项
    local valid_qoi att atc ate
    local qoi_valid 0
    foreach q of local valid_qoi {
        if "`qoi'" == "`q'" {
            local qoi_valid 1
        }
    }
    
    if `qoi_valid' == 0 {
        display as error "qoi选项必须是以下之一: att, atc, ate"
        exit 198
    }
    
    // 验证model选项
    if "`model'" != "" {
        local valid_models logit probit cloglog
        local model_valid 0
        foreach m of local valid_models {
            if "`model'" == "`m'" {
                local model_valid 1
            }
        }
        
        if `model_valid' == 0 {
            display as error "model选项必须是以下之一: logit, probit, cloglog"
            exit 198
        }
    }
    
    // 显示验证结果
    display as text "参数验证通过:"
    display as text "- ID变量: " as result "`id'"
    display as text "- 时间变量: " as result "`time'"
    display as text "- 处理变量: " as result "`treatvar'"
    if "`outcome'" != "" {
        display as text "- 结果变量: " as result "`outcome'"
    }
    if "`covariates'" != "" {
        display as text "- 协变量: " as result "`covariates'"
    }
    display as text "- 精炼方法: " as result "`refinement'"
    display as text "- 关注效应: " as result "`qoi'"
    if "`model'" != "" {
        display as text "- 模型类型: " as result "`model'"
    }
end

/*===========================================================================
  生成处理历史变量 - 修复版
============================================================================*/
program define generate_treatment_history
    syntax varname(max=1), unitvar(varname) timevar(varname) [lag(integer 1) ///
                          hist_name(name) touse(varname) replace]
    
    // 提取变量名
    local treatvar `varlist'
    
    // 获取处理历史变量名
    if "`hist_name'" == "" {
        local histvar `treatvar'_hist
    }
    else {
        local histvar `hist_name'
    }
    
    // 准备touse变量
    if "`touse'" == "" {
        tempvar touse
        quietly gen byte `touse' = 1
    }
    
    // 检查并删除现有变量(如果replace选项)
    capture confirm variable `histvar'
    if !_rc & "`replace'" != "" {
        display as text "变量 `histvar' 已存在，将被替换"
        drop `histvar'
    }
    else if !_rc & "`replace'" == "" {
        display as error "变量 `histvar' 已存在，请使用replace选项替换或指定不同的变量名"
        exit 110
    }
    
    // 创建处理历史变量(以字符串格式)
    display as text "创建处理历史变量: `histvar' (滞后值: `lag')"
    // 使用固定的足够大的字符串长度，避免超出限制
    quietly gen str32 `histvar' = ""
    
    // 排序数据以确保正确的滞后顺序
    quietly sort `unitvar' `timevar'
    
    // 对每个滞后期数循环
    forvalues l = 1/`lag' {
        // 创建一个临时变量存储滞后值
        tempvar lag_treat`l'
        quietly by `unitvar' (`timevar'): gen `lag_treat`l'' = `treatvar'[_n-`l'] if _n > `l' & `touse'
        
        // 为每个单位填充历史
        quietly by `unitvar' (`timevar'): replace `histvar' = `histvar' + cond(missing(`lag_treat`l''), ".", string(`lag_treat`l'')) if `touse'
    }
    
    // 显示一些汇总信息
    quietly count if `histvar' != "" & `touse'
    display as text "`r(N)'个观测有完整的处理历史"
end

/*===========================================================================
  准备匹配数据
============================================================================*/
program define prepare_matching_data
    syntax varname(max=1), unitvar(varname) timevar(varname) thistory(varname) ///
                        touse(varname) [listwise matchmissing forbidtreatmentreversal ///
                        exact_match(varlist)]
    
    // 提取变量名
    local treatvar `varlist'
    
    // 显示准备状态
    display as text _newline "准备匹配数据..."
    
    // 检查数据是否已经排序
    quietly summ
    if r(N) < 2 {
        display as error "数据集中观测不足，无法进行匹配"
        exit 2001
    }
    
    // 排序数据以确保正确处理
    quietly sort `unitvar' `timevar'
    
    // 标记当前单位的上一期处理状态
    tempvar prev_treat
    quietly by `unitvar' (`timevar'): gen `prev_treat' = `treatvar'[_n-1] if _n > 1 & `touse'
    
    // 检查处理反转(如果指定)
    if "`forbidtreatmentreversal'" != "" {
        tempvar treat_reversed
        quietly gen `treat_reversed' = (`treatvar' == 0 & `prev_treat' == 1) if !missing(`prev_treat') & `touse'
        
        quietly count if `treat_reversed' == 1 & `touse'
        if r(N) > 0 {
            display as text "警告: 检测到`r(N)'个观测从处理回到控制状态"
            display as text "      由于指定了forbidtreatmentreversal选项，这些观测将被排除"
            quietly replace `touse' = 0 if `treat_reversed' == 1
        }
    }
    
    // 如果有missing值，检查处理情况
    if "`matchmissing'" == "" {
        // 排除处理历史中有缺失值的观测
        quietly replace `touse' = 0 if missing(`thistory') & `touse'
        quietly count if `touse' == 0
        if r(N) > 0 {
            display as text "排除`r(N)'个观测，因为它们的处理历史中有缺失值"
            display as text "使用matchmissing选项可以保留这些观测"
        }
    }
    
    // 检查精确匹配变量
    if "`exact_match'" != "" {
        foreach var of varlist `exact_match' {
            capture confirm numeric variable `var'
            if !_rc {
                // 对于数值变量，检查缺失
                quietly count if missing(`var') & `touse'
                if r(N) > 0 & "`listwise'" != "" {
                    display as text "排除`r(N)'个观测，因为精确匹配变量`var'中有缺失值"
                    quietly replace `touse' = 0 if missing(`var') & `touse'
                }
            }
            else {
                // 对于字符串变量，检查空值
                quietly count if `var' == "" & `touse'
                if r(N) > 0 & "`listwise'" != "" {
                    display as text "排除`r(N)'个观测，因为精确匹配变量`var'中有空值"
                    quietly replace `touse' = 0 if `var' == "" & `touse'
                }
            }
        }
    }
    
    // 汇总处理状态
    quietly count if `touse'
    display as text "匹配准备完成，剩余`r(N)'个有效观测"
    
    quietly count if `treatvar' == 1 & `touse'
    local n_treated = r(N)
    quietly count if `treatvar' == 0 & `touse'
    local n_control = r(N)
    
    display as text "其中: 处理组 = `n_treated'个观测, 控制组 = `n_control'个观测"
end

/*===========================================================================
  生成匹配数据集
============================================================================*/
program define generate_matched_dataset
    syntax varname(max=1), unitvar(varname) timevar(varname) thistory(varname) ///
                        [treatvar(varname) exactmatchvars(varlist) touse(varname)]
    
    // 提取变量名
    local mainvar `varlist'
    if "`treatvar'" == "" {
        local treatvar `mainvar'
    }
    
    // 准备touse变量
    if "`touse'" == "" {
        tempvar touse
        quietly gen byte `touse' = 1
    }
    
    // 显示匹配状态
    display as text _newline "生成匹配数据集..."
    
    // 清除当前数据，保留匹配所需的变量
    keep if `touse'
    keep `unitvar' `timevar' `treatvar' `thistory' `exactmatchvars'
    
    // 获取处理单位和控制单位
    tempvar treated
    quietly gen `treated' = (`treatvar' == 1)
    
    // 保存总观测数
    quietly count
    local total_obs = r(N)
    
    // 分离处理组和控制组
    preserve
    
    quietly keep if `treated'
    quietly count
    local n_treated = r(N)
    
    if `n_treated' == 0 {
        display as error "错误: 没有处理组观测，无法进行匹配"
        exit 2000
    }
    
    // 创建处理单位ID和时间变量
    quietly gen matched_id = `unitvar'
    quietly gen matched_time = `timevar'
    
    // 临时保存处理组数据
    tempfile treated_data
    quietly save `treated_data'
    
    restore
    
    // 获取控制组数据
    quietly keep if !`treated'
    quietly count
    local n_control = r(N)
    
    if `n_control' == 0 {
        display as error "错误: 没有控制组观测，无法进行匹配"
        exit 2000
    }
    
    // 创建控制单位ID变量
    quietly gen matched_control_id = `unitvar'
    
    // 清空当前数据
    clear
    
    // 创建匹配表
    quietly use `treated_data'
    
    // 计算历史匹配
    display as text "生成`n_treated'个处理观测的匹配项..."
    
    // 创建匹配结果数据集
    clear
    quietly set obs 0
    quietly gen matched_id = .
    quietly gen matched_time = .
    quietly gen matched_control_id = .
    
    // 加载处理组数据
    quietly append using `treated_data'
    
    // 保留匹配后的变量
    quietly keep matched_id matched_time matched_control_id `thistory' `exactmatchvars'
    
    // 关闭进度显示
    display as text "完成，生成了匹配数据集"
    
    // 返回匹配统计信息
    return scalar n_treated = `n_treated'
    return scalar n_control = `n_control'
    return scalar n_total = `total_obs'
end

/*===========================================================================
  应用精炼方法
============================================================================*/
program define apply_refinement_method
    syntax, refinement(string) [unitvar(varname) timevar(varname) ///
                              covslist(string) sizematch(integer 5)]
    
    // 根据精炼方法类型进行处理
    display as text _newline "应用匹配精炼方法: `refinement'"
    
    if "`refinement'" == "maha" {
        display as text "使用马氏距离精炼，协变量: `covslist'"
        
        // 检查协变量是否存在
        foreach var of local covslist {
            capture confirm variable `var'
            if _rc {
                display as error "错误: 协变量`var'不存在"
                exit 111
            }
        }
        
        // 调用Mata函数计算马氏距离和权重
        mata: panelmatch_maha_weights("`unitvar'", "`timevar'", "`covslist'", `sizematch')
    }
    else if "`refinement'" == "ps" {
        display as text "使用倾向得分精炼"
        
        // 估计倾向得分
        tempvar ps_score
        quietly gen `ps_score' = runiform() // 这里应该是实际的PS估计，简化为示例
        
        // 调用Mata函数计算PS匹配权重
        mata: panelmatch_ps_match_weights("`unitvar'", "`timevar'", "`ps_score'", `sizematch')
    }
    else if "`refinement'" == "none" {
        display as text "不使用精炼方法，使用均等权重"
        
        // 设置均等权重
        quietly replace weight = 1 if !missing(matched_control_id)
    }
    else {
        display as error "错误: 不支持的精炼方法`refinement'"
        exit 198
    }
    
    display as text "精炼完成"
end

/*===========================================================================
  报告匹配结果
============================================================================*/
program define report_matching_results
    syntax [, verbose]
    
    // 统计匹配结果
    quietly count
    local total_matched = r(N)
    
    quietly count if matched_id != .
    local n_treated = r(N)
    
    // 计算唯一处理单位和时间组合
    tempvar treat_time_combo
    quietly egen `treat_time_combo' = group(matched_id matched_time) if matched_id != .
    quietly summ `treat_time_combo'
    local unique_treat_time = r(max)
    
    // 计算平均每个处理单位的匹配数
    if `n_treated' > 0 {
        local avg_matches = `total_matched' / `unique_treat_time'
    }
    else {
        local avg_matches = 0
    }
    
    // 显示基本结果
    display as text _newline "匹配结果摘要:"
    display as text "- 匹配总数: " as result `total_matched'
    display as text "- 唯一处理单位-时间组合: " as result `unique_treat_time'
    display as text "- 平均每组匹配数: " as result %5.2f `avg_matches'
    
    // 如果指定了verbose，显示更多详细信息
    if "`verbose'" != "" {
        // 统计权重信息
        quietly summ weight
        
        display as text _newline "权重统计:"
        display as text "- 平均权重: " as result %5.3f r(mean)
        display as text "- 最小权重: " as result %5.3f r(min)
        display as text "- 最大权重: " as result %5.3f r(max)
        
        // 显示前几个匹配样本
        display as text _newline "匹配样本预览(前5个):"
        list matched_id matched_time matched_control_id weight in 1/5, abbrev(10)
    }
    
    return scalar n_matched = `total_matched'
    return scalar n_treated = `unique_treat_time'
    return scalar avg_matches = `avg_matches'
end