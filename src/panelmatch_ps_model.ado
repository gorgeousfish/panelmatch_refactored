/*===========================================================================
 PanelMatch倾向得分模型函数 - 专用版
 
 版本: 1.0.0
 日期: 2025-03-04
 作者: [作者信息]
 
 说明: 本文件包含PanelMatch命令中使用的倾向得分模型估计相关函数
=============================================================================*/

/*===========================================================================
  倾向得分模型估计
============================================================================*/
program define estimate_ps_model
    syntax varlist(min=1 max=1), UNitvar(varname) TImevar(varname) ///
                             [PSCovs(varlist) PSFormula(string) ///
                             PSMODel(string) Touse(varname) ///
                             PSOUtput(name) REPlace SAVEps(name)]
    
    // 提取变量名
    local treatvar `varlist'
    
    // 准备touse变量
    if "`touse'" == "" {
        tempvar touse
        quietly gen byte `touse' = 1
    }
    
    // 确定模型类型
    if "`psmodel'" == "" {
        local psmodel "logit"
    }
    
    // 检查模型类型是否有效
    if !inlist("`psmodel'", "logit", "probit", "cloglog") {
        display as error "倾向得分模型必须是以下之一: logit, probit, cloglog"
        exit 198
    }
    
    // 确定协变量
    if "`pscovs'" == "" & "`psformula'" == "" {
        display as error "必须指定pscovs或psformula"
        exit 198
    }
    
    // 如果同时指定了两者，优先使用psformula
    if "`psformula'" != "" & "`pscovs'" != "" {
        display as text "注意: 同时指定了psformula和pscovs，将使用psformula"
    }
    
    // 确定倾向得分变量名
    if "`psoutput'" == "" {
        tempvar psoutput
    }
    else {
        // 检查变量是否存在
        capture confirm variable `psoutput'
        if !_rc & "`replace'" == "" {
            display as error "变量`psoutput'已存在，请使用replace选项或指定其他变量名"
            exit 110
        }
        else if !_rc {
            drop `psoutput'
        }
    }
    
    // 开始估计
    display as text _newline "估计倾向得分模型..."
    display as text "- 模型类型: " as result "`psmodel'"
    if "`psformula'" != "" {
        display as text "- 使用公式: " as result "`psformula'"
    }
    else {
        display as text "- 使用协变量: " as result "`pscovs'"
    }
    
    // 设置模型公式
    local model_formula `treatvar'
    if "`psformula'" != "" {
        local model_formula `model_formula' `psformula'
    }
    else {
        local model_formula `model_formula' `pscovs'
    }
    
    // 拟合模型
    if "`psmodel'" == "logit" {
        quietly logit `model_formula' if `touse'
    }
    else if "`psmodel'" == "probit" {
        quietly probit `model_formula' if `touse'
    }
    else if "`psmodel'" == "cloglog" {
        quietly cloglog `model_formula' if `touse'
    }
    
    // 检查模型是否收敛
    if e(converged) != 1 {
        display as error "倾向得分模型未收敛，请检查模型设置或使用不同的模型"
        exit 430
    }
    
    // 显示模型结果
    quietly estimates store ps_model
    estimates table ps_model, b(%9.4f) se(%9.4f) stats(N r2_p)
    
    // 计算并保存倾向得分
    quietly predict `psoutput' if `touse', pr
    
    // 如果saveps选项被指定，创建新变量
    if "`saveps'" != "" {
        capture confirm variable `saveps'
        if !_rc & "`replace'" != "" {
            quietly replace `saveps' = `psoutput' if `touse'
            display as text "倾向得分已保存到变量: `saveps' (已替换)"
        }
        else if !_rc {
            display as error "变量`saveps'已存在，请使用replace选项或指定其他变量名"
            exit 110
        }
        else {
            quietly gen `saveps' = `psoutput' if `touse'
            display as text "倾向得分已保存到变量: `saveps'"
        }
    }
    
    // 汇总倾向得分
    quietly summarize `psoutput' if `touse'
    display as text _newline "倾向得分摘要:"
    display as text "- 观测数: " as result r(N)
    display as text "- 平均值: " as result %9.4f r(mean)
    display as text "- 标准差: " as result %9.4f r(sd)
    display as text "- 最小值: " as result %9.4f r(min)
    display as text "- 最大值: " as result %9.4f r(max)
    
    // 返回模型信息
    return scalar N = e(N)
    return scalar r2_p = e(r2_p)
    return local cmd = e(cmd)
    return local depvar = e(depvar)
    return local psvariable = "`psoutput'"
end

/*===========================================================================
  检查倾向得分公共支撑
============================================================================*/
program define check_ps_overlap, rclass
    syntax varname(max=1), Treatvar(varname) [Touse(varname) ///
                               THREshold(real 0.1) GRaph NOLabel]
    
    // 提取变量名
    local psvar `varlist'
    
    // 准备touse变量
    if "`touse'" == "" {
        tempvar touse
        quietly gen byte `touse' = 1
    }
    
    // 显示标题
    display as text _newline "检查倾向得分公共支撑..."
    
    // 按处理状态计算倾向得分分布
    tempvar ps_treated ps_control
    quietly gen `ps_treated' = `psvar' if `treatvar' == 1 & `touse'
    quietly gen `ps_control' = `psvar' if `treatvar' == 0 & `touse'
    
    // 计算处理组和对照组的PS分位数
    quietly summarize `ps_treated' if `touse', detail
    local t_p01 = r(p1)
    local t_p05 = r(p5)
    local t_p10 = r(p10)
    local t_p25 = r(p25)
    local t_p50 = r(p50)
    local t_p75 = r(p75)
    local t_p90 = r(p90)
    local t_p95 = r(p95)
    local t_p99 = r(p99)
    local t_mean = r(mean)
    local t_sd = r(sd)
    local t_min = r(min)
    local t_max = r(max)
    local t_n = r(N)
    
    quietly summarize `ps_control' if `touse', detail
    local c_p01 = r(p1)
    local c_p05 = r(p5)
    local c_p10 = r(p10)
    local c_p25 = r(p25)
    local c_p50 = r(p50)
    local c_p75 = r(p75)
    local c_p90 = r(p90)
    local c_p95 = r(p95)
    local c_p99 = r(p99)
    local c_mean = r(mean)
    local c_sd = r(sd)
    local c_min = r(min)
    local c_max = r(max)
    local c_n = r(N)
    
    // 显示倾向得分分布
    display as text _newline "倾向得分分布:"
    display as text "----------------------------------------"
    display as text "      |  处理组  |  控制组  | 差异"
    display as text "----------------------------------------"
    display as text "样本量 | " _continue
    display as result %8.0f `t_n' _continue
    display as text " | " _continue
    display as result %8.0f `c_n' _continue
    display as text " | " _continue
    display as result %8.0f `t_n' - `c_n'
    
    display as text "均值  | " _continue
    display as result %8.4f `t_mean' _continue
    display as text " | " _continue
    display as result %8.4f `c_mean' _continue
    display as text " | " _continue
    display as result %8.4f `t_mean' - `c_mean'
    
    display as text "标准差 | " _continue
    display as result %8.4f `t_sd' _continue
    display as text " | " _continue
    display as result %8.4f `c_sd' _continue
    display as text " | " _continue
    display as result %8.4f `t_sd' - `c_sd'
    
    display as text "最小值 | " _continue
    display as result %8.4f `t_min' _continue
    display as text " | " _continue
    display as result %8.4f `c_min' _continue
    display as text " | " _continue
    display as result %8.4f `t_min' - `c_min'
    
    display as text "p10   | " _continue
    display as result %8.4f `t_p10' _continue
    display as text " | " _continue
    display as result %8.4f `c_p10' _continue
    display as text " | " _continue
    display as result %8.4f `t_p10' - `c_p10'
    
    display as text "p25   | " _continue
    display as result %8.4f `t_p25' _continue
    display as text " | " _continue
    display as result %8.4f `c_p25' _continue
    display as text " | " _continue
    display as result %8.4f `t_p25' - `c_p25'
    
    display as text "p50   | " _continue
    display as result %8.4f `t_p50' _continue
    display as text " | " _continue
    display as result %8.4f `c_p50' _continue
    display as text " | " _continue
    display as result %8.4f `t_p50' - `c_p50'
    
    display as text "p75   | " _continue
    display as result %8.4f `t_p75' _continue
    display as text " | " _continue
    display as result %8.4f `c_p75' _continue
    display as text " | " _continue
    display as result %8.4f `t_p75' - `c_p75'
    
    display as text "p90   | " _continue
    display as result %8.4f `t_p90' _continue
    display as text " | " _continue
    display as result %8.4f `c_p90' _continue
    display as text " | " _continue
    display as result %8.4f `t_p90' - `c_p90'
    
    display as text "最大值 | " _continue
    display as result %8.4f `t_max' _continue
    display as text " | " _continue
    display as result %8.4f `c_max' _continue
    display as text " | " _continue
    display as result %8.4f `t_max' - `c_max'
    display as text "----------------------------------------"
    
    // 检查公共支撑
    local min_treated = max(`t_min', `threshold')
    local max_treated = min(`t_max', 1 - `threshold')
    local min_control = max(`c_min', `threshold')
    local max_control = min(`c_max', 1 - `threshold')
    
    local support_min = max(`min_treated', `min_control')
    local support_max = min(`max_treated', `max_control')
    
    display as text _newline "公共支撑区间(阈值=`threshold'):"
    display as text "- 处理组: [" as result %6.4f `min_treated' as text ", " as result %6.4f `max_treated' as text "]"
    display as text "- 控制组: [" as result %6.4f `min_control' as text ", " as result %6.4f `max_control' as text "]"
    display as text "- 公共支撑: [" as result %6.4f `support_min' as text ", " as result %6.4f `support_max' as text "]"
    
    // 计算落在公共支撑外的观测
    tempvar outside_support
    quietly gen `outside_support' = (`psvar' < `support_min' | `psvar' > `support_max') if `touse'
    
    quietly count if `outside_support' == 1 & `treatvar' == 1 & `touse'
    local n_out_treated = r(N)
    local pct_out_treated = 100 * `n_out_treated' / `t_n'
    
    quietly count if `outside_support' == 1 & `treatvar' == 0 & `touse'
    local n_out_control = r(N)
    local pct_out_control = 100 * `n_out_control' / `c_n'
    
    quietly count if `outside_support' == 1 & `touse'
    local n_out_total = r(N)
    local pct_out_total = 100 * `n_out_total' / (`t_n' + `c_n')
    
    display as text _newline "公共支撑外的观测:"
    display as text "- 处理组: " as result `n_out_treated' as text " 个观测 (" as result %5.1f `pct_out_treated' as text "%)"
    display as text "- 控制组: " as result `n_out_control' as text " 个观测 (" as result %5.1f `pct_out_control' as text "%)"
    display as text "- 总计: " as result `n_out_total' as text " 个观测 (" as result %5.1f `pct_out_total' as text"%)"
    
    // 如果请求图形，绘制PS分布图
    if "`graph'" != "" {
        tempvar treat_label
        if "`nolabel'" == "" {
            quietly gen `treat_label' = "控制组" if `treatvar' == 0
            quietly replace `treat_label' = "处理组" if `treatvar' == 1
        }
        else {
            quietly gen `treat_label' = "0" if `treatvar' == 0
            quietly replace `treat_label' = "1" if `treatvar' == 1
        }
        
        // 绘制直方图
        histogram `psvar' if `touse', by(`treat_label', note("")) ///
            start(0) width(0.05) fcol(navy) lcol(white) ///
            xtitle("倾向得分") ytitle("密度") ///
            title("倾向得分分布") subtitle("按处理状态分组")
    }
    
    // 返回结果
    return scalar n_treated = `t_n'
    return scalar n_control = `c_n'
    return scalar p_treated = `t_mean'
    return scalar p_control = `c_mean'
    return scalar support_min = `support_min'
    return scalar support_max = `support_max'
    return scalar n_outside = `n_out_total'
    return scalar pct_outside = `pct_out_total'
end