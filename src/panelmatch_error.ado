/*===========================================================================
 PanelMatch错误处理模块 - 专用版
 
 版本: 1.0.0
 日期: 2023-03-04
 作者: [作者信息]
 
 说明: 本文件包含PanelMatch命令中使用的错误处理相关函数
=============================================================================*/

/*===========================================================================
  错误代码定义
============================================================================*/
// 100-199: 输入参数错误
// 200-299: 数据结构错误
// 300-399: 变量错误
// 400-499: 匹配错误
// 500-599: 估计错误
// 900-999: 其他错误

/*===========================================================================
  显示错误信息
============================================================================*/
program define display_error
    args error_code error_message
    
    // 设置错误颜色
    local error_color "error"
    
    // 显示错误信息
    display as `error_color' _newline "错误 [代码 `error_code']: `error_message'"
    
    // 根据错误代码提供更多信息
    if `error_code' >= 100 & `error_code' < 200 {
        display as `error_color' "这是一个输入参数错误。请检查命令语法和参数。"
    }
    else if `error_code' >= 200 & `error_code' < 300 {
        display as `error_color' "这是一个数据结构错误。请检查数据格式和结构。"
    }
    else if `error_code' >= 300 & `error_code' < 400 {
        display as `error_color' "这是一个变量错误。请检查变量名称、类型或值。"
    }
    else if `error_code' >= 400 & `error_code' < 500 {
        display as `error_color' "这是一个匹配错误。请检查匹配条件和选项。"
    }
    else if `error_code' >= 500 & `error_code' < 600 {
        display as `error_color' "这是一个估计错误。请检查模型设置和数据。"
    }
    else {
        display as `error_color' "这是一个一般错误。"
    }
    
    // 显示分隔线
    display as `error_color' "----------------------------------------"
end

/*===========================================================================
  检查变量
============================================================================*/
program define check_variable, rclass
    syntax varname(max=1) [, TYPE(string) RANGE(string) NONmissing OPTional LABEL(string)]
    
    // 提取变量名
    local var `varlist'
    
    // 设置变量标签
    if "`label'" == "" {
        local label "`var'"
    }
    
    // 检查变量是否存在
    capture confirm variable `var'
    if _rc != 0 {
        if "`optional'" != "" {
            display as text "注意: 可选变量 `label' 不存在"
            return scalar exists = 0
            return scalar type_ok = .
            return scalar range_ok = .
            return scalar nonmiss_ok = .
            exit
        }
        else {
            display_error 301 "变量 `label' 不存在"
            exit 111
        }
    }
    
    // 变量存在
    return scalar exists = 1
    
    // 检查变量类型
    if "`type'" != "" {
        local var_type: type `var'
        local type_ok = 0
        
        foreach t in `type' {
            if "`t'" == "numeric" & substr("`var_type'", 1, 3) != "str" {
                local type_ok = 1
            }
            else if "`t'" == "string" & substr("`var_type'", 1, 3) == "str" {
                local type_ok = 1
            }
            else if "`t'" == "`var_type'" {
                local type_ok = 1
            }
            else if "`t'" == "binary" & substr("`var_type'", 1, 3) != "str" {
                qui levelsof `var', local(values)
                local n_values: word count `values'
                if `n_values' <= 2 local type_ok = 1
            }
        }
        
        if `type_ok' == 0 {
            display_error 302 "变量 `label' 类型不正确。期望: `type', 实际: `var_type'"
            exit 108
        }
        
        return scalar type_ok = 1
    }
    else {
        return scalar type_ok = .
    }
    
    // 检查变量范围
    if "`range'" != "" {
        local range_ok = 0
        
        // 仅对数值变量检查范围
        capture confirm numeric variable `var'
        if _rc == 0 {
            qui summarize `var'
            local min = r(min)
            local max = r(max)
            
            // 解析范围表达式
            gettoken min_range max_range: range, parse(",")
            local max_range = subinstr("`max_range'", ",", "", .)
            
            // 检查最小值
            if "`min_range'" != "" & "`min_range'" != "." {
                if `min' < real("`min_range'") {
                    display_error 303 "变量 `label' 的最小值 (`min') 低于允许的最小值 (`min_range')"
                    exit 125
                }
            }
            
            // 检查最大值
            if "`max_range'" != "" & "`max_range'" != "." {
                if `max' > real("`max_range'") {
                    display_error 304 "变量 `label' 的最大值 (`max') 高于允许的最大值 (`max_range')"
                    exit 125
                }
            }
            
            local range_ok = 1
        }
        
        return scalar range_ok = `range_ok'
    }
    else {
        return scalar range_ok = .
    }
    
    // 检查缺失值
    if "`nonmissing'" != "" {
        qui count if missing(`var')
        local n_miss = r(N)
        
        if `n_miss' > 0 {
            display_error 305 "变量 `label' 有 `n_miss' 个缺失值"
            exit 416
        }
        
        return scalar nonmiss_ok = 1
    }
    else {
        return scalar nonmiss_ok = .
    }
end

/*===========================================================================
  检查数据结构
============================================================================*/
program define check_data_structure, rclass
    syntax varlist(min=2 max=2) [, BALanced MINobs(integer 0) PANELsetup]
    
    // 提取变量
    tokenize `varlist'
    local id_var `1'
    local time_var `2'
    
    // 显示数据检查信息
    display as text _newline "检查数据结构..."
    
    // 检查ID和时间变量
    check_variable `id_var', type(numeric) nonmissing
    check_variable `time_var', type(numeric) nonmissing
    
    // 执行基本分析
    preserve
    
    // 计算唯一ID数
    qui duplicates report `id_var'
    local n_ids = r(unique_value)
    
    // 计算唯一时间点数
    qui duplicates report `time_var'
    local n_periods = r(unique_value)
    
    // 检查总观测数
    qui count
    local n_obs = r(N)
    
    // 显示基本信息
    display as text "- 唯一ID数: " as result `n_ids'
    display as text "- 唯一时间点数: " as result `n_periods'
    display as text "- 总观测数: " as result `n_obs'
    
    // 检查是否为平衡面板
    if "`balanced'" != "" {
        display as text "- 检查面板是否平衡..."
        
        tempvar n_obs_per_id
        qui bysort `id_var': gen `n_obs_per_id' = _N
        qui summ `n_obs_per_id'
        local min_obs_per_id = r(min)
        local max_obs_per_id = r(max)
        
        if `min_obs_per_id' != `max_obs_per_id' {
            display_error 201 "面板数据不平衡。每个ID的观测数在 `min_obs_per_id' 到 `max_obs_per_id' 之间"
            exit 460
        }
        else {
            display as text "  ✓ 面板数据平衡，每个ID有 `min_obs_per_id' 个��测"
        }
    }
    
    // 检查最小观测数
    if `minobs' > 0 {
        display as text "- 检查每个ID是否至少有 `minobs' 个观测..."
        
        tempvar n_obs_per_id
        qui bysort `id_var': gen `n_obs_per_id' = _N
        qui count if `n_obs_per_id' < `minobs'
        local n_below_min = r(N)
        
        if `n_below_min' > 0 {
            display_error 202 "有 `n_below_min' 个观测所属的ID观测数少于要求的最小值 `minobs'"
            exit 460
        }
        else {
            display as text "  ✓ 所有ID至少有 `minobs' 个观测"
        }
    }
    
    // 设置面板数据
    if "`panelsetup'" != "" {
        display as text "- 设置面板数据结构..."
        
        qui xtset `id_var' `time_var'
        local balanced = e(balanced)
        local delta = e(tdelta)
        local tmin = e(tmin)
        local tmax = e(tmax)
        local gaps = e(gaps)
        
        display as text "  ✓ 面板数据已设置"
        display as text "  - 面板平衡: " as result "`balanced'"
        display as text "  - 时间间隔: " as result "`delta'"
        display as text "  - 时间范围: " as result "`tmin' 到 `tmax'"
        display as text "  - 有时间缺口: " as result "`gaps'"
    }
    
    restore
    
    // 返回分析结果
    return scalar n_ids = `n_ids'
    return scalar n_periods = `n_periods'
    return scalar n_obs = `n_obs'
end