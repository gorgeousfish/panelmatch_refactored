/*===========================================================================
 PanelMatch命令 - 重构版
 
 版本: 1.0.1
 日期: 2025-03-04
 作者: [作者信息]
 
 说明: 本程序实现了面板数据匹配方法，用于因果推断分析
=============================================================================*/

/*===========================================================================
  PanelMatch主命令
============================================================================*/
program define panelmatch, eclass

    // 初始化，确保辅助函数可用 - 使用adopath机制
    cap program drop validate_panelmatch_opts
    cap program drop generate_treatment_history
    cap program drop prepare_matching_data
    cap program drop generate_matched_dataset
    cap program drop apply_refinement_method
    cap program drop report_matching_results

    // 定义命令语法
    syntax varlist(min=1 max=1) [if] [in], ///
        UNitvar(varname) TImevar(varname) ///
        Lag(integer) ///
        [ REFinement(string) SIZEmatch(integer 5) COVSformula(string) ///
          MATCHmissing QOI(string) LEAD(numlist >0 integer) OUTCome(varname) ///
          FORBIDtreatmentreversal EXACTmatchvars(varlist) LISTWISEdelete ///
          VERBose SAVEMatched(string) REPlace ]
    
    // 设置默认值
    if "`refinement'" == "" local refinement none
    if "`qoi'" == "" local qoi att
    if "`sizematch'" == "" local sizematch 5
    
    // 检查必填选项
    if "`unitvar'" == "" | "`timevar'" == "" {
        display as error "必须指定unitvar和timevar选项。"
        exit 198
    }
    
    // 验证选项值
    validate_panelmatch_opts, refinement(`refinement') qoi(`qoi')
    
    // 显示PanelMatch版本和选项信息
    display as text _newline "执行PanelMatch命令 (v1.0.1)"
    display as text "================================================="
    display as text "参数设置:"
    display as text "- 单位ID变量: " as result "`unitvar'"
    display as text "- 时间变量: " as result "`timevar'"
    display as text "- 处理变量: " as result "`varlist'"
    display as text "- 滞后期数: " as result "`lag'"
    display as text "- 精炼方法: " as result "`refinement'"
    if "`covsformula'" != "" {
        display as text "- 协变量公式: " as result "`covsformula'"
    }
    display as text "- 最大匹配数: " as result "`sizematch'"
    if "`lead'" != "" {
        display as text "- Lead值: " as result "`lead'"
    }
    if "`outcome'" != "" {
        display as text "- 结果变量: " as result "`outcome'"
    }
    display as text "================================================="
    
    // 标记要使用的样本
    marksample touse
    
    // 提取处理变量
    local treatvar `varlist'
    
    // 检查数据集状态
    qui describe, short
    if r(changed) {
        display as text "注意: 当前数据集已修改但未保存"
    }
    if r(N) == 0 {
        display as error "错误: 当前数据集为空"
        exit 2000
    }
    
    // 生成处理历史变量
    tempvar thistory
    generate_treatment_history `treatvar', unitvar(`unitvar') timevar(`timevar') ///
                            lag(`lag') touse(`touse') hist_name(`thistory') replace
    
    // 准备匹配数据
    prepare_matching_data `treatvar', unitvar(`unitvar') timevar(`timevar') ///
                          thistory(`thistory') touse(`touse') ///
                          `listwise' `matchmissing' `forbidtreatmentreversal' ///
                          exact_match(`exactmatchvars')
    
    // 进行匹配并生成匹配数据集
    preserve
    
    // 生成匹配结果
    generate_matched_dataset `treatvar', unitvar(`unitvar') timevar(`timevar') ///
                            thistory(`thistory') treatvar(`treatvar') ///
                            exactmatchvars(`exactmatchvars') touse(`touse')
    
    // 添加权重变量
    quietly gen weight = 1
    
    // 应用精炼方法计算权重
    if "`refinement'" == "maha" {
        if "`covsformula'" == "" {
            display as error "使用马氏距离匹配时必须指定covsformula选项。"
            exit 198
        }
        
        // 提取协变量列表
        local covs_list = subinstr("`covsformula'", "+", " ", .)
        
        // 应用马氏距离精炼
        apply_refinement_method, refinement(`refinement') unitvar(`unitvar') ///
                               timevar(`timevar') covslist(`covs_list') ///
                               sizematch(`sizematch')
    }
    else if "`refinement'" == "ps" {
        // 应用倾向得分精炼
        apply_refinement_method, refinement(`refinement') unitvar(`unitvar') ///
                               timevar(`timevar') sizematch(`sizematch')
    }
    else {
        // 使用无精炼方法(均等权重)
        apply_refinement_method, refinement(`refinement')
    }
    
    // 报告匹配结果
    report_matching_results, `verbose'
    
    // 保存匹配结果(如果指定)
    if "`savematched'" != "" {
        // 检查是否存在同名文件
        capture confirm file "`savematched'.dta"
        if !_rc & "`replace'" == "" {
            display as error "文件 `savematched'.dta 已存在。使用replace选项覆盖。"
            exit 602
        }
        
        // 保存数据
        quietly save "`savematched'", `replace'
        
        if "`verbose'" != "" {
            display as text "匹配结果已保存到文件: `savematched'.dta"
        }
    }
    
    // 存储结果
    return scalar n_matched = r(N)
    
    // 如果有lead和outcome选项，计算每个lead期的效应
    if "`lead'" != "" & "`outcome'" != "" {
        // 准备lead列表
        numlist "`lead'"
        local lead_list `r(numlist)'
        local n_leads: word count `lead_list'
        
        // 创建结果矩阵
        matrix b = J(1, `n_leads', .)
        matrix V = J(`n_leads', `n_leads', 0)
        matrix colnames b = `lead_list'
        matrix rownames V = `lead_list'
        matrix colnames V = `lead_list'
        
        // 对每个lead计算效应
        local i = 1
        foreach l of local lead_list {
            // 创建lead结果变量
            tempvar lead_outcome
            quietly gen `lead_outcome' = .
            
            // 填充lead结果
            quietly {
                // 获取处理单位在lead期的结果
                foreach v of varlist matched_id matched_time {
                    sort `v'
                }
                
                levelsof matched_id, local(ids)
                levelsof matched_time, local(times)
                
                // 对每个处理单位计算lead结果
                foreach id of local ids {
                    foreach t of local times {
                        if matched_id == `id' & matched_time == `t' {
                            // 寻找lead期的结果
                            local lead_t = `t' + `l'
                            
                            // 在原始数据中查找lead结果
                            sum `outcome' if `unitvar' == `id' & `timevar' == `lead_t' & `touse'
                            
                            if r(N) > 0 {
                                replace `lead_outcome' = r(mean) ///
                                    if matched_id == `id' & matched_time == `t'
                            }
                        }
                    }
                }
                
                // 对控制单位也做相同操作
                foreach id of local ids {
                    foreach t of local times {
                        if matched_id == `id' & matched_time == `t' {
                            // 获取控制单位ID
                            local control_id = matched_control_id
                            
                            // 寻找lead期的结果
                            local lead_t = `t' + `l'
                            
                            // 在原始数据中查找lead结果
                            sum `outcome' if `unitvar' == `control_id' & `timevar' == `lead_t' & `touse'
                            
                            if r(N) > 0 {
                                replace `lead_outcome' = r(mean) ///
                                    if matched_control_id == `control_id' & matched_time == `t'
                            }
                        }
                    }
                }
            }
            
            // 计算该lead期的效应
            tempname att se
            mata: calculate_att("matched_id", "matched_time", "matched_control_id", ///
                              "weight", "`lead_outcome'", "`att'", "`se'")
            
            // 存储结果
            matrix b[1, `i'] = scalar(`att')
            matrix V[`i', `i'] = scalar(`se')^2
            
            local ++i
        }
        
        // 存储面板匹配结果
        ereturn post b V
        ereturn local cmd "panelmatch"
        ereturn local cmdline `"panelmatch `0'"'
        ereturn local treatvar "`treatvar'"
        ereturn local outcome "`outcome'"
        ereturn local refinement "`refinement'"
        ereturn local lead "`lead'"
        ereturn scalar lag = `lag'
        ereturn scalar n_matched = r(N)
        
        // 显示结果
        display as text _newline "Panel匹配效应估计结果 (QOI: `qoi')" _continue
        display as text " - 使用 `refinement' 精炼方法"
        display as text "处理变量: `treatvar'  结果变量: `outcome'"
        display as text "样本大小: " as result e(n_matched)
        display as text _newline "估计效应:"
        
        foreach l of local lead_list {
            local i_col = colnumb(b, "`l'")
            local est = b[1, `i_col']
            local se = sqrt(V[`i_col', `i_col'])
            local t = `est' / `se'
            local p = 2 * normal(-abs(`t'))
            
            display as text "Lead `l': " _continue
            display as result %9.6f `est' _continue
            display as text " (SE: " _continue
            display as result %9.6f `se' _continue
            display as text ", t: " _continue
            display as result %8.3f `t' _continue
            display as text ", p: " _continue
            display as result %8.3f `p' _continue
            display as text ")"
        }
    }
    
    restore
end 