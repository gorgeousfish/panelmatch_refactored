/*===========================================================================
 PanelMatch Mata函数集 - 重构版
 
 版本: 1.0.1
 日期: 2025-03-04
 作者: [作者信息]
 
 说明: 本文件包含PanelMatch和PanelEstimate命令使用的Mata函数
=============================================================================*/

// 启动mata环境并设置版本兼容性
version 14.0  // 降低版本要求，确保兼容更多Stata版本
mata:

/*===========================================================================
  马氏距离计算和权重分配函数
===========================================================================*/
/**
 * panelmatch_maha_weights()函数 - 计算基于马氏距离的匹配权重
 * 
 * 参数:
 *   unit_id - 单位ID变量名
 *   time_id - 时间ID变量名
 *   covs_list - 协变量列表(空格分隔)
 *   size_match - 每个处理单位最多匹配的控制单位数量
 */
void panelmatch_maha_weights(string scalar unit_id, string scalar time_id, 
                 string scalar covs_list, real scalar size_match)
{
    // 将协变量列表转换为向量
    string vector cov_vars
    cov_vars = tokens(covs_list)
    
    // 初始化各变量视图
    real vector matched_id, matched_time, matched_control_id
    st_view(matched_id=., ., "matched_id")
    st_view(matched_time=., ., "matched_time")
    st_view(matched_control_id=., ., "matched_control_id")
    
    // 创建权重变量(如果不存在)
    if (!st_varindex("weight")) st_addvar("float", "weight")
    real vector weights
    st_view(weights=., ., "weight")
    
    // 获取唯一的处理单位-时间组合
    real matrix unique_matched
    unique_matched = uniqrows((matched_id, matched_time))
    
    // 初始化权重为0
    weights[.] = J(rows(weights), 1, 0)
    
    // 对每个处理单位-时间组计算马氏距离和权重
    real scalar num_groups
    num_groups = rows(unique_matched)
    
    for (g=1; g<=num_groups; g++) {
        // 当前处理单位和时间
        real scalar current_id, current_time
        current_id = unique_matched[g,1]
        current_time = unique_matched[g,2]
        
        // 找到当前组的行索引
        real vector group_rows
        group_rows = selectindex((matched_id :== current_id) :& (matched_time :== current_time))
        
        if (rows(group_rows) == 0) continue
        
        // 获取处理单位和对应控制单位的协变量
        real scalar treated_unit_id, treated_unit_time
        treated_unit_id = current_id
        treated_unit_time = current_time
        
        // 获取处理单位的协变量
        real matrix treated_covs
        treated_covs = J(1, cols(cov_vars), .)
        
        for (j=1; j<=cols(cov_vars); j++) {
            string scalar treat_var_name
            treat_var_name = cov_vars[j]
            
            real colvector var_data
            st_view(var_data=., ., treat_var_name)
            
            real colvector unit_data, time_data, where
            st_view(unit_data=., ., unit_id)
            st_view(time_data=., ., time_id)
            
            // 找到处理单位对应行
            where = selectindex((unit_data :== treated_unit_id) :& (time_data :== treated_unit_time))
            
            if (rows(where) > 0) {
                treated_covs[1,j] = var_data[where[1]]
            }
        }
        
        // 获取控制单位的协变量和ID
        real scalar num_controls
        num_controls = rows(group_rows)
        
        real matrix control_ids, control_covs
        control_ids = J(num_controls, 1, .)
        control_covs = J(num_controls, cols(cov_vars), .)
        
        real scalar valid_controls
        valid_controls = 0
        
        for (i=1; i<=num_controls; i++) {
            real scalar control_unit_id
            control_unit_id = matched_control_id[group_rows[i]]
            
            real colvector has_all_covs
            has_all_covs = J(1, 1, 1)
            
            for (j=1; j<=cols(cov_vars); j++) {
                string scalar control_var_name
                control_var_name = cov_vars[j]
                
                real colvector ctrl_var_data, ctrl_unit_data, ctrl_time_data, ctrl_where
                st_view(ctrl_var_data=., ., control_var_name)
                st_view(ctrl_unit_data=., ., unit_id)
                st_view(ctrl_time_data=., ., time_id)
                
                // 找到控制单位对应行
                ctrl_where = selectindex((ctrl_unit_data :== control_unit_id) :& (ctrl_time_data :== current_time))
                
                if (rows(ctrl_where) > 0) {
                    control_covs[valid_controls+1, j] = ctrl_var_data[ctrl_where[1]]
                }
                else {
                    has_all_covs[1] = 0
                    break
                }
            }
            
            if (has_all_covs[1]) {
                valid_controls = valid_controls + 1
                control_ids[valid_controls] = control_unit_id
            }
        }
        
        if (valid_controls == 0) continue
        
        // 裁剪控制单位矩阵为有效大小
        if (valid_controls < num_controls) {
            control_ids = control_ids[1::valid_controls]
            control_covs = control_covs[1::valid_controls,]
        }
        
        // 计算处理单位和控制单位协变量差异
        real matrix diffs
        diffs = J(valid_controls, cols(cov_vars), .)
        
        for (i=1; i<=valid_controls; i++) {
            diffs[i,] = control_covs[i,] :- treated_covs
        }
        
        // 计算协变量协方差矩阵
        real matrix S
        S = quadcross(diffs, diffs) :/ (valid_controls - 1)
        
        // 添加一个小的对角线增量以确保可逆
        real scalar epsilon
        epsilon = 1e-10
        
        for (j=1; j<=cols(S); j++) {
            S[j,j] = S[j,j] + epsilon
        }
        
        // 计算马氏距离
        real matrix S_inv, m_dists
        S_inv = invsym(S)
        m_dists = J(valid_controls, 1, .)
        
        for (i=1; i<=valid_controls; i++) {
            m_dists[i] = sqrt((diffs[i,] * S_inv * diffs[i,]'))
        }
        
        // 根据距离排序
        real vector sorted_indexes
        sorted_indexes = order(m_dists, 1)
        
        // 选择最接近的size_match个控制单位
        real scalar n_matches
        n_matches = min((size_match, valid_controls))
        
        // 分配权重
        for (i=1; i<=n_matches; i++) {
            real scalar idx
            idx = group_rows[sorted_indexes[i]]
            weights[idx] = 1.0 / n_matches
        }
    }
}

/*===========================================================================
  倾向得分匹配函数
===========================================================================*/
/**
 * panelmatch_ps_match_weights()函数 - 计算基于倾向得分的匹配权重
 * 
 * 参数:
 *   unit_id - 单位ID变量名
 *   time_id - 时间ID变量名
 *   ps_score - 倾向得分变量名
 *   size_match - 每个处理单位最多匹配的控制单位数量
 */
void panelmatch_ps_match_weights(string scalar unit_id, string scalar time_id, 
                    string scalar ps_score, real scalar size_match)
{
    // 初始化各变量视图
    real vector matched_id, matched_time, matched_control_id
    st_view(matched_id=., ., "matched_id")
    st_view(matched_time=., ., "matched_time")
    st_view(matched_control_id=., ., "matched_control_id")
    
    // 创建权重变量(如果不存在)
    if (!st_varindex("weight")) st_addvar("float", "weight")
    real vector weights
    st_view(weights=., ., "weight")
    
    // 获取唯一的处理单位-时间组合
    real matrix unique_matched
    unique_matched = uniqrows((matched_id, matched_time))
    
    // 初始化权重为0
    weights[.] = J(rows(weights), 1, 0)
    
    // 对每个处理单位-时间组计算PS距离和权重
    real scalar num_groups
    num_groups = rows(unique_matched)
    
    for (g=1; g<=num_groups; g++) {
        // 当前处理单位和时间
        real scalar current_id, current_time
        current_id = unique_matched[g,1]
        current_time = unique_matched[g,2]
        
        // 找到当前组的行索引
        real vector group_rows
        group_rows = selectindex((matched_id :== current_id) :& (matched_time :== current_time))
        
        if (rows(group_rows) == 0) continue
        
        // 获取处理单位的PS
        real scalar treated_unit_id, treated_unit_time, treated_ps
        treated_unit_id = current_id
        treated_unit_time = current_time
        
        // 获取PS值
        real colvector ps_data, unit_data, time_data, where
        st_view(ps_data=., ., ps_score)
        st_view(unit_data=., ., unit_id)
        st_view(time_data=., ., time_id)
        
        // 找到处理单位对应行
        where = selectindex((unit_data :== treated_unit_id) :& (time_data :== treated_unit_time))
        
        if (rows(where) == 0) continue
        
        treated_ps = ps_data[where[1]]
        
        // 获取控制单位的PS
        real scalar num_controls
        num_controls = rows(group_rows)
        
        real matrix control_ids, control_ps
        control_ids = J(num_controls, 1, .)
        control_ps = J(num_controls, 1, .)
        
        real scalar valid_controls
        valid_controls = 0
        
        for (i=1; i<=num_controls; i++) {
            real scalar control_unit_id
            control_unit_id = matched_control_id[group_rows[i]]
            
            // 找到控制单位对应行
            real colvector ctrl_where
            ctrl_where = selectindex((unit_data :== control_unit_id) :& (time_data :== current_time))
            
            if (rows(ctrl_where) > 0) {
                valid_controls = valid_controls + 1
                control_ids[valid_controls] = control_unit_id
                control_ps[valid_controls] = ps_data[ctrl_where[1]]
            }
        }
        
        if (valid_controls == 0) continue
        
        // 裁剪控制单位矩阵为有效大小
        if (valid_controls < num_controls) {
            control_ids = control_ids[1::valid_controls]
            control_ps = control_ps[1::valid_controls]
        }
        
        // 计算PS距离
        real vector ps_dists
        ps_dists = abs(control_ps :- treated_ps)
        
        // 根据PS距离排序
        real vector sorted_indexes
        sorted_indexes = order(ps_dists, 1)
        
        // 选择PS最接近的size_match个控制单位
        real scalar n_matches
        n_matches = min((size_match, valid_controls))
        
        // 分配权重
        for (i=1; i<=n_matches; i++) {
            real scalar idx
            idx = group_rows[sorted_indexes[i]]
            weights[idx] = 1.0 / n_matches
        }
    }
}

/*===========================================================================
  处理效应计算函数
===========================================================================*/
/**
 * panelmatch_calculate_att()函数 - 计算处理效应(ATT)和标准误
 * 
 * 参数:
 *   matched_id_var - 处理单位ID变量名
 *   matched_time_var - 处理时间变量名
 *   matched_control_id_var - 匹配的控制单位ID变量名
 *   weight_var - 权重变量名
 *   outcome_var - 结果变量名
 *   att_name - 返回ATT估计值的局部宏名称
 *   se_name - 返回标准误的局部宏名称
 */
void panelmatch_calculate_att(string scalar matched_id_var, 
                 string scalar matched_time_var,
                 string scalar matched_control_id_var, 
                 string scalar weight_var,
                 string scalar outcome_var,
                 string scalar att_name,
                 string scalar se_name)
{
    // 初始化各变量视图
    real vector matched_id, matched_time, matched_control_id, weights, outcome
    st_view(matched_id=., ., matched_id_var)
    st_view(matched_time=., ., matched_time_var)
    st_view(matched_control_id=., ., matched_control_id_var)
    st_view(weights=., ., weight_var)
    st_view(outcome=., ., outcome_var)
    
    // 获取唯一的处理单位-时间组合
    real matrix unique_matched
    unique_matched = uniqrows((matched_id, matched_time))
    
    // 初始化结果
    real scalar num_groups, att, att_var
    num_groups = rows(unique_matched)
    att = 0
    att_var = 0
    
    // 计算每个组的处理效应差异
    real vector group_effects, group_weights
    group_effects = J(num_groups, 1, .)
    group_weights = J(num_groups, 1, .)
    
    for (g=1; g<=num_groups; g++) {
        // 当前处理单位和时间
        real scalar current_id, current_time
        current_id = unique_matched[g,1]
        current_time = unique_matched[g,2]
        
        // 找到当前组的行索引
        real vector group_rows
        group_rows = selectindex((matched_id :== current_id) :& (matched_time :== current_time))
        
        if (rows(group_rows) == 0) continue
        
        // 获取处理单位的结果
        real scalar treated_outcome
        treated_outcome = outcome[group_rows[1]]
        
        // 获取加权控制单位结果
        real scalar ctrl_sum, weight_sum
        ctrl_sum = 0
        weight_sum = 0
        
        for (i=1; i<=rows(group_rows); i++) {
            real scalar row_idx, control_weight
            row_idx = group_rows[i]
            control_weight = weights[row_idx]
            
            if (control_weight > 0) {
                ctrl_sum = ctrl_sum + outcome[row_idx] * control_weight
                weight_sum = weight_sum + control_weight
            }
        }
        
        // 如果有有效的控制单位
        if (weight_sum > 0) {
            real scalar weighted_ctrl_outcome, effect
            weighted_ctrl_outcome = ctrl_sum / weight_sum
            effect = treated_outcome - weighted_ctrl_outcome
            
            group_effects[g] = effect
            group_weights[g] = 1  // 等权重，可以基于处理单位的特性调整
        }
        else {
            group_effects[g] = .
            group_weights[g] = 0
        }
    }
    
    // 过滤掉缺失值
    real vector valid_idx
    valid_idx = selectindex(group_weights :> 0)
    
    if (rows(valid_idx) == 0) {
        st_local(att_name, "0")
        st_local(se_name, ".")
        return
    }
    
    group_effects = group_effects[valid_idx]
    group_weights = group_weights[valid_idx]
    
    // 标准化权重
    group_weights = group_weights :/ sum(group_weights)
    
    // 计算加权平均处理效应
    att = sum(group_effects :* group_weights)
    
    // 计算效应变异
    real vector demeaned_effects
    demeaned_effects = group_effects :- att
    
    real scalar se
    se = sqrt(sum((demeaned_effects:^2) :* group_weights) / (rows(valid_idx) - 1))
    
    // 返回结果
    st_local(att_name, strofreal(att))
    st_local(se_name, strofreal(se))
}

// 结束mata环境
end