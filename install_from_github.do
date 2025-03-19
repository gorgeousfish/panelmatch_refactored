/*=============================================================================
 从GitHub安装并测试PanelMatch修复版
 
 此脚本用于从GitHub安装PanelMatch修复版并进行基本测试
=============================================================================*/

// 清除工作环境
clear all
set more off

// 开始记录日志
capture log close
log using "install_test.log", replace

// 显示版本信息
display as text "Stata版本: " c(stata_version)
display as text "操作系统: " c(os)
display as text "机器类型: " c(machine_type)
display as text "========================================"

// 从GitHub安装包
display as text _newline "开始从GitHub安装PanelMatch修复版..."

// 使用githubuser仓库
capture noisily net install panelmatch_refactored, from(https://raw.githubusercontent.com/gorgeousfish/panelmatch_refactored/main) replace

if _rc == 0 {
    display as text "PanelMatch修复版安装成功!"
}
else {
    display as error "安装失败，错误代码: " _rc
    exit _rc
}

// 检查命令是否可用
display as text _newline "验证安装的命令..."
which panelmatch
which panelestimate

// 显示采用的寻找路径
display as text "Stata寻找命令的路径:"
adopath

// 创建测试数据
display as text _newline "创建测试数据..."

// 设置种子以确保结果可重复
set seed 12345

// 创建模拟数据
clear
set obs 500
gen country_id = ceil(_n/25)  // 20个国家
gen year = mod(_n-1, 25) + 1990  // 25年(1990-2014)

// 创建经济指标
gen gdp = 100 + 2*year - 1990 + 50*rnormal() + 200*country_id
gen inflation = 5 + 0.5*rnormal()*country_id - 0.2*(year-1990) + rnormal()*2
gen unemployment = 8 + 0.3*rnormal()*country_id - 0.1*(year-1990) + rnormal()*1.5

// 创建处理变量(政策改革)
gen reform_prob = normal(-3 + 0.02*gdp - 0.5*inflation + 0.3*unemployment)
gen reform = 0
replace reform = 1 if reform_prob > 0.7 & year >= 1995

// 确保有处理变化
replace reform = 0 if year < 1995

// 创建结果变量(经济增长率)
gen growth = 3 + 0.0001*gdp - 0.2*inflation - 0.1*unemployment + 1.5*reform + rnormal()*1.2

// 显示数据摘要
summarize
display as text "测试数据创建成功"

// 设置面板数据格式
xtset country_id year

// 测试基本匹配
display as text _newline "测试基本匹配..."
capture noisily {
    panelmatch reform, unitvar(country_id) timevar(year) lag(3) ///
                      refinement(none) savematched("matched_test.dta") replace
}

if _rc == 0 {
    display as text "基本匹配测试成功!"
    use "matched_test.dta", clear
    summarize matched_id matched_time matched_control_id weight
}
else {
    display as error "基本匹配测试失败，错误代码: " _rc
}

// 清理临时文件
capture erase "matched_test.dta"

// 结束日志
display as text _newline "PanelMatch修复版安装和测试完成"
display as text "========================================"
log close