using Distributed
@everywhere using MutationEntropy
using CairoMakie

"""
计算多个突变体的突变熵(ME)

此函数基于给定的alpha值计算多个突变体的突变熵，并可选地绘制结果
"""
function test_mutants_with_alpha(
    mutations::Vector{String}, 
    data_path::String="data", 
    alpha::Float64=2.0, 
    residue_range=83:231, 
    num_rounds::Int=20, 
    verbose::Bool=true,
    plot_results::Bool=true
)
    # 计算野生型应变值(只需计算一次作为参考)
    println("计算野生型应变值...")
    wt_avg_S_values = MutationEntropy.collect_strains(
        data_path,                # 数据基本目录
        "thermonulease",          # 蛋白质名称
        alpha=alpha,              # 应变计算的alpha参数
        residue_range=residue_range, # 分析的残基范围
        num_rounds=num_rounds,    # 平均的轮次数
        verbose=verbose           # 显示进度信息
    )
    
    # 处理每个突变体
    for (i, mutant) in enumerate(mutations)
        println("处理突变体 $i/$(length(mutations)): $mutant")
        
        # 计算突变体应变值
        mutant_avg_S_values = MutationEntropy.collect_strains(
            data_path,            # 数据基本目录
            mutant,               # 突变名称
            alpha=alpha,          # 应变计算的alpha参数
            residue_range=residue_range, # 分析的残基范围
            num_rounds=num_rounds, # 平均的轮次数
            verbose=verbose       # 显示进度信息
        )
        
        # 计算突变熵
        mutant_ME = MutationEntropy.calculate_ME(mutant_avg_S_values, wt_avg_S_values)
        
        # 如果需要就绘制结果
        if plot_results
            println("为 $mutant 绘制图表...")
            fig = MutationEntropy.plot_MEvsDist(mutant_ME, data_path, mutant)
            
            # 创建figs目录（如果不存在）
            mkpath("figs")
            
            # 保存图表
            save_path = joinpath("figs", "ME_$(mutant)_alpha_$(alpha).png")
            save(save_path, fig)
            println("已保存图表到: $save_path")
        end
    end
    
    println("已完成所有突变体的测试，alpha = $alpha")
end

# 示例: 测试多个突变体
# 定义要测试的突变体列表
mutations = ["a140e", "a140g", "a140v", "a142c", "a142f", "a142g", "a142v", "a151g", "a151t", "a151v"]
# 设置alpha值
alpha = 1.0
# 运行测试
test_mutants_with_alpha(mutations, "data", alpha)
