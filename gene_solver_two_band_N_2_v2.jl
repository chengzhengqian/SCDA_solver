# we use the new symbolic solver
# we first enter the development mode
# ]activate "/home/chengzhengqian/share_workspace/czq_julia_package/Wick"
# we add convert in MathExpr.jl, check it now
# using Revise
using Wick
using SSA
using MathExpr
include("./gene_solver_two_band_common.jl")
engine=ExprEngine()
__default__engine__map__["default"]=engine

function proj_full(t,engine::ExprEngine)
    if(t<=2)
        p_t=gene_P_full(t,engine)
        return sum([engine(Symbol("u$(i)"))*p_t[i] for i in 1:18])
    end
end

# the calTree,uniopMap, need to update, we first use as it right now, 
N_orbital,N_spin_orbital,N_time_step,da,input,input_args,total_input_args,ssatape,uniopMap,termMap,calTree=init_para(2,engine)

# uniopMap=initUniOpMap(input,ssatape)

p=proj_full(1,engine)*proj_full(2,engine)
nn=gene_nn(N_time_step,engine)

@time cal(p,"p",calTree)        # first run 0.078836

for i in 1:length(nn)
    cal(p*nn[i],"nn_$(i)",calTree)
end
# i=1;j=1;orb=1
for i in 1:da.N
    for j in 1:da.N
        for orb in 1:N_orbital
            idx_up=defaultSpatialIndex(orb,1)
            idx_down=defaultSpatialIndex(orb,2)
            n_orb_i_j=0.5*(op(da,"dm",[idx_up,i,idx_up,j],engine(1.0))+op(da,"dm",[idx_down,i,idx_down,j],engine(1.0)))
            cal(p*n_orb_i_j,"g_$(orb)_$(i)_$(j)",calTree)
        end
    end
end

outputSyms=Symbol.(["p",["nn_$(i)" for i in 1:4]..., ["g_$(orb)_$(i)_$(j)" for orb in 1:N_orbital  for j in 1:da.N for i in 1:da.N]...])

func_no_diff_v2=SSACompiledFunc(ssatape,outputSyms)
saveSSAFunc(func_no_diff_v2,"./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff_v2")
func_no_diff_v2=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff_v2")
func_no_diff=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff")

test_input=rand(26)
func_no_diff_v2(test_input...)-func_no_diff(test_input...)

# typeof(p)
cal_diff_with_G(p,"p",engine)

for i in 1:length(nn)
    cal_diff_with_G(p*nn[i],"nn_$(i)",engine)
end

for i in 1:da.N
    for j in 1:da.N
        for orb in 1:N_orbital
            idx_up=defaultSpatialIndex(orb,1)
            idx_down=defaultSpatialIndex(orb,2)
            n_orb_i_j=0.5*(op(da,"dm",[idx_up,i,idx_up,j],engine(1.0))+op(da,"dm",[idx_down,i,idx_down,j],engine(1.0)))
            cal_diff_with_G(p*n_orb_i_j,"g_$(orb)_$(i)_$(j)",engine)
        end
    end
end

outputSymsDiff=Symbol.([gene_diff_name("p")...,vcat([gene_diff_name("nn_$(i)") for i in 1:4]...)...,  vcat([gene_diff_name("g_$(orb)_$(i)_$(j)") for orb in 1:N_orbital  for j in 1:da.N for i in 1:da.N]...)...])
outputSymsWithDiff=[outputSyms...,outputSymsDiff...]

func_with_diff_v2=SSACompiledFunc(ssatape,outputSymsWithDiff)
saveSSAFunc(func_with_diff_v2,"./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff_v2")
func_with_diff_v2=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff_v2")


func_with_diff=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff")

test_input=rand(26)
@time result_1=func_with_diff(test_input...)
@time result_2=func_with_diff_v2(test_input...)
sum(abs.(result_1-result_2))

