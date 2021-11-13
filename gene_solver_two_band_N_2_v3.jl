# we first evaluate the coefficient only

# ]activate "/home/chengzhengqian/share_workspace/czq_julia_package/Wick"

# using Revise
using Wick
using SSA
using MathExpr
using JITFunc
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


p=proj_full(1,engine)*proj_full(2,engine)
nn=gene_nn(N_time_step,engine)

p=evalWick(p,ssatape)

@time cal(p,"p",calTree)        # first run 0.078836

for i in 1:length(nn)
    cal(p*nn[i],"nn_$(i)",calTree)
end

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


func_no_diff_v3=SSACompiledFunc(ssatape,outputSyms)
saveSSAFunc(func_no_diff_v3,"./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff_v3")
func_no_diff_v3=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff_v3")
func_no_diff=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff")



test_input=rand(26)
func_no_diff_v3(test_input...)-func_no_diff(test_input...)
