# run on cluster

# !! to run the code, one should use
# ulimit -s unlimited

using Wick
using SymEngine
using CZQUtils

include("./gene_solver_two_band_common.jl")

"""
Notice that this is different from N=2,3
as the SPD has the following structure
P1*n2*P2*n3*P2*n2*P4
so momentum [0.5,n2,n3,n2]      #  similar to N=3 as [n1,n2,n1]
but local [P1,P2,P2,P1]
So we use u_i_1, ..,u_i_18  to store teh data
proj, we use a different way to parameterize local projector
using u1,..,u18
# here, we just have three time step, but how only have P1 and P2, as P3=1
!!!!!!!!!!!!!!!!!Notice u_1_.. is the outside
"""
function proj_full(t)
    if(t==1 || t==4)
        tag=1
    elseif(t==2 || t==3)
        tag=2
    else
        error("t should be in 1...4\n")
    end
    p_t=gene_P_full(t)
    return simplify(sum([symbols("u_$(tag)_$(i)")*p_t[i] for i in 1:18]))
end


N_orbital,N_spin_orbital,N_time_step,da,input,input_args,total_input_args,ssatape,uniopMap,calTree=init_para(4)

# we first evaluate coefficient
p=proj_full(1)*proj_full(2)*proj_full(3)*proj_full(4)
#p=simplify(p) 
p=evalWick(p,ssatape)           # we should not call simplify in this case
nn=gene_nn(N_time_step)
# we could use  a larger step to gc, as we have ~200G 
cal(p,"p",calTree)

for i in 1:length(nn)
    cal(p*nn[i],"nn_$(i)",calTree;num_step_to_gc=20000)
    [GC.gc() for _ in 1:2]      # ensure garbage collection
end

#  as we have larger memory now, we could reduce the gc times
for i in 1:da.N
    for j in 1:da.N
        for orb in 1:N_orbital
            idx_up=defaultSpatialIndex(orb,1)
            idx_down=defaultSpatialIndex(orb,2)
            n_orb_i_j=0.5*(op(da,"dm",[idx_up,i,idx_up,j])+op(da,"dm",[idx_down,i,idx_down,j]))
            cal(p*n_orb_i_j,"g_$(orb)_$(i)_$(j)",calTree;num_step_to_gc=20000)
        end
    end
    [GC.gc() for _ in 1:2]      # ensure garbage collection
end

# it seems that we should implement as lisp like symbolic engine
# maybe just the Expr, but just use a uniform interface.
# dump(:(1+2))
outputSyms=Symbol.(["p",["nn_$(i)" for i in 1:4]..., ["g_$(orb)_$(i)_$(j)" for orb in 1:N_orbital  for j in 1:da.N for i in 1:da.N]...])

func_no_diff_v1=SSACompiledFunc(ssatape,outputSyms)
saveSSAFunc(func_no_diff_v1,"./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff_v1")
func_no_diff_v1=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff_v1")
func_no_diff=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff")

# test_input=rand(68)             #  16*2+18+18
# sum(abs.(func_no_diff_v1(test_input...)-func_no_diff(test_input...)))

cal_diff_with_G(p,"p",num_step_to_gc=20000)

for i in 1:length(nn)
    cal_diff_with_G(p*nn[i],"nn_$(i)",num_step_to_gc=20000)
    [GC.gc() for _ in 1:2]      # ensure garbage collection
end
for i in 1:da.N
    for j in 1:da.N
        for orb in 1:N_orbital
            idx_up=defaultSpatialIndex(orb,1)
            idx_down=defaultSpatialIndex(orb,2)
            n_orb_i_j=0.5*(op(da,"dm",[idx_up,i,idx_up,j])+op(da,"dm",[idx_down,i,idx_down,j]))
            cal_diff_with_G(p*n_orb_i_j,"g_$(orb)_$(i)_$(j)",num_step_to_gc=20000)
        end
    end
    [GC.gc() for _ in 1:2]      # ensure garbage collection
end

outputSymsDiff=Symbol.([gene_diff_name("p")...,vcat([gene_diff_name("nn_$(i)") for i in 1:4]...)...,  vcat([gene_diff_name("g_$(orb)_$(i)_$(j)") for orb in 1:N_orbital  for j in 1:da.N for i in 1:da.N]...)...])
outputSymsWithDiff=[outputSyms...,outputSymsDiff...]

func_with_diff_v1=SSACompiledFunc(ssatape,outputSymsWithDiff)
saveSSAFunc(func_with_diff_v1,"./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff_v1")
func_with_diff_v1=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff_v1")

