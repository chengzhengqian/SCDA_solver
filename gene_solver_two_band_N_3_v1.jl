# this is for run the cluster

using Wick
using SymEngine
using CZQUtils

include("./gene_solver_two_band_common.jl")

"""
proj, we use a different way to parameterize local projector
using u1,..,u18
# here, we just have three time step, but how only have P1 and P2, as P3=1
"""
function proj_full(t)
    if(t<=2)
        p_t=gene_P_full(t)
        return simplify(sum([symbols("u$(i)")*p_t[i] for i in 1:18]))
    end
end


N_orbital,N_spin_orbital,N_time_step,da,input,input_args,total_input_args,ssatape,uniopMap,calTree=init_para(3)


# generate ops                    

p=proj_full(1)*proj_full(2)
p=evalWick(p,ssatape)
nn=gene_nn(N_time_step)

# compute outputs
cal(p,"p",calTree)
for i in 1:length(nn)
    cal(p*nn[i],"nn_$(i)",calTree)
end
for i in 1:da.N
    for j in 1:da.N
        for orb in 1:N_orbital
            idx_up=defaultSpatialIndex(orb,1)
            idx_down=defaultSpatialIndex(orb,2)
            n_orb_i_j=0.5*(op(da,"dm",[idx_up,i,idx_up,j])+op(da,"dm",[idx_down,i,idx_down,j]))
            cal(p*n_orb_i_j,"g_$(orb)_$(i)_$(j)",calTree)
        end
    end
end


outputSyms=Symbol.(["p",["nn_$(i)" for i in 1:4]..., ["g_$(orb)_$(i)_$(j)" for orb in 1:N_orbital  for j in 1:da.N for i in 1:da.N]...])

func_no_diff_v1=SSACompiledFunc(ssatape,outputSyms)
saveSSAFunc(func_no_diff_v1,"./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff_v1")
func_no_diff_v1=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff_v1")
func_no_diff=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff")

test_input=rand(36)
sum(abs.(func_no_diff_v1(test_input...)-func_no_diff(test_input...)))
# check the result

# now, we start compute the derivative againt G
# we just accumulate the new outputs to the previous calTree
cal_diff_with_G(p,"p")
for i in 1:length(nn)
    cal_diff_with_G(p*nn[i],"nn_$(i)")
end
for i in 1:da.N
    for j in 1:da.N
        for orb in 1:N_orbital
            idx_up=defaultSpatialIndex(orb,1)
            idx_down=defaultSpatialIndex(orb,2)
            n_orb_i_j=0.5*(op(da,"dm",[idx_up,i,idx_up,j])+op(da,"dm",[idx_down,i,idx_down,j]))
            cal_diff_with_G(p*n_orb_i_j,"g_$(orb)_$(i)_$(j)")
        end
    end
end


outputSymsDiff=Symbol.([gene_diff_name("p")...,vcat([gene_diff_name("nn_$(i)") for i in 1:4]...)...,  vcat([gene_diff_name("g_$(orb)_$(i)_$(j)") for orb in 1:N_orbital  for j in 1:da.N for i in 1:da.N]...)...])
outputSymsWithDiff=[outputSyms...,outputSymsDiff...]

func_with_diff_v1=SSACompiledFunc(ssatape,outputSymsWithDiff)
saveSSAFunc(func_with_diff_v1,"./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff_v1")
func_with_diff_v1=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff_v1")

func_with_diff=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff")

test_input=rand(36)
@time result_1=func_with_diff(test_input...)
@time result_2=func_with_diff_v1(test_input...)
sum(abs.(result_1-result_2))

# now, we start with third stage, the compoent form
# we first init paras

N_orbital,N_spin_orbital,N_time_step,da,input,input_args,total_input_args,ssatape,uniopMap,calTree=init_para(3)

proj1=gene_P_full(1)
proj2=gene_P_full(2)
N_terms=length(proj1)

# create project matrix
proj_matrix=Array{Any,2}(undef,N_terms,N_terms)
for i in 1:N_terms
    for j in 1:N_terms
        proj_matrix[i,j]=proj1[i]*proj2[j]
    end
end


cal_component(1,"p",proj_matrix)
for n in 1:length(nn)
    cal_component(nn[n],"nn_$(n)",proj_matrix)
end        
for i in 1:da.N
    for j in 1:da.N
        for orb in 1:N_orbital
            idx_up=defaultSpatialIndex(orb,1)
            idx_down=defaultSpatialIndex(orb,2)
            n_orb_i_j=0.5*(op(da,"dm",[idx_up,i,idx_up,j])+op(da,"dm",[idx_down,i,idx_down,j]))
            cal_component(n_orb_i_j,"g_$(orb)_$(i)_$(j)",proj_matrix)
        end
    end
end

output_p=component_symbol(:p)
output_nn=vcat([component_symbol("nn_$(n)") for n in 1:length(nn)]...)
output_g=vcat([component_symbol("g_$(orb)_$(i)_$(j)") for orb in 1:N_orbital for j in 1:da.N for i in 1:da.N]...)
outputSyms=[output_p...,output_nn...,output_g...]

func_component_v1=SSACompiledFunc(ssatape,outputSyms)
saveSSAFunc(func_component_v1,"./gene/two_band_N_$(N_time_step)_spin_symmetric_component_v1")
func_component_v1=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_component_v1")

func_component=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_component")

N_component_size=N_terms^2
test_G=rand(2*da.N^2)
test_u=rand(18)

test_para=[test_G...,test_u...]

result1=func_component_v1(test_para...)
result2=func_component(test_para...)
sum(abs.(result1-result2))











