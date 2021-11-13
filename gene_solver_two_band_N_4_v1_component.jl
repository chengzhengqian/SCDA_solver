# we generate the component form in this file, the  no diff and with diff is in the file ....v1.jl

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

# for component form, we have two function, P2, and P1, we first treat the coponent for P2xP2, in this case, we only need the variational parameter for P1, but for consistency, we formally list all them, like what we have done in N=3
N_orbital,N_spin_orbital,N_time_step,da,input,input_args,total_input_args,ssatape,uniopMap,calTree=init_para(4)

proj1=gene_P_full(1)
proj2=gene_P_full(2)
proj3=gene_P_full(3)
proj4=gene_P_full(4)


N_terms=length(proj1)
# P2xP2 (P2xP3) part
proj_matrix=Array{Any,2}(undef,N_terms,N_terms)
for i in 1:N_terms
    for j in 1:N_terms
        proj_matrix[i,j]=proj2[i]*proj3[j]
    end
end

p_remain=proj_full(1)*proj_full(4) # for N=3, p_remain is 1
nn=gene_nn(N_time_step)

cal_component(p_remain,"p",proj_matrix)

#  the remain code need to be run


for n in 1:length(nn)
    cal_component(p_remain*nn[n],"nn_$(n)",proj_matrix;num_step_to_gc=20000)
end        

for i in 1:da.N
    for j in 1:da.N
        for orb in 1:N_orbital
            idx_up=defaultSpatialIndex(orb,1)
            idx_down=defaultSpatialIndex(orb,2)
            n_orb_i_j=0.5*(op(da,"dm",[idx_up,i,idx_up,j])+op(da,"dm",[idx_down,i,idx_down,j]))
            cal_component(p_remain*n_orb_i_j,"g_$(orb)_$(i)_$(j)",proj_matrix;num_step_to_gc=20000)
        end
    end
end

#  for P2xP3 part
output_p=component_symbol(:p)
output_nn=vcat([component_symbol("nn_$(n)") for n in 1:length(nn)]...)
output_g=vcat([component_symbol("g_$(orb)_$(i)_$(j)") for orb in 1:N_orbital for j in 1:da.N for i in 1:da.N]...)
outputSyms=[output_p...,output_nn...,output_g...]

func_component_t_2_3_v1=SSACompiledFunc(ssatape,outputSyms)
saveSSAFunc(func_component_t_2_3_v1,"./gene/two_band_N_$(N_time_step)_spin_symmetric_component_t_2_3_v1")
func_component_t_2_3_v1=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_component_t_2_3_v1")

test_input=rand(68)
@time func_component_t_2_3_v1(test_input...)



# now, we start with 1,4 component.


N_terms=length(proj1)
# P1xP1 (P1xP4) part
proj_matrix=Array{Any,2}(undef,N_terms,N_terms)
for i in 1:N_terms
    for j in 1:N_terms
        proj_matrix[i,j]=proj1[i]*proj4[j]
    end
end

p_remain=proj_full(2)*proj_full(3) # for N=3, p_remain is 1
nn=gene_nn(N_time_step)

N_orbital,N_spin_orbital,N_time_step,da,input,input_args,total_input_args,ssatape,uniopMap,calTree=init_para(4)
cal_component(p_remain,"p",proj_matrix)
#  the remain code need to be run
for n in 1:length(nn)
    cal_component(p_remain*nn[n],"nn_$(n)",proj_matrix;num_step_to_gc=20000)
end        
for i in 1:da.N
    for j in 1:da.N
        for orb in 1:N_orbital
            idx_up=defaultSpatialIndex(orb,1)
            idx_down=defaultSpatialIndex(orb,2)
            n_orb_i_j=0.5*(op(da,"dm",[idx_up,i,idx_up,j])+op(da,"dm",[idx_down,i,idx_down,j]))
            cal_component(p_remain*n_orb_i_j,"g_$(orb)_$(i)_$(j)",proj_matrix;num_step_to_gc=20000)
        end
    end
end

#  for P1xP4 part
output_p=component_symbol(:p)
output_nn=vcat([component_symbol("nn_$(n)") for n in 1:length(nn)]...)
output_g=vcat([component_symbol("g_$(orb)_$(i)_$(j)") for orb in 1:N_orbital for j in 1:da.N for i in 1:da.N]...)
outputSyms=[output_p...,output_nn...,output_g...]

func_component_t_1_4_v1=SSACompiledFunc(ssatape,outputSyms)
saveSSAFunc(func_component_t_1_4_v1,"./gene/two_band_N_$(N_time_step)_spin_symmetric_component_t_1_4_v1")
func_component_t_1_4_v1=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_component_t_1_4_v1")
