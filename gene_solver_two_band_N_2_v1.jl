# we update the way how we compute the derivative with G, also, we first evaluate the

using Wick
using SymEngine
using CZQUtils

include("./gene_solver_two_band_common.jl")

function proj_full(t)
    p_t=gene_P_full(t)
    simplify(sum([symbols("u$(i)")*p_t[i] for i in 1:18]))
end


# for N=2
N_orbital=2
N_spin_orbital=2*N_orbital
N_time_step=2
da=DA(N_spin_orbital,N_time_step)


input,input_args=initInputMultiBandSpinSymmetric(da)
ssatape=SSATape()
total_input_args=[input_args...,[Symbol("u$(i)") for i in 1:18]...]
setInputArgs(ssatape,total_input_args) #  this is important
uniopMap=initUniOpMap(input,ssatape)
calTree=[uniopMap,ssatape]      # to mimic the prevouis API

p=proj_full(1)*proj_full(2)
nn=gene_nn(N_time_step)

# the difference is that we first evaluate the coefficinet
p=evalWick(p,ssatape)
cal(p,"p",calTree)
# cal(p,"p",calTree; num_step_to_gc=10) 
# control the gc
for i in 1:length(nn)
    cal(p*nn[i],"nn_$(i)",calTree)
end


# we should average over the spin up and down.
# orb=1
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

# we first compute the normal function
func_no_diff_v1=SSACompiledFunc(ssatape,outputSyms)
saveSSAFunc(func_no_diff_v1,"./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff_v1")
func_no_diff_v1=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff_v1")
func_no_diff=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_no_diff")

test_input=rand(26)
func_no_diff_v1(test_input...)-func_no_diff(test_input...)

# now, we add diff using the operator differentiation
# the order is ∂p∂g1_11,∂g1_21,...
# k=1;l=1;orb_=1;
# we organize the code a little bit

function cal_diff_with_G(O,target)
    for k in 1:da.N
        for l in 1:da.N
            for orb_ in 1:N_orbital
                tag="G_$(orb_)_$(k)_$(l)"
                print("$(tag)\n")
                idx_up=defaultSpatialIndex(orb_,1)
                idx_down=defaultSpatialIndex(orb_,2)
                n_orb_k_l_up=op(da,"dm",[idx_up,k,idx_up,l])
                n_orb_k_l_dn=op(da,"dm",[idx_down,k,idx_down,l])
                ∂O∂n_orb_k_l=opDiff(O,n_orb_k_l_up)+opDiff(O,n_orb_k_l_dn)
                cal(∂O∂n_orb_k_l,"∂$(target)∂$(tag)",calTree)
            end
        end
    end
end

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


function gene_diff_name(target)
    [ "∂$(target)∂G_$(orb)_$(i)_$(j)"  for orb in 1:N_orbital  for j in 1:da.N for i in 1:da.N]
end

outputSymsDiff=Symbol.([gene_diff_name("p")...,vcat([gene_diff_name("nn_$(i)") for i in 1:4]...)...,  vcat([gene_diff_name("g_$(orb)_$(i)_$(j)") for orb in 1:N_orbital  for j in 1:da.N for i in 1:da.N]...)...])
outputSymsWithDiff=[outputSyms...,outputSymsDiff...]


func_with_diff_v1=SSACompiledFunc(ssatape,outputSymsWithDiff)
saveSSAFunc(func_with_diff_v1,"./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff_v1")
func_with_diff_v1=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff_v1")

func_with_diff=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff")

test_input=rand(26)
@time result_1=func_with_diff(test_input...)
@time result_2=func_with_diff_v1(test_input...)
sum(abs.(result_1-result_2))
# now, everything agrees

# now, we generate the component form using the new version

N_orbital=2
N_spin_orbital=2*N_orbital
N_time_step=2
da=DA(N_spin_orbital,N_time_step)

input,input_args=initInputMultiBandSpinSymmetric(da)
ssatape=SSATape()
total_input_args=[input_args...,[Symbol("u$(i)") for i in 1:18]...]
setInputArgs(ssatape,total_input_args) #  this is important
uniopMap=initUniOpMap(input,ssatape)
calTree=[uniopMap,ssatape]      # to mimic the prevouis API

proj1=gene_P_full(1)
proj2=gene_P_full(2)
N_terms=length(proj1)
proj_matrix=Array{Any,2}(undef,N_terms,N_terms)
# this holds the matrix of P_i(1)*P_j(2) and P(t)=Σ U_i*P_i(t)
# this has no effect, as it contains on parameters
for i in 1:N_terms
    for j in 1:N_terms
        proj_matrix[i,j]=evalWick(proj1[i]*proj2[j],ssatape)
    end
end

"""
We wrap them into a function
"""
function cal_component(O,name)
    for i in 1:N_terms
        for j in 1:N_terms
            cal(proj_matrix[i,j]*O,"$(name)_$(i)_$(j)",calTree)
        end
    end
end

cal_component(1,"p")
for n in 1:length(nn)
    cal_component(nn[n],"nn_$(n)")
end        

for i in 1:da.N
    for j in 1:da.N
        for orb in 1:N_orbital
            idx_up=defaultSpatialIndex(orb,1)
            idx_down=defaultSpatialIndex(orb,2)
            n_orb_i_j=0.5*(op(da,"dm",[idx_up,i,idx_up,j])+op(da,"dm",[idx_down,i,idx_down,j]))
            cal_component(n_orb_i_j,"g_$(orb)_$(i)_$(j)")
        end
    end
end


function component_symbol(sym)
    [Symbol("$(sym)_$(i)_$(j)") for j in 1:N_terms for i in 1:N_terms]
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
test_G=rand(8)
test_u=rand(18)
test_para=[test_G...,test_u...]

result1=func_component_v1(test_para...)
result2=func_component(test_para...)
sum(abs.(result1-result2))
