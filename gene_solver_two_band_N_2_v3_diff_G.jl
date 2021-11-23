# store_info=[ssatape,input_args,total_input_args,outputSyms]
# use one without xxx_diff_G.jl to generate the ssatape and save it.
# here, we just to the new diffTape API to generate the diff expression
# we should load the definition before deserialize
# using Wick
using SSA
using MathExpr
using Serialization
# include("./gene_solver_two_band_common.jl")

N_time_step=2
ssatape,input_args,total_input_args,outputSyms=deserialize("./gene/store_info_two_band_N_$(N_time_step)_v3.jls")
diffMap=initDefaultDiffDict(ssatape)

diffNodes,diffSyms=diffTape(outputSyms,input_args,diffMap,ssatape)

outputSymsWithDiff=[outputSyms...,diffSyms...]

func_with_ssadiff_v3=SSACompiledFunc(ssatape,outputSymsWithDiff)
saveSSAFunc(func_with_ssadiff_v3,"./gene/two_band_N_$(N_time_step)_spin_symmetric_with_ssadiff_v3")
func_with_ssadiff_v3=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_ssadiff_v3")

func_with_diff=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff")

test_input=rand(26)
@time result_1=func_with_diff(test_input...)
@time result_2=func_with_ssadiff_v3(test_input...)
sum(abs.(result_1-result_2))

# now, they agrees

