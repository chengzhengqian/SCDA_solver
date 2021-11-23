using SSA
using MathExpr
using Serialization


N_time_step=3
ssatape,input_args,total_input_args,outputSyms=deserialize("./gene/store_info_two_band_N_$(N_time_step)_v3.jls")
diffMap=initDefaultDiffDict(ssatape)

diffNodes,diffSyms=diffTape(outputSyms,input_args,diffMap,ssatape)

outputSymsWithDiff=[outputSyms...,diffSyms...]

func_with_ssadiff_v3=SSACompiledFunc(ssatape,outputSymsWithDiff)
saveSSAFunc(func_with_ssadiff_v3,"./gene/two_band_N_$(N_time_step)_spin_symmetric_with_ssadiff_v3")
func_with_ssadiff_v3=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_ssadiff_v3")

func_with_diff=loadSSAFunc("./gene/two_band_N_$(N_time_step)_spin_symmetric_with_diff")


test_input=rand(36)
@time result_1=func_with_diff(test_input...)
@time result_2=func_with_ssadiff_v3(test_input...)
sum(abs.(result_1-result_2))

# now, they agrees at N=3
# we try run it for N=4 at cluster
