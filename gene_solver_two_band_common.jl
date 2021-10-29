# this is the common file for the local CDA solver.

"""
compute a given multi-operator and set the output name in ssatape
info contains uniopMap,ssatape
"""
function cal(mop,name,info;num_step_to_gc=5000)
    uniopMap,ssatape=info
    evalTape(evalWick(mop,uniopMap,ssatape;num_step_to_gc=num_step_to_gc),Symbol(name),ssatape)
end

"""
X operator for a given orbital and time step
"""
function x_ops(da,orb,t)
    idx_up=defaultSpatialIndex(orb,1)
    idx_dn=defaultSpatialIndex(orb,2)
    n_up=op(da,"dm",[idx_up,t,idx_up,t])
    n_dn=op(da,"dm",[idx_dn,t,idx_dn,t])
    p_up=1-n_up
    p_dn=1-n_dn
    Xe=p_up*p_dn
    Xup=n_up*p_dn
    Xdn=p_up*n_dn
    Xd=n_up*n_dn
    Xe,Xup,Xdn,Xd
end

"""
a_i_spin1†a_j_spin2
hopping between two spin orbital, 
"""
function hopping(i,spin1,j,spin2,t)
    idx1=defaultSpatialIndex(i,spin1)
    idx2=defaultSpatialIndex(j,spin2)
    op(da,"dm",[idx1,t,idx2,t])
end

"""
generate the projector for the Hamiltonian
as we want to control the density, we don't assume any symmetry about the projector
"""
function gene_P_full(t)
    X1e,X1up,X1dn,X1d=x_ops(da,1,t)
    X2e,X2up,X2dn,X2d=x_ops(da,2,t)
    n1up1dn=hopping(1,1,1,2,t)
    n1dn1up=hopping(1,2,1,1,t)
    n2up2dn=hopping(2,1,2,2,t)
    n2dn2up=hopping(2,2,2,1,t)
    n1up2dn=hopping(1,1,2,2,t)
    n1dn2up=hopping(1,2,2,1,t)
    n2up1dn=hopping(2,1,1,2,t)
    n2dn1up=hopping(2,2,1,1,t)
    p1=X1e*X2e
    p2=X1up*X2e
    p3=X1dn*X2e
    p4=X1e*X2up
    p5=X1e*X2dn
    p6=X1up*X2up
    p8=X1dn*X2dn
    # p7=Basic("0.5")*(X1up*X2dn+X1dn*X2up
    #                  +n1up1dn*n2dn2up+n1dn1up*n2up2dn)
    p7=X1up*X2dn
    p9=X1dn*X2up
    # p9=Basic("0.5")*(X1up*X2dn+X1dn*X2up
    #                  -n1up1dn*n2dn2up-n1dn1up*n2up2dn)
    # p10=Basic("0.5")*(X1d*X2e+X1e*X2d
    #                   +n1up2dn*n1dn2up+n2up1dn*n2dn1up)
    # p11=Basic("0.5")*(X1d*X2e+X1e*X2d
    #                   -n1up2dn*n1dn2up-n2up1dn*n2dn1up)
    p10=X1d*X2e
    p11=X1e*X2d
    p12=X1d*X2up
    p13=X1d*X2dn
    p14=X1up*X2d
    p15=X1dn*X2d
    p16=X1d*X2d
    p17=n1up1dn*n2dn2up+n1dn1up*n2up2dn
    p18=n1up2dn*n1dn2up+n2up1dn*n2dn1up
    [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18]
end

"""
generate the two particle density operator for  Hamiltonian
see the notes for details
"""
function gene_nn(t)
    n1up=hopping(1,1,1,1,t)
    n1dn=hopping(1,2,1,2,t)
    n2up=hopping(2,1,2,1,t)
    n2dn=hopping(2,2,2,2,t)
    n1up1dn=hopping(1,1,1,2,t)
    n1dn1up=hopping(1,2,1,1,t)
    n2up2dn=hopping(2,1,2,2,t)
    n2dn2up=hopping(2,2,2,1,t)
    n1up2dn=hopping(1,1,2,2,t)
    n1dn2up=hopping(1,2,2,1,t)
    n2up1dn=hopping(2,1,1,2,t)
    n2dn1up=hopping(2,2,1,1,t)
    nn_1=n1up*n1dn+n2up*n2dn  # coupling with U
    nn_2=n1up*n2dn+n1dn*n2up  # coupling with U'=U-2J
    nn_3=n1up*n2up+n1dn*n2dn  # coupling with U'-J=U-3J
    nn_4=(n1dn1up*n2up2dn+n2dn2up*n1up1dn
          +n1dn2up*n1up2dn+n2dn1up*n2up1dn) #  coupling with -J
    [nn_1,nn_2,nn_3,nn_4]
end

