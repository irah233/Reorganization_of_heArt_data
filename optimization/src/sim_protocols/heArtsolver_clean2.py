from dolfin import *
from dolfin_adjoint import *
from ufl import cofac, indices
import numpy as np
import os
from pyadjoint.placeholder import Placeholder
from mpi4py import MPI as pyMPI

def heArtsolver_clean(Tmax_ctrls, Vtargets, Ptargets, mesh, facetboundaries, edgeboundaries, f0, s0, n0, SimDet, IODet, u_expt=None, isannotate=True, w_init=None):

        comm = mesh.mpi_comm()
        rank = comm.Get_rank()

        parameters["form_compiler"]["representation"]="uflacs"
        parameters["form_compiler"]["quadrature_degree"]=4

        Velem = VectorElement("CG", mesh.ufl_cell(), 2, quad_scheme="default")
        Qelem = FiniteElement("CG", mesh.ufl_cell(), 1, quad_scheme="default")
        Relem = FiniteElement("Real", mesh.ufl_cell(), 0, quad_scheme="default")
        VRelem = MixedElement([Relem, Relem, Relem, Relem, Relem])
        Rspace = FunctionSpace(mesh, "R", 0)

        W = FunctionSpace(mesh, MixedElement([Velem, Qelem, Relem, VRelem]))
        w = Function(W, annotate=isannotate)
        wtest = TestFunction(W)
        dw = TrialFunction(W)

        (u, p, pendo, c) = split(w)
        (u_, p_, pendo_, c_) = w.split(annotate=isannotate, deepcopy=True)
        (v, q, qendo, cq) = TestFunctions(W)

        X = SpatialCoordinate(mesh)
        N = FacetNormal (mesh)
        d = len(u)
        I = Identity(d)             # Identity tensor
        Fmat = I + grad(u)             # Deformation gradient
        Cmat = Fmat.T*Fmat                   # Right Cauchy-Green tensor

        JJ = det(Fmat)
        Emat = 0.5*(Cmat - I)
        dx_ = dx(mesh)

        # Passive Material parameters
        bff = SimDet["GiccioneParams"]["Passive params"]["bff"]
        bfx = SimDet["GiccioneParams"]["Passive params"]["bfx"]
        bxx = SimDet["GiccioneParams"]["Passive params"]["bxx"]
        Cparam = SimDet["GiccioneParams"]["Passive params"]["Cparam"]

        Eff = inner(f0, Emat*f0)
        Ess = inner(s0, Emat*s0)
        Enn = inner(n0, Emat*n0)
        Efs = inner(f0, Emat*s0)
        Efn = inner(f0, Emat*n0)
        Ens = inner(n0, Emat*s0)
        Esf = inner(s0, Emat*f0)
        Enf = inner(n0, Emat*f0)
        Esn = inner(s0, Emat*n0)
        
        QQ = bff*Eff**2.0 + bxx*(Ess**2.0 + Enn**2.0 + Ens**2.0 +  Esn**2.0) + bfx*(Efs**2.0 + Esf**2.0 + Efn**2.0 + Enf**2.0)
        
        Wp = Cparam/2.0*(exp(QQ) -  1.0) - 1.0*p*(JJ - 1.0)
       
        # Volume constraints
        dsendo = ds(SimDet["LVendoid"], domain = mesh, subdomain_data = facetboundaries)
        area = assemble(Constant(1.0) * dsendo, annotate=isannotate)#, form_compiler_parameters={"representation":"uflacs"})
        V0 = Expression(("vol"), vol=0.0, degree=0, annotate=isannotate)
        V0.vol = Vtargets[0]
        x = u + X
        n = cofac(Fmat)*N
        
        V_u = -Constant(1.0/3.0) * inner(x, n)
        Wvol = (Constant(1.0)/area * pendo  * V0 * dsendo) - (pendo * V_u *dsendo)
 
        # Include active stress
        i, j = indices(2)
        Mij = f0[i]*f0[j]
        Tmax = Function(Rspace, annotate=isannotate)
        Ca0max = SimDet["GiccioneParams"]["Active params"]["Ca0max"]
        Ca0 = SimDet["GiccioneParams"]["Active params"]["Ca0"]
        B = SimDet["GiccioneParams"]["Active params"]["B"]
        lr = SimDet["GiccioneParams"]["Active params"]["lr"]
        l0 = SimDet["GiccioneParams"]["Active params"]["l0"]
        lmbda = sqrt(dot(f0, Cmat*f0))
        ls = lmbda*lr
        ls_l0 = conditional(le(ls, l0+.002), 0.002, ls - l0)
        denom = sqrt(exp(B*(ls_l0)) - 1)
        ECa = Ca0max/denom
        #Pactive = Tmax*as_tensor(Mij, (i,j))
        Pactive = (Tmax*Ca0**2.0)/(Ca0**2.0 + ECa**2.0)*as_tensor(Mij, (i,j))
        F4 = inner(Fmat*Pactive, grad(v))*dx_

        # Rigid BC
        Wrigid = inner(as_vector([c[0], c[1], 0.0]), u) + \
             inner(as_vector([0.0, 0.0, c[2]]), cross(X, u)) + \
             inner(as_vector([c[3], 0.0, 0.0]), cross(X, u)) + \
             inner(as_vector([0.0, c[4], 0.0]), cross(X, u))
        F5 = derivative(Wrigid, w, wtest)*dx_

        # Weak form
        F1 = derivative(Wp, w, wtest)*dx_
        F2 = derivative(Wvol, w, wtest)
        Ftotal = F1 + F2 + F4 + F5

        Jac1 = derivative(F1, w, dw) 
        Jac2 = derivative(F2, w, dw) 
        Jac4 = derivative(F4, w, dw) 
        Jac5 = derivative(F5, w, dw)
        Jac = Jac1 + Jac2 + Jac4 + Jac5

        topid = SimDet["topid"] 
        bctop = DirichletBC(W.sub(0).sub(2), Expression(("0.0"), degree = 2), facetboundaries, topid)
        bcs = [bctop]
 
        problem = NonlinearVariationalProblem(Ftotal, w, bcs=bcs, J=Jac)
        nsolver = NonlinearVariationalSolver(problem)
        nsolver.parameters['nonlinear_solver'] = 'newton'                                                                   
        nsolver.parameters["newton_solver"]["maximum_iterations"] = 50
        nsolver.parameters["newton_solver"]["absolute_tolerance"] = 1e-9
        nsolver.parameters["newton_solver"]["relative_tolerance"] = 1e-8

        if(w_init is not None):
            w.assign(w_init, annotate=isannotate)

        j = 0
        pp = Placeholder(w)
        w0 = Function(W)
        w0.vector().get_local()[:] = 0
        pp.set_value(w0)

        # Volume definition
        vol_form = -Constant(1.0/3.0) * inner(det(Fmat)*dot(inv(Fmat).T, N), X + u)*dsendo

        LVVarray = []
        LVParray = []

        # Preloading
        nload = 10; EDV = Vtargets[0]; unV = 62;
        for p in range(0, nload):
            Vtarget = unV + (EDV - unV)/(nload - 1)*p
            V0.vol = Vtarget
            nsolver.solve(annotate=isannotate)
            (u_, p_, pendo_, c_) = w.split(annotate=isannotate, deepcopy=True)
            LVV = assemble(vol_form, annotate=isannotate)
            LVP = GetLVP(pendo_)
#
            LVVarray.append(LVV)
            LVParray.append(LVP)

            if(rank == 0):
                print("LVP = ", LVP, " LVV = ", LVV, flush=True)


        j = 0
        for p in range(0,len(Ptargets)):

            Tmax.assign(Tmax_ctrls[p], annotate=isannotate)
            V0.vol = Vtargets[p]
            nsolver.solve(annotate=isannotate)
            (u_, p_, pendo_, c_) = w.split(annotate=isannotate, deepcopy=True)
            Ptarget_ = Ptargets[p]
            j += assemble((pendo_*0.0075 - Constant(Ptarget_))**2*dx_, annotate=isannotate)/assemble(Constant(1.0)*dx_, annotate=isannotate)
            LVV = assemble(vol_form, annotate=isannotate)
            LVP = GetLVP(pendo_)

            LVVarray.append(LVV)
            LVParray.append(LVP)

            Tmax_ = GetTmax(Tmax_ctrls[p])

            if(rank == 0):
                print(" Tmax = ", Tmax_, " LVP = ", LVP, " LVV = ", LVV, flush=True)

            #Ptarget += 2


        return w, j, LVParray, LVVarray

def GetLVP(pendo):

    comm = pendo.function_space().mesh().mpi_comm()

    try:
        val_local = pendo.vector().get_local()[0]*0.0075
    except IndexError:
        val_local = 0.0
    
    
    pressure = MPI.sum(comm, val_local)
    
    return pressure

def GetTmax(Tmax_ctrl):

    comm = Tmax_ctrl.function_space().mesh().mpi_comm()

    try:
        val_local = Tmax_ctrl.vector().get_local()[0]
    except IndexError:
        val_local = 0.0
    
    
    Tmax = MPI.sum(comm, val_local)
    
    return Tmax




