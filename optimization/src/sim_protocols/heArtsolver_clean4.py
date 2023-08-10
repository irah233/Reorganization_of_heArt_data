from dolfin import *
from dolfin_adjoint import *
from ufl import cofac, indices
import numpy as np
import os
from pyadjoint.placeholder import Placeholder
from mpi4py import MPI as pyMPI
import math

def heArtsolver_clean(Tmax_ctrls, Vtargets, Ptargets, mesh, facetboundaries, edgeboundaries, f0, s0, n0, eC0, eL0, eR0, matid, SimDet, IODet, isannotate=True, w_init=None, eCCtargets=None, eLLtargets=None ):

        comm = mesh.mpi_comm()
        rank = comm.Get_rank()

        parameters["form_compiler"]["representation"]="uflacs"
        parameters["form_compiler"]["quadrature_degree"]=4

        Velem = VectorElement("CG", mesh.ufl_cell(), 2, quad_scheme="default")
        Qelem = FiniteElement("CG", mesh.ufl_cell(), 1, quad_scheme="default")
        Relem = FiniteElement("Real", mesh.ufl_cell(), 0, quad_scheme="default")
        #VRelem = MixedElement([Relem, Relem, Relem, Relem, Relem])
        VRelem = MixedElement([Relem, Relem, Relem])#, Relem, Relem])
        Telem2 = TensorElement("DG", mesh.ufl_cell(), 0, shape=2*(3,), quad_scheme='default')

        Rspace = Tmax_ctrls[0].function_space()
        TFspace = FunctionSpace(mesh, Telem2)
        DGspace = FunctionSpace(mesh, "DG", 0)
        W = FunctionSpace(mesh, MixedElement([Velem, Qelem, Relem, VRelem]))

        w = Function(W, annotate=isannotate)
        w_prev = Function(W, annotate=isannotate)
        F_ED = Function(TFspace, annotate=isannotate)
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
        dx_ = dx(mesh, subdomain_data=matid)

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
        #Wrigid = inner(as_vector([c[0], c[1], 0.0]), u) + \
             #inner(as_vector([0.0, 0.0, c[2]]), cross(X, u)) + \
             #inner(as_vector([c[3], 0.0, 0.0]), cross(X, u)) + \
             #inner(as_vector([0.0, c[4], 0.0]), cross(X, u))

        Wrigid = inner(as_vector([c[0], c[1], 0.0]), u) + \
             inner(as_vector([0.0, 0.0, c[2]]), cross(X, u))# + \
             #inner(as_vector([c[3], 0.0, 0.0]), cross(X, u)) + \
             #inner(as_vector([0.0, c[4], 0.0]), cross(X, u))
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
        if("abs_tol" in SimDet.keys()):
            nsolver.parameters["newton_solver"]["absolute_tolerance"] = SimDet["abs_tol"]
        else:
            nsolver.parameters["newton_solver"]["absolute_tolerance"] = 1e-9
        if("rel_tol" in SimDet.keys()):
            nsolver.parameters["newton_solver"]["relative_tolerance"] = SimDet["rel_tol"]
        else:
            nsolver.parameters["newton_solver"]["relative_tolerance"] = 1e-8
        if("max_newton_iter" in SimDet.keys()):
            nsolver.parameters["newton_solver"]["maximum_iterations"] = SimDet["max_newton_iter"]
        else:
            nsolver.parameters["newton_solver"]["maximum_iterations"] = 40


        if(w_init is not None):
            w.assign(w_init, annotate=isannotate)

        j = 0
        pp = Placeholder(w)
        w0 = Function(W)
        w0.vector().get_local()[:] = 0
        pp.set_value(w0)

        # Volume definition
        vol_form = -Constant(1.0/3.0) * inner(det(Fmat)*dot(inv(Fmat).T, N), X + u)*dsendo
        init_LVV = assemble(vol_form, annotate=isannotate)
        print("Initial LV volume = ", init_LVV, flush=True)

        LVVarray = []
        LVParray = []
        eCCarray = []
        eLLarray = []
        eRRarray = []

        outdir = IODet["outputfolder"]
        #fdata = File(outdir+"disp.pvd")
        # Preloading
        nload = SimDet["nLoadSteps"]; EDV = Vtargets[0]; unV = init_LVV;
        for p in range(0, nload):
            Vtarget = unV + (EDV - unV)/(nload - 1)*p
            V0.vol = Vtarget
            nsolver.solve(annotate=isannotate)
            (u_, p_, pendo_, c_) = w.split(annotate=isannotate, deepcopy=True)
            u_.rename("disp", "disp")
            fdata << u_
            #File("u_disp"+".pvd") << u_     ### preloading deformation
            LVV = assemble(vol_form, annotate=isannotate)
            LVP = GetLVP(pendo_)

            LVVarray.append(LVV)
            LVParray.append(LVP)

            if(rank == 0):
                print("LVP = ", LVP, " LVV = ", LVV, flush=True)

            w_prev.assign(w, annotate=False)

            
        # Get End-Diastole Deformation Gradient
        F_ED.assign(project(Fmat, TFspace))

        j = 0

        Tmax_prev = Function(Rspace)
        Tmax_temp = Function(Rspace)

        # Starting of PV loop
        for p in range(0,len(Ptargets)):

            Tmax.assign(Tmax_ctrls[p], annotate=isannotate)
            Tmax_avg, Tmax_std, Tmax_arr = GetTmax(Tmax_ctrls[p])

            V0.vol = Vtargets[p]
         
            # Re-run using more cycles:
            ncyc = 1
            if("Nintermediate_stp" in SimDet.keys()):
                ncyc = SimDet["Nintermediate_stp"]

            for cyc in range(0,ncyc):
                Tmax_temp = (cyc+1)/ncyc*(Tmax_ctrls[p] - Tmax_prev) + Tmax_prev
                if(cyc == ncyc -1):
                    Tmax.assign(Tmax_ctrls[p], annotate=True)
                    w.assign(w_prev, annotate=True)
                    nsolver.solve(annotate=True)
                    w_prev.assign(w, annotate=True)

                else:
                    Tmax.assign(Tmax_temp, annotate=True)
                    w.assign(w_prev, annotate=True)
                    nsolver.solve(annotate=True)
                    w_prev.assign(w, annotate=True)


            (u_, p_, pendo_, c_) = w.split(annotate=isannotate, deepcopy=True)
            Ptarget_ = Ptargets[p]

            # Compute minimization functional with Ptarget
            if("Klvp" in SimDet.keys()):
                Klvp = SimDet["Klvp"]
            else:
                Klvp = 1e1
            j += Klvp*assemble(((pendo_*0.0075 - Constant(Ptarget_))/Constant(Ptarget_))**2*dx_, annotate=isannotate)/(assemble(Constant(1.0)*dx_, annotate=isannotate))

            # Add regularization penalty
            if("Kreg" in SimDet.keys()):
                Kreg = SimDet["Kreg"]
            else:
                Kreg = 1e-4
            avg_Tmax = assemble(Tmax_ctrls[p]*dx_, annotate=isannotate)/assemble(Constant(1.0)*dx_, annotate=isannotate)
            j += Kreg*assemble(((Tmax_ctrls[p] - avg_Tmax)/avg_Tmax)**2*dx_, annotate=isannotate)/assemble(Constant(1.0)*dx_, annotate=isannotate)

            # Compute minimization functional with eCCtarget
            if("Kstr" in SimDet.keys()):
                Kstr = SimDet["Kstr"]
            if("Kell" in SimDet.keys()):
                Kell = SimDet["Kell"]
            if("Kecc" in SimDet.keys()):
                Kecc = SimDet["Kecc"]
            else:
                Kstr = 0
                Kell = 1e1
                Kecc = 1e1

            eCC_, eLL_, eRR_ = GetStrainForm(eC0, eR0, eL0, Fmat, DGspace, F_ED, isannotate)
            if not (eCCtargets is None):
                eCCtarget_ = eCCtargets.T[p]
                for k in range(0, len(eCCtarget_)):
                    j += Kecc*assemble((eCC_ - Constant(eCCtarget_[k]))*dx_(k+1), annotate=isannotate)**2/(assemble(Constant(1.0)*dx_(k+1), annotate=isannotate)**2)

            if not (eLLtargets is None):
                eLLtarget_ = eLLtargets.T[p]
                for k in range(0, len(eLLtarget_)):
                    j += Kell*assemble((eLL_ - Constant(eLLtarget_[k]))*dx_(k+1), annotate=isannotate)**2/(assemble(Constant(1.0)*dx_(k+1), annotate=isannotate)**2)

            # Get strain
            eCC, eRR, eLL = GetStrain(eC0, eR0, eL0, Fmat, DGspace, F_ED, matid, dx_, isannotate)


            LVV = assemble(vol_form, annotate=isannotate)
            LVP = GetLVP(pendo_)

            LVVarray.append(LVV)
            LVParray.append(LVP)
            eCCarray.append(eCC)
            eLLarray.append(eLL)
            eRRarray.append(eRR)

            Tmax_prev.assign(Tmax_ctrls[p], annotate=isannotate)
########################### 
            #Tmaxarray1 = assemble(Tmax_arr[p]*dx_, annotate=isannotate)
            Tmaxarray2 = assemble(Tmax_ctrls[p]*dx_, annotate=isannotate)
###########################

            if(rank == 0):
                print(" Tmax_avg = ", Tmax_avg, "+/-", Tmax_std,  " LVP = ", LVP, " LVV = ", LVV, flush=True)

            w_prev.assign(w, annotate=False)

        j = j/len(Ptargets)

        return w, j, LVParray, LVVarray, eCCarray, eLLarray, eRRarray, Tmaxarray2 

def GetLVP(pendo):

    comm = pendo.function_space().mesh().mpi_comm()

    try:
        val_local = pendo.vector().get_local()[0]*0.0075
    except IndexError:
        val_local = 0.0
    
    
    pressure = MPI.sum(comm, val_local)
    
    return pressure

def GetTmax(Tmax_):

    comm = MPI.comm_world
    rank = comm.Get_rank()
    N = Tmax_.vector().size()

    try:
        val_local_sum = np.sum(Tmax_.vector().get_local()[:])
    except IndexError:
        val_local_sum = 0.0

    Tmax_avg = MPI.sum(comm, val_local_sum)/N
    Tmax_arr = Tmax_.vector().get_local()[:]

    try:
        val_local_diff_arr = (Tmax_.vector().get_local()[:] - Tmax_avg)
        val_local_sum_diff_sq = np.sum(np.multiply(val_local_diff_arr, val_local_diff_arr))
    except IndexError:
        val_local_sum_diff_sq = 0.0

    Tmax_std = math.sqrt(MPI.sum(comm, val_local_sum_diff_sq)/N)

    return Tmax_avg, Tmax_std, Tmax_arr

def GetStrainForm(eC0, eR0, eL0, Fmat, Fspace, F_ref, isannotate):

    i, j, k = indices(3)
    F = Fmat*inv(F_ref)
    I = Identity(3)
    Emat = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))

    eCC = 100*project(inner(eC0, Emat*eC0), Fspace, annotate=isannotate)
    eLL = 100*project(inner(eL0, Emat*eL0), Fspace, annotate=isannotate)
    eRR = 100*project(inner(eR0, Emat*eR0), Fspace, annotate=isannotate)

    return eCC, eLL, eRR

def GetStrain(eC0, eR0, eL0, Fmat, Fspace, F_ref, matid, dx_, isannotate):

    max_matid = int(GetMaxVal(matid))
    min_matid = int(GetMinVal(matid))

    eCC, eLL, eRR = GetStrainForm(eC0, eR0, eL0, Fmat, Fspace, F_ref, isannotate)

    eCC_array = [assemble(eCC*dx_(p), annotate=isannotate)/assemble(Constant(1.0)*dx_(p), annotate=isannotate) \
            for p in range(min_matid, max_matid+1)]
    eLL_array = [assemble(eLL*dx_(p), annotate=isannotate)/assemble(Constant(1.0)*dx_(p), annotate=isannotate) \
            for p in range(min_matid, max_matid+1)]
    eRR_array = [assemble(eRR*dx_(p), annotate=isannotate)/assemble(Constant(1.0)*dx_(p), annotate=isannotate) \
            for p in range(min_matid, max_matid+1)]

    return eCC_array, eRR_array, eLL_array

def GetMaxVal(func):

    comm = func.mesh().mpi_comm()
    maxval = comm.allreduce(max(func.array()), op=pyMPI.MAX)

    return (maxval)
 
def GetMinVal(func):

    comm = func.mesh().mpi_comm()
    minval = comm.allreduce(min(func.array()), op=pyMPI.MIN)

    return (minval)
 
