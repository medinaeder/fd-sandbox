import matplotlib.pyplot as plt
from firedrake import *
import numpy as np

mesh = Mesh("mesh/sphere.msh")
E = 1.0 
nu = 0.45
# Elasticity parameters
mu, lmbda = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))
def PK1(r,mu,lmbda):

    I = Identity(3)             # Identity tensor
    F = variable(I + grad(r))             # Deformation gradient
    C = F.T*F                   # Right Cauchy-Green tensor

    # Invariants of deformation tensors
    Ic = tr(C)
    J  = det(F)

    # Stored strain energy density (compressible neo-Hookean model)
    psi_ = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2 # PEF

    return diff(psi_, F)



V = VectorFunctionSpace(mesh,"CG",1)
R = FunctionSpace(mesh, "R", 0)
W = V*R
z = Function(W)
zt  = TestFunction(W)

u,r = split(z)
v,q = split(zt)

X = SpatialCoordinate(mesh)
n = FacetNormal(mesh)
F = Identity(3) + grad(u) 
J = det(F)
trFinv = transpose(inv(F))
scale = Constant(1.0)
dvol_cav  = r*Constant(-1./3)*dot(X+u,trFinv*n)*J*ds(2)-r*scale*Constant(-1./3)*dot(X,Identity(3)*n)*ds(2)
F1 = inner(PK1(u, mu, lmbda), grad(v))*dx 
F2 = derivative(dvol_cav,z,zt)
bc = DirichletBC(W.sub(0), Constant((0,0,0)), 3)
R = F1+F2 
sp = {"mat_type":"matfree",
      "ksp_type":"fgmres",
      "pc_type":"fieldsplit",
      "pc_fieldsplit_type":"schur",
      "pc_fieldsplit_schur_fact_type":"full",
      "pc_fieldsplit_0_fields":"0",
      "pc_fieldsplit_1_fields":"1",
      "fieldsplit_0_ksp_type": "preonly",
      "fieldsplit_0_pc_type": "python",
      "fieldsplit_0_pc_python_type": "firedrake.AssembledPC",
      "fieldsplit_0_assembled_pc_type": "lu",
      "fieldsplit_0_assembled_pc_factor_mat_solver_type": "mumps",
      "fieldsplit_1_ksp_type": "gmres",
      "fieldsplit_1_ksp_monitor_true_residual": None,
      "fieldsplit_1_ksp_max_it": 1,
      "fieldsplit_1_ksp_convergence_test": "skip",
      "fieldsplit_1_pc_type": "none"}

problem = NonlinearVariationalProblem(R, z, bcs=bc)
solver = NonlinearVariationalSolver(problem,solver_parameters = sp)
ss = np.linspace(1.01,2,20)
fil= File("test.pvd")
for s in ss:
    scale.assign(s)
    solver.solve()
    u,r = z.split()
    fil.write(u)
# Delta Volume 
# Drop the coordinate multiply by r
# if I do \int r*Constant(0.01)*dx # increase in off 0.01 volume?
# u,p = split(zz) or is it zz.split()
# derivative(PSI,U,z)--> residual
#print("Total Volume of the cavity", assemble(vol_cav))




import sys; sys.exit()
print(V.dim())
v = TestFunction(V)
N = FacetNormal(mesh)
bc = DirichletBC(V, Constant((0,0,0)), 3)

u0 = Function(V)
F0 = grad(u0)+Identity(3)
J0 = det(F0)
p0 = Constant(0.0)

u = Function(V)
u.rename("displacements")
R =  inner(PK1(u,mu,lmbda),grad(v))*dx
R+= p0*J0*dot(dot(inv(F0.T),N),v)*ds(2)

u1 = Function(V)
u1.rename("displacements_2")
R1 = inner(PK1(u1,mu,lmbda),grad(v))*dx
R1+= p0*dot(N,v)*ds(2)
sp  =  {
               "mat_type": "aij",
               "snes_max_it": 10,
               "snes_atol": 1.0e-8,
               "snes_rtol": 1.0e-10,
               "snes_max_linear_solve_fail": 200,
               "snes_linesearch_type": "l2",
               #"snes_monitor": None,
               "snes_converged_reason": None,
               "ksp_type": "preonly",
               #"ksp_monitor_true_residual": None,
               #"ksp_monitor_cancel": None,
               #"ksp_converged_reason": None,
               "ksp_max_it": 100,
               "pc_type": "lu",
               #"pc_hypre_type": "boomeramg",
               "pc_factor_mat_solver_type": "mumps",
               }




outfile = File("output/results_p2.pvd")
diff = []

er = Function(V)

# Maybe set up a nonlinear variational problem
# Then perform a solve and assign
problem1 = NonlinearVariationalProblem(R, u, bcs=bc)
solver1 = NonlinearVariationalSolver(problem1,solver_parameters = sp)
problem2 = NonlinearVariationalProblem(R1, u1, bcs=bc)
solver2 = NonlinearVariationalSolver(problem2,solver_parameters = sp)

for i in range(100):
    print(i)
    solver1.solve()
    solver2.solve()
    #solve(R==0, u, bcs=bc, solver_parameters=sp)
    #solve(R1==0, u1, bcs=bc, solver_parameters=sp)
    u0.assign(u)
    p0.assign(1e-3*(i+1))
    er.assign(u-u1)
    
    diff.append(assemble(dot(u-u1,u-u1)*dx)/assemble(Constant(1)*dx(domain=mesh)))
    outfile.write(u,u1,er)
    

plt.plot(diff)
plt.show()

