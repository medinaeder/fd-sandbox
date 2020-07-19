from firedrake import *

mesh = Mesh("mesh/act.msh")
E = 1.0
nu = 0.4
# Elasticity parameters
mu1, lmbda1 = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))
E0 = 0.1
mu2, lmbda2 = Constant(E0/(2*(1 + nu))), Constant(E0*nu/((1 + nu)*(1 - 2*nu)))


lmbda1 = lmbda1+5./6*mu1
mu1 = 4./3*mu1
lmbda2 = lmbda2+5./6*mu2
mu2 = 4./3*mu2
def PK1_BW(r,mu,lmbda):

    I = Identity(3)             # Identity tensor
    F = variable(I + grad(r))             # Deformation gradient
    C = F.T*F                   # Right Cauchy-Green tensor

    # Invariants of deformation tensors
    Ic = tr(C)
    J  = det(F)

    # Stored strain energy density (compressible neo-Hookean model)
    psi_ = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2 # PEF

    return diff(psi_, F)

def PK1_new(r,mu,lmbda):

    I = Identity(3)             # Identity tensor
    F = variable(I + grad(r))             # Deformation gradient
    C = F.T*F                   # Right Cauchy-Green tensor

    # Invariants of deformation tensors
    Ic = tr(C)
    J  = det(F)

    # Stored strain energy density (compressible neo-Hookean model)
    alpha = 1 + mu/lmbda - mu/(4*lmbda) 
    
    psi_ = (mu/2)*(Ic-3) + lmbda/2*(J-alpha)**2 - mu/2*ln(Ic+1) 

    return diff(psi_, F)


PK1 = PK1_new
V = VectorFunctionSpace(mesh,"CG",2)
v = TestFunction(V)
u = Function(V)
u0 = Function(V)
N = FacetNormal(mesh)
bc = DirichletBC(V, Constant((0,0,0)), 3)

F0 = grad(u0)+Identity(3)
J0 = det(F0)
p0 = Constant(0.0)

R =  inner(PK1(u,mu1,lmbda1),grad(v))*dx(1)+inner(PK1(u,mu2,lmbda2),grad(v))*dx(2)

f = p0*J0*dot(dot(inv(F0.T),N),v)*ds(4)
R+=f



sp  =  {
               "mat_type": "aij",
               "snes_max_it": 10,
               "snes_atol": 1.0e-8,
               "snes_rtol": 1.0e-10,
               "snes_max_linear_solve_fail": 200,
               "snes_linesearch_type": "l2",
               "snes_monitor": None,
               "snes_converged_reason": None,
               "ksp_type": "gmres",
               "ksp_monitor_true_residual": None,
               "ksp_view": None,
               "ksp_monitor_cancel": None,
               "ksp_converged_reason": None,
               "ksp_max_it": 100,
               "pc_type": "hypre",
               "pc_hypre_type": "boomeramg",
               "pc_factor_mat_solver_type": "mumps",
               }

outfile = File("output_new/results.pvd")
u.rename("displacements")
for i in range(100):
    solve(R==0, u, bcs=bc, solver_parameters=sp)
    u0.assign(u)
    p0.assign(2e-4*(i+1))

    outfile.write(u)

