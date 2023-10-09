def run_pwa_sim(pDens, dFile, wFile, num_diags, z_offd = 0, z_offw = 0):
    import numpy as np
    from scipy.constants import c, e, m_e


# Import the relevant structures in FBPIC
    from fbpic.main import Simulation
    from fbpic.lpa_utils.bunch import add_particle_bunch_file
    from fbpic.lpa_utils.bunch import add_particle_bunch_openPMD
    from fbpic.lpa_utils.laser.laser_profiles import GaussianLaser
    from fbpic.openpmd_diag import FieldDiagnostic, ParticleDiagnostic, set_periodic_checkpoint, restart_from_checkpoint

# ----------
# Parameters
# ----------

# Whether to use the GPU
    use_cuda = True

# Order of the stencil for z derivatives in the Maxwell solver.
# Use -1 for infinite order, i.e. for exact dispersion relation in
# all direction (adviced for single-GPU/single-CPU simulation).
# Use a positive number (and multiple of 2) for a finite-order stencil
# (required for multi-GPU/multi-CPU with MPI). A large `n_order` leads
# to more overhead in MPI communications, but also to a more accurate
# dispersion relation for electromagnetic waves. (Typically,
# `n_order = 32` is a good trade-off.)
# See https://arxiv.org/abs/1611.05712 for more information.
    n_order = 32

# The simulation box
    Nz = 800         # Number of gridpoints along z
    zmax = 110.e-6    # Right end of the simulation box (meters)
    zmin = -300.e-6   # Left end of the simulation box (meters)
    Nr = 100          # Number of gridpoints along r
    rmax = 350.e-6    # Length of the box along r (meters)
    Nm = 2           # Number of modes used

# The simulation timestep
    dt = (zmax-zmin)/Nz/c   # Timestep (seconds)

# The particles
    p_zmin = 120.e-6  # Position of the beginning of the plasma (meters)
    p_zmax = 20.e-3 # Position of the end of the plasma (meters)
    p_rmax = 400.e-6  # Maximal radial position of the plasma (meters)

    #n_e = 1.e16*1.e6 # Density (electrons.meters^-3)
    n_e = pDens


    p_nz = 2         # Number of particles per cell along z
    p_nr = 2         # Number of particles per cell along r
    p_nt = 4         # Number of particles per cell along theta

# The bunch
    sig_r = 30e-6
    sig_z = 30e-6
    n_emit = 1.e-6
    gamma0 = 19569
    sig_gamma = 20
    n_physical_particles = 1000e-12/e
    n_macroparticles = 10**6

    z0_beam = 0.0
    zf_b = 0.0

    tf_b = np.abs(z0_beam - zf_b)/c

# The moving window
    v_window = c       # Speed of the window

# The diagnostics and the checkpoints/restarts
    
    save_checkpoints = False # Whether to write checkpoint files
    checkpoint_period = 500  # Period for writing the checkpoints
    use_restart = False      # Whether to restart from a previous checkpoint
    track_electrons = False  # Whether to track and write particle ids

# The density profile

    Plasma_start = 100e-6
    Plasma_end   = 20.e-3
    RampGradient = 0.333e-3
    Upramp_center = 2.25e-3
    Downramp_center = 16.75e-3

    def dens_func( z, r ):
   
        n = np.ones_like(z)
    # Make plasma shape
        n = np.where( z<Plasma_end, 1/(1+np.exp(-(z-Upramp_center)/RampGradient) )/(1+np.exp((z-Downramp_center)/RampGradient)), n )
    # Supress density before the ramp
        n = np.where( z<Plasma_start, 0, n )
    # Supress density after the ramp
        n = np.where( z>Plasma_end, 0, n )
  
        return(n)

# The interaction length of the simulation (meters)
    L_interact = 21.e-3 # increase to simulate longer distance!
# Interaction time (seconds) (to calculate number of PIC iterations)
    T_interact = ( L_interact + (zmax-zmin) ) / v_window
# (i.e. the time it takes for the moving window to slide across the plasma)

# ---------------------------
# Carrying out the simulation
# ---------------------------

# NB: The code below is only executed when running the script,
# (`python lwfa_script.py`), but not when importing it (`import lwfa_script`).


    # Initialize the simulation object
    sim = Simulation( Nz, zmax, Nr, rmax, Nm, dt, zmin=zmin,
        n_order=n_order, use_cuda=use_cuda,
        boundaries={'z':'open', 'r':'open'})
        # 'r': 'open' can also be used, but is more computationally expensive

    # Create the plasma electrons
    elec = sim.add_new_species( q=-e, m=m_e, n=n_e,
        dens_func=dens_func, p_zmin=p_zmin, p_zmax=p_zmax, p_rmax=p_rmax,
        p_nz=p_nz, p_nr=p_nr, p_nt=p_nt )

    # Load initial fields
    # Create a Gaussian laser profile

    with open(dFile) as file:
        macroNd = len(file.readlines())

    with open(wFile) as file:
        macroNw = len(file.readlines())

    
    chargeMagd = 2*10**-15
    weightd = chargeMagd/e
    massd = weightd*m_e

    
    chargeMagw = 2*10**-15
    weightw = chargeMagw/e
    massw = weightw*m_e

    dbunch = add_particle_bunch_file(sim, -e, m_e, dFile, macroNd*weightd, z_off = z_offd) 
    wbunch = add_particle_bunch_file(sim, -e, m_e, wFile, macroNw*weightw,z_off= z_offw)

    
    #bunch = add_particle_bunch_openPMD(sim, -e, m_e, eFile)


    if use_restart is False:
        # Track electrons if required (species 0 correspond to the electrons)
        if track_electrons:
            elec.track( sim.comm )
            dbunch.track( sim.comm)
            wbunch.track( sim.comm)
    else:
        # Load the fields and particles from the latest checkpoint file
        restart_from_checkpoint( sim )

    # Configure the moving window
    sim.set_moving_window( v=v_window )

    # Number of iterations to perform
    N_step = int(T_interact/sim.dt)

    diag_period = 200
    #diag_period = N_step/(num_diags-1)-1         # Period of the diagnostics in number of timesteps

    # Add diagnostics

    sim.diags = [ FieldDiagnostic( diag_period, sim.fld, comm=sim.comm ),
                  ParticleDiagnostic( diag_period, {"electrons_plasma" : elec, "electrons_drive" : dbunch, "electrons_witness" : wbunch},
                    select={"uz" : [1., None ]}, comm=sim.comm ) ]

    # Add checkpoints
    if save_checkpoints:
        set_periodic_checkpoint( sim, checkpoint_period )

    ### Run the simulation
    sim.step( N_step )

    
if __name__ == '__main__':
    #fl = 'fbNewPz.txt' #works

    import sys

    if len(sys.argv) > 1:
        dr = sys.argv[1]
        wt = sys.argv[2]
        pDens = sys.argv[3]
    else:
        dr = '/home/mstobbe/FBPIC_Sims/drive_test.txt'
        wt = '/home/mstobbe/FBPIC_Sims/witness_test.txt'

    run_pwa_sim(pDens, dr,wt,30,0,0)

