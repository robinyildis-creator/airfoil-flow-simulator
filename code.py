    import streamlit as st
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.ndimage as ndimage

    # ==========================================
    # 1. LBM COMPUTATIONAL ENGINE
    # ==========================================
    class LBMEngine:
        def __init__(self, nx, ny, tau, v_in):
            self.nx, self.ny = nx, ny
            self.tau = tau
            self.v_in = v_in
            
            # D2Q9 Constants
            self.c = np.array([[0,0], [1,0], [0,1], [-1,0], [0,-1], [1,1], [-1,1], [-1,-1], [1,-1]])
            self.w = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])
            self.noslip = [0, 3, 4, 1, 2, 7, 8, 5, 6]
            
            # Initialization
            self.rho = np.ones((nx, ny))
            self.u = np.zeros((nx, ny, 2))
            self.u[:,:,0] = v_in
            self.f = self.get_equilibrium(self.rho, self.u)

        def get_equilibrium(self, rho, u):
            # Identifiera formen på indata (kan vara (nx, ny) eller bara (ny,))
            shape = rho.shape
            
            # Beräkna cu och usq baserat på indata
            cu = 3.0 * np.dot(u, self.c.T)
            usq = 1.5 * np.sum(u**2, axis=-1)
            
            # Skapa feq med matchande dimensioner + 9 för hastighetsvektorerna
            feq = np.zeros((*shape, 9))
            
            for i in range(9):
                # Vi använder "..." (Ellipsis) för att hantera 1D eller 2D automatiskt
                feq[..., i] = self.w[i] * rho * (1 + cu[..., i] + 0.5 * cu[..., i]**2 - usq)
            return feq

        def step(self, obstacle_mask):
            # Streaming step
            for i in range(9):
                self.f[:,:,i] = np.roll(self.f[:,:,i], self.c[i], axis=(0,1))
            
            # Collision step & Boundary conditions
            f_old = self.f.copy()
            f_boundary = self.f[obstacle_mask, :]
            self.f[obstacle_mask, :] = f_boundary[:, self.noslip]
            
            # Update macro variables
            self.rho = np.sum(self.f, axis=2)
            self.u = np.dot(self.f, self.c) / self.rho[:,:,np.newaxis]
            
            # Inlet and Outlet (Drichlet/Neumann approximation)
            self.u[0, :, 0] = self.v_in
            self.f[0, :, :] = self.get_equilibrium(self.rho[0, :], self.u[0, :, :])
            self.f[-1, :, :] = self.f[-2, :, :]
            
            # Final collision
            feq = self.get_equilibrium(self.rho, self.u)
            self.f += -(1/self.tau) * (self.f - feq)
            
            # Force Calculation (Lift & Drag)
            fx, fy = 0, 0
            for i in range(1, 9):
                diff = (f_old[obstacle_mask, i] + self.f[obstacle_mask, self.noslip[i]])
                fx += np.sum(diff * self.c[i,0])
                fy += np.sum(diff * self.c[i,1])
            return fx, fy

    # ==========================================
    # 2. AERODYNAMICS & GEOMETRY
    # ==========================================
    def naca_0012(x):
        """NACA 0012 Thickness distribution formula"""
        return 0.6 * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)

    def get_wing_mask(nx, ny, chord, aoa_deg):
        x_points = np.linspace(0, 1, int(chord))
        y_points = naca_0012(x_points)
        
        # Create coordinate pairs
        wing_coords = []
        for xi, yi in zip(x_points, y_points):
            wing_coords.append([xi * chord, yi * chord])
            wing_coords.append([xi * chord, -yi * chord])
        
        # Rotation Matrix for Angle of Attack
        theta = np.radians(-aoa_deg)
        rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        rotated_wing = np.array(wing_coords) @ rot.T
        
        # Map to Grid
        mask = np.zeros((nx, ny), dtype=bool)
        cx, cy = nx//4, ny//2
        for px, py in rotated_wing:
            ix, iy = int(px + cx), int(py + cy)
            if 0 <= ix < nx and 0 <= iy < ny:
                mask[ix, iy] = True
        return ndimage.binary_fill_holes(mask)

    # ==========================================
    # 3. INTERACTIVE DASHBOARD
    # ==========================================
    st.set_page_config(page_title="Rob Aero", layout="wide")
    st.title("🦅 Aero: NACA Airfoil Wind Tunnel")

    with st.sidebar:
        st.header("Simulation Parameters")
        aoa = st.slider("Angle of Attack (Degrees)", -15, 20, 5)
        viscosity = st.slider("Viscosity (Tau)", 0.51, 1.0, 0.6)
        v_in = st.slider("Flow Velocity", 0.02, 0.15, 0.1)
        chord = 60

    # Initialize Engine
    if 'engine' not in st.session_state:
        st.session_state.engine = LBMEngine(400, 100, viscosity, v_in)

    engine = st.session_state.engine
    engine.tau = viscosity
    engine.v_in = v_in
    mask = get_wing_mask(400, 100, chord, aoa)

    # Visualization
    col1, col2 = st.columns([3, 1])
    plot_placeholder = col1.empty()
    metric_placeholder = col2.empty()

    for t in range(1000):
        fx, fy = engine.step(mask)
        
        if t % 25 == 0:
            # Calculate Vorticity for Visuals
            vorticity = (np.gradient(engine.u[:,:,1], axis=0) - 
                        np.gradient(engine.u[:,:,0], axis=1))
            
            fig, ax = plt.subplots(figsize=(12, 4))
            # Background: Vorticity field
            ax.imshow(vorticity.T, cmap='RdBu_r', origin='lower', vmin=-0.02, vmax=0.02, alpha=0.9)
            # Foreground: Streamlines (Rök)
            Y, X = np.mgrid[0:100, 0:400]
            ax.streamplot(X, Y, engine.u[:,:,0].T, engine.u[:,:,1].T, color='black', linewidth=0.5, density=1.5)
            ax.contour(mask.T, colors='black', levels=[0.5]) # Draw wing
            ax.set_axis_off()
            plot_placeholder.pyplot(fig)
            plt.close()
            
            # Calculate Engineering Metrics
            cl = fy / (0.5 * v_in**2 * chord)
            cd = fx / (0.5 * v_in**2 * chord)
            
            with metric_placeholder:
                st.metric("Lift Coefficient (Cl)", f"{cl:.3f}")
                st.metric("Drag Coefficient (Cd)", f"{cd:.3f}")
                st.metric("Efficiency (L/D)", f"{abs(cl/cd):.2f}" if cd != 0 else "0")
